import geopandas as gpd
import laspy
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from shapely.geometry import Point
import pandas as pd
import alphashape
from scipy.spatial import Delaunay

def extract_alpha_shape(polygon_points, alpha):
    points_2d = [(p[0], p[1]) for p in polygon_points]
    alpha_shape = alphashape.alphashape(points_2d, alpha)

    return alpha_shape

def calculate_mesh_surface_area(points):
    tri = Delaunay(points[:, :2])
    triangles = tri.simplices

    total_area = 0.0
    for tri_indices in triangles:
        p1, p2, p3 = points[tri_indices]
        area = 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1))
        total_area += area

    return total_area

# wczytanie danych
polygons = gpd.read_file("1939.geojson").to_crs(epsg=2180)
las = laspy.read("1939.las")
x, y, z = las.xyz[:, 0], las.xyz[:, 1], las.xyz[:, 2]
pred_class = las['pred_class']
pred_ID = las['pred_ID']

# filtracja punktów hałd
heap_mask = pred_class == 1
heap_points = np.vstack((x[heap_mask], y[heap_mask], z[heap_mask])).T
pred_ID = pred_ID[heap_mask]

# utworzenie powierzchni bazowej terenu
cell_size = 0.1
cell_area = cell_size ** 2
xmin, xmax = x.min(), x.max()
ymin, ymax = y.min(), y.max()
nx = int((xmax - xmin) / cell_size)
ny = int((ymax - ymin) / cell_size)

grid_x = np.linspace(xmin, xmax, nx)
grid_y = np.linspace(ymin, ymax, ny)
grid_X, grid_Y = np.meshgrid(grid_x, grid_y)
surface = RegularGridInterpolator(
    (grid_y, grid_x),
    griddata((x, y), z, (grid_X, grid_Y), method='linear'),
    bounds_error=False,
    fill_value=0.0
)

# obliczanie statystyk
results = []
alpha_shapes = []

for idx, poly in polygons.iterrows():
    pid = poly["pred_ID"]
    print(f"Przetwarzanie poligonu {pid}...")

    # selekcja punktów z danego poligonu po pred_ID
    mask = pred_ID == pid
    selected_points = heap_points[mask]
    if selected_points.shape[0] == 0:
        continue

    # obliczenie indeksów komórek siatki dla punktów
    i = np.floor((selected_points[:, 0] - xmin) / cell_size).astype(int)
    j = np.floor((selected_points[:, 1] - ymin) / cell_size).astype(int)

    # dataframe do grupowania punktów po komórkach
    df = pd.DataFrame({
        'x': selected_points[:, 0],
        'y': selected_points[:, 1],
        'z': selected_points[:, 2],
        'i': i,
        'j': j
    })

    # grupowanie punktów po komórkach i pobieranie średnich współrzędnych
    sampled = df.groupby(['i', 'j'], as_index=False).mean()

    # interpolacja wysokości terenu
    interp_pts = np.column_stack((sampled['y'], sampled['x']))
    base_heights = surface(interp_pts)
    base_heights = np.nan_to_num(base_heights, nan=0.0)

    # obliczanie różnicy wysokości
    height_diff = sampled['z'].values - base_heights
    height_diff = np.clip(height_diff, 0, None)

    # obliczenie objętości
    volume = np.sum(height_diff) * cell_area

    # obliczenie powierzchni pokrycia poligonu punktami -> zliczenie powierzchni komórek siatki terenu, w obrębie których zawierają się punkty
    unique_cells = df[['i', 'j']].drop_duplicates()
    x_centers = xmin + (unique_cells['i'].values + 0.5) * cell_size
    y_centers = ymin + (unique_cells['j'].values + 0.5) * cell_size
    centers = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x_centers, y_centers), crs=polygons.crs)

    within_mask = centers.within(poly.geometry)
    count_tiles_inside = within_mask.sum()

    coverage_area = count_tiles_inside * cell_area

    # obliczenie powierzchni 3D
    surface_area = calculate_mesh_surface_area(selected_points)

    # wyodrębnianie obrysu hałdy
    alpha_shape = extract_alpha_shape(selected_points, 0.5)

    if isinstance(alpha_shape, Point):
        continue
    elif alpha_shape.geom_type == 'Polygon':
        alpha_shapes.append({
            "pred_ID": pid,
            "geometry": alpha_shape
        })

    results.append({
        "polygon_id": pid,
        "volume_m3": round(volume, 3),
        "coverage_area_m2": round(coverage_area, 3),
        "surface_area_3D_m2": round(surface_area, 3),
    })


# eskport wyników do pliku CSV
result_file = "statystyki_hald.csv"
df = pd.DataFrame(results)
df.to_csv(result_file, index=False)
print(f"Wyniki zapisane do pliku {result_file}")

# eksport obrysów hałd do pliku shapefile
obrysy_file = "obrysy_hald.shp" 
alpha_shape_gdf = gpd.GeoDataFrame(alpha_shapes, geometry="geometry", crs=polygons.crs)
alpha_shape_gdf.to_file(obrysy_file)
print(f"Obrysy zapisane do pliku {obrysy_file}")