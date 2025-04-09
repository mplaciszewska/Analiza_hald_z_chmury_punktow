# Analiza hałd na podstawie danych LiDAR

## Dane wejściowe

- **1939.las** – chmura punktów w układzie ETRF2000-PL/CS92 (EPSG:2180), zawierająca dodatkowe pola klasyfikacji:
  - `pred_class` – klasyfikacja punktów (`0` – brak hałdy, `1` – hałda),
  - `pred_ID` – identyfikator przypisany do konkretnej hałdy.
- **1939.geojson** – obrysy hałd w formacie GeoJSON, w układzie WGS84 (EPSG:4326).

## Opis

Skrypt przetwarza dane z chmury punktów, tworzy na jej podstawie powierzchnię bazową terenu, a następnie dla każdej hałdy (`pred_ID`) oblicza:

- **Objętość hałdy** - względem interpolowanej powierzchni bazowej,
- **Powierzchnię pokrycia poligonu punktami** – czyli obszar zajmowany przez punkty hałdy w siatce,
- **Powierzchnię 3D hałdy** – wyznaczaną na podstawie triangulacji Delaunaya.

Dodatkowo wykonywane jest **automatyczne wyodrębnienie obrysu hałdy** z punktów przy użyciu algorytmu **Alpha Shape**.

## Wyniki

Skrypt generuje dwa pliki wyjściowe:

- `statystyki_hald.csv` – tabela ze statystykami dla każdej hałdy:
  - `polygon_id` – identyfikator hałdy,
  - `volume_m3` – objętość hałdy (m³),
  - `coverage_area_m2` – powierzchnia pokrycia (m²),
  - `surface_area_3D_m2` – powierzchnia 3D hałdy (m²).

- `obrysy_hald.shp` – plik Shapefile zawierający wygenerowane obrysy hałd.

## Uruchomienie

1. Wymagane biblioteki:

```bash
pip install geopandas laspy numpy scipy shapely pandas alphashape
```

2. Uruchomienie skryptu
```bash
python zadanie.py --geojson <ścieżka_do_pliku.geojson> --las <ścieżka_do_pliku.las>
``` 
