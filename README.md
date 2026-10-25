# planetx-bodies-terrain

Height tiles of other bodies for [PlanetX Lite](https://github.com/Valery35/planetx-lite)
and the PlanetX plugin for QGIS.

The workflow `.github/workflows/tiles.yml` runs `scripts/tiles.py` once
by hand. It downloads a global height grid of the body from USGS
Astrogeology, checks it against a known summit, cuts it into tiles and
pushes them to `main` in parts.

| Path | Source | Terms |
|---|---|---|
| `mars/{z}/{x}/{y}.png`, levels 0-7 | Mars MGS MOLA DEM 463 m, NASA, USGS (`Mars_MGS_MOLA_DEM_mosaic_global_463m.tif`), heights above the areoid | Public domain; credit NASA MGS MOLA, USGS |
| `moon/{z}/{x}/{y}.png`, levels 0-7 | Moon LRO LOLA DEM 118 m, NASA, USGS (`Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif`), heights above the sphere of 1737.4 km | Public domain; credit NASA LRO LOLA, USGS |

The tiles are XYZ Web Mercator of the body sphere, 256 × 256 PNG in the
Terrarium encoding of the Mapzen terrain tiles: height in metres =
R × 256 + G + B / 256 - 32768.
