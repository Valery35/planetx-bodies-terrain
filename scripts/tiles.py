# PlanetX Lite - глобус PlanetX без QGIS.
# Copyright (C) 2026 ООО «Информ++». Лицензия Apache 2.0.
# SPDX-License-Identifier: Apache-2.0
"""Тайлы высот Terrarium тела из глобальной сетки высот.

python3 tiles.py сетка.tif папка уровень_от уровень_до

Сетка - глобальный растр в простой цилиндрической проекции (MOLA,
LOLA USGS), охват 360° по долготе и 180° по широте. Высоты - значение
пикселя, умноженное на масштаб полосы GDAL, плюс смещение. Тайлы
256 × 256 в Web Mercator тела, разметка XYZ, папка/{z}/{x}/{y}.png,
кодировка Terrarium: (R · 256 + G + B / 256) - 32768 метров, как
тайлы высот Земли Mapzen. Значение пикселя тайла - билинейно
из сетки, прочитанной с прореживанием GDAL (среднее) под шаг тайла.
"""
import math
import os
import sys

import numpy as np
from osgeo import gdal
from PIL import Image

SIZE = 256


def terrarium(h):
    v = np.clip(h + 32768.0, 0.0, 65535.996)
    r = np.floor(v / 256.0)
    g = np.floor(v - r * 256.0)
    b = np.floor((v - r * 256.0 - g) * 256.0)
    return np.dstack([r, g, b]).astype(np.uint8)


def tile_lats(z, y):
    n = 2 ** z
    t = (y + (np.arange(SIZE) + 0.5) / SIZE) / n
    return np.degrees(np.arctan(np.sinh(math.pi * (1.0 - 2.0 * t))))


def tile_lons(z, x):
    n = 2 ** z
    return (x + (np.arange(SIZE) + 0.5) / SIZE) / n * 360.0 - 180.0


class Grid:
    def __init__(self, path):
        self.ds = gdal.Open(path)
        self.band = self.ds.GetRasterBand(1)
        self.w, self.h = self.ds.RasterXSize, self.ds.RasterYSize
        gt = self.ds.GetGeoTransform()
        # Глобальная сетка: ширина - 360°, высота - 180°. Начало долгот
        # по геопривязке, в метрах или в градусах.
        deg_x = 360.0 / self.w
        deg_y = 180.0 / self.h
        self.lon0 = gt[0] / gt[1] * deg_x
        self.lat0 = gt[3] / -gt[5] * deg_y
        self.deg_x, self.deg_y = deg_x, deg_y
        self.scale = self.band.GetScale() or 1.0
        self.offset = self.band.GetOffset() or 0.0
        self.nodata = self.band.GetNoDataValue()
        print("grid %d x %d, lon0 %.4f, lat0 %.4f, scale %g, offset %g, nodata %s"
              % (self.w, self.h, self.lon0, self.lat0, self.scale,
                 self.offset, self.nodata), flush=True)

    def strip(self, lat_hi, lat_lo, step):
        """Полоса сетки от широты lat_hi до lat_lo, прорежённая в step
        раз: высоты, широта первой строки и шаг в градусах."""
        r0 = max(0, int(math.floor((self.lat0 - lat_hi) / self.deg_y)) - 2)
        r1 = min(self.h, int(math.ceil((self.lat0 - lat_lo) / self.deg_y)) + 2)
        rows = r1 - r0
        bw = max(1, self.w // step)
        bh = max(1, rows // step)
        a = self.band.ReadAsArray(0, r0, self.w, rows, buf_xsize=bw, buf_ysize=bh,
                                  resample_alg=gdal.GRIORA_Average).astype(np.float64)
        if self.nodata is not None:
            a[a == self.nodata] = np.nan
        a = a * self.scale + self.offset
        a[~np.isfinite(a)] = 0.0
        dy = rows * self.deg_y / bh
        dx = self.w * self.deg_x / bw
        top = self.lat0 - r0 * self.deg_y
        return a, top, dy, dx

    def sample(self, strip, lats, lons):
        a, top, dy, dx = strip
        bh, bw = a.shape
        fy = (top - lats) / dy - 0.5
        fx = ((lons - self.lon0) % 360.0) / dx - 0.5
        y0 = np.clip(np.floor(fy).astype(int), 0, bh - 1)
        y1 = np.clip(y0 + 1, 0, bh - 1)
        ty = np.clip(fy - np.floor(fy), 0.0, 1.0)
        x0 = np.floor(fx).astype(int) % bw
        x1 = (x0 + 1) % bw
        tx = fx - np.floor(fx)
        a00, a01 = a[y0][:, x0], a[y0][:, x1]
        a10, a11 = a[y1][:, x0], a[y1][:, x1]
        ty = ty[:, None]
        return (a00 * (1 - tx) + a01 * tx) * (1 - ty) + (a10 * (1 - tx) + a11 * tx) * ty


def build(path, out, z_from, z_to):
    grid = Grid(path)
    for z in range(z_from, z_to + 1):
        n = 2 ** z
        # Прореживание: пиксель сетки не мельче половины пикселя тайла.
        step = max(1, int(grid.w / (SIZE * n * 2)))
        high, low = -1e9, 1e9
        for y in range(n):
            lats = tile_lats(z, y)
            strip = grid.strip(lats[0], lats[-1], step)
            for x in range(n):
                h = grid.sample(strip, lats, tile_lons(z, x))
                high, low = max(high, float(h.max())), min(low, float(h.min()))
                folder = os.path.join(out, str(z), str(x))
                os.makedirs(folder, exist_ok=True)
                Image.fromarray(terrarium(h), "RGB").save(
                    os.path.join(folder, "%d.png" % y), optimize=True)
        print("level %d: %d tiles, heights %.0f..%.0f m" % (z, n * n, low, high),
              flush=True)


def check(path, lat, lon, low, high, name):
    """Высота известной точки в пределах low..high, иначе ошибка."""
    grid = Grid(path)
    strip = grid.strip(lat + 1.0, lat - 1.0, 1)
    h = float(grid.sample(strip, np.array([lat]), np.array([lon]))[0, 0])
    print("%s: %.0f m" % (name, h), flush=True)
    if not low <= h <= high:
        sys.exit("%s: %.0f m вне %g..%g" % (name, h, low, high))


if __name__ == "__main__":
    if sys.argv[1] == "check":
        check(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]),
              float(sys.argv[5]), float(sys.argv[6]), sys.argv[7])
    else:
        build(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
