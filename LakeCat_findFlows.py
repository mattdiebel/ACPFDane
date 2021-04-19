# LatCat functions for defining watershed topology
# Inputs:
# zone_file is the watershed raster
# fdr_file is the flow direction raster

import os
import numpy as np
import pandas as pd
import re
import arcpy
import tempfile


def rollArray(a, d):
    if len(d) == 4:
        a = a[0, :]
        new = np.roll(np.roll(a, d[0], axis=0), d[1], axis=1)
        new[d[2], :] = a[d[2], :]
        new[:, d[3]] = a[:, d[3]]
    if len(d) == 3:
        new = np.roll(a[0, :], d[0], axis=d[1])
        if d[1] == 0:
            new[d[2], :] = a[0, d[2], :]
        if d[1] == 1:
            new[:, d[2]] = a[0, :, d[2]]
    return np.expand_dims(new, axis=0)


def makeFlows(arr, shiftd, fdr, path, nd):
    # cells change value after shift * cells not equal to NoData
    iso = np.not_equal(arr, shiftd) * np.not_equal(shiftd, nd) * np.not_equal(arr, nd)
    pth = np.equal(fdr, path)  # True when equal to path value
    val = iso * pth * arr
    shiftval = iso * pth * shiftd
    idx = np.not_equal(val, shiftd)
    fromcom = val[idx]
    tocom = shiftval[idx]
    fromcom = fromcom[fromcom > 0]
    tocom = tocom[tocom > 0]
    # don't load-in the entire array to the DF, just connection vals
    df = pd.DataFrame({'TOCOMID': tocom,
                       'FROMCOMID': fromcom,
                       'move': path})
    return df.drop_duplicates(['FROMCOMID', 'TOCOMID'])


def compAll(arr, fdr, moves, from_to, nd):
    for move in moves:
        flow = makeFlows(arr, rollArray(np.copy(arr), moves[move][0]), fdr, moves[move][1], nd)
        from_to = pd.concat([from_to, flow])
    return from_to


def expand(window, size=1):
    r, c = window
    return ((r[0] - size, r[1] + size), (c[0] - size, c[1] + size))


def check_window(window, w, h):
    r, c = window
    return ((max(0, r[0]), min(h, r[1])), (max(0, c[0]), min(w, c[1])))


def lower_left_coord(r, window):
    xmin = r.extent.XMin
    ymax = r.extent.YMin
    cell_size = r.meanCellHeight
    lower_x = xmin + window[1][0] * cell_size
    lower_y = ymax - window[0][1] * cell_size
    return arcpy.Point(lower_x, lower_y)


def chunk_windows(r, max_ram=250000000):
    nbytes = int(re.findall(r'\d+', r.pixelType)[0])
    pixel_size = nbytes / 8
    chunk_size, _ = divmod(max_ram, pixel_size)
    r_h, r_w = r.height, r.width
    if chunk_size >= r_h * r_w:
        yield (0, 0), ((0, r_h), (0, r_w))
    else:
        b_h, b_w = [128, 128]
        d, _ = map(int, divmod(chunk_size, r_w * b_h))
        chunk_height = d * b_h
        d, m = map(int, divmod(r_h, chunk_height))
        n_chunks = d + int(m > 0)
        for i in range(n_chunks):
            row = i * chunk_height
            # height = min(chunk_height, r_h - row)
            yield (i, 0), ((row, row + chunk_height), (0, r_w))


def findFlows(zone_file, fdr_file):
    moves = {'up': [(-1, 0, -1), 4], 'left': [(-1, 1, -1), 1], 'down': [(1, 0, 0), 64],
             'right': [(1, 1, 0), 16], 'downRight': [(1, 1, 0, 0), 32],
             'downLeft': [(1, -1, 0, -1), 128], 'upRight': [(-1, 1, -1, 0), 8],
             'upLeft': [(-1, -1, -1, -1), 2]}
    flows = pd.DataFrame()
    with tempfile.TemporaryDirectory() as td:
        temp_z = arcpy.CopyRaster_management(zone_file, os.path.join(td, "z.tif"))
        temp_f = arcpy.CopyRaster_management(fdr_file, os.path.join(td, "fdr.tif"))
        z = arcpy.Raster(temp_z)
        f = arcpy.Raster(temp_f)
        for _, w in chunk_windows(z):  # currently defaults to 250MB
            nd = int(z.noDataValue)
            new_w = check_window(expand(w, 2), z.width, z.height)
            ll = lower_left_coord(z, window=new_w)
            ncols = new_w[0][1] - new_w[0][0]
            nrows = new_w[1][1] - new_w[1][0]
            data = arcpy.RasterToNumPyArray(z, lower_left_corner=ll, ncols=ncols, nrows=nrows)
            data = data.reshape((1, nrows, ncols))
            f_r = arcpy.RasterToNumPyArray(f, lower_left_corner=ll, ncols=ncols, nrows=nrows)
            f_r = f_r.reshape((1, nrows, ncols))
            flows = pd.concat([flows, compAll(data, f_r, moves, flows, nd)])
        del z, f
    return flows.drop_duplicates(['FROMCOMID', 'TOCOMID'])
