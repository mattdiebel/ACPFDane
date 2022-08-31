# Import modules and scripts
import arcpy
from arcpy import env
from arcpy.sa import *
from copy import copy
import sys, string, os, time, re
import numpy as np
import pandas as pd
arcpy.CheckOutExtension("Spatial")
aprx = arcpy.mp.ArcGISProject('CURRENT')
m = aprx.activeMap

# Set environments
arcpy.env.extent = os.path.join(path, "dem")
arcpy.env.snapRaster = os.path.join(path, "dem")
arcpy.env.overwriteOutput = True

# Condition DEM with cutLines and stormsewers
def conditionDEM():
    cut_zmin = arcpy.sa.ZonalStatistics("cutLines", "OBJECTID", "dem", "MINIMUM", "DATA")
    stormsewer_zmin = arcpy.sa.ZonalStatistics("stormsewers", "OBJECTID", "dem", "MINIMUM", "DATA")
    comb_zmin = arcpy.sa.CellStatistics([cut_zmin, stormsewer_zmin], "MINIMUM", "DATA")
    DEMCut = Con(IsNull(comb_zmin), "dem", comb_zmin)
    DEMCut.save(os.path.join(path, "demCut"))
    DEMFill = Fill(DEMCut)
    DEMFill.save(os.path.join(path, "demFill"))
    arcpy.env.parallelProcessingFactor = "100%"
    D8FlowDir = FlowDirection(DEMFill, "", "")    
    D8Accumulation = FlowAccumulation(D8FlowDir, "", "INTEGER")    
    arcpy.env.parallelProcessingFactor = ""
    D8FlowDir.save(os.path.join(path,"D8FlowDir"))
    D8Accumulation.save(os.path.join(path,"D8FlowAcc"))
    m.addDataFromPath(os.path.join(path, "D8FlowDir"))
    m.addDataFromPath(os.path.join(path, "D8FlowAcc"))
    m.addDataFromPath(os.path.join(path, "demFill"))
    m.addDataFromPath(os.path.join(path, "demCut"))
    
def FlowPaths(AreaThreshold):
    number_cells = float(AreaThreshold) * 4046 / 4
    FlowNetRas = Con(os.path.join(path, "D8FlowAcc"), 1, "", "VALUE >= %s" % number_cells)
    Link = StreamLink(FlowNetRas, os.path.join(path, "D8FlowDir"))
    StreamToFeature(Link, os.path.join(path, "D8FlowDir"), os.path.join(path, "flowPaths"), "NO_SIMPLIFY")
    flowPathsCopy = arcpy.management.CopyFeatures(os.path.join(path, "flowPaths"))
    arcpy.AddField_management(flowPathsCopy, "grid_code_to", "LONG")
    arcpy.CalculateField_management(flowPathsCopy, "grid_code_to", "!grid_code!") 
    arcpy.JoinField_management(os.path.join(path, "flowPaths"), "to_node", flowPathsCopy, "from_node", ["grid_code_to"])
    m.removeLayer(m.listLayers("flowPaths_CopyFeatures*")[0])

def Depressions(MinimumSize):
    stormsewer_zmin = arcpy.sa.ZonalStatistics("stormsewers", "OBJECTID", "dem", "MINIMUM", "DATA")
    demCutNoSS = Con(IsNull(stormsewer_zmin), os.path.join(path, "demCut"), os.path.join(path, "demFill"))
    FillReg = Minus(os.path.join(path, "demFill"), demCutNoSS)
    AllSinks = Con(FillReg, 1, "", "VALUE > 0.01")
    Depressions = arcpy.RasterToPolygon_conversion(AllSinks, os.path.join(path, "Depressions"), "NO_SIMPLIFY")
    MinSize = float(MinimumSize) * 4046
    arcpy.MakeFeatureLayer_management(os.path.join(path, "Depressions"), "small_sinks", "\"Shape_Area\" < %s" % MinSize)
    arcpy.DeleteFeatures_management("small_sinks")
    arcpy.AddField_management(os.path.join(path, "Depressions"), "Depress_ID")
    arcpy.CalculateField_management(os.path.join(path, "Depressions"), "Depress_ID", "10000 + !OBJECTID!", "PYTHON3")
    gSSURGO_copy = arcpy.CopyRaster_management("gSSURGO", "gSSURGO_copy", "", 0, 0)
    arcpy.AddField_management(gSSURGO_copy,"HydFlt","FLOAT","8","0")
    arcpy.CalculateField_management(gSSURGO_copy,"HydFlt", '!Hydric!', "PYTHON")
    HydricPct = Lookup(gSSURGO_copy, "HydFlt")
    HydricTable = ZonalStatisticsAsTable(os.path.join(path, "Depressions"), "Depress_ID", HydricPct, "in_memory" + "\\HydricTable", "#", "MEAN")
    arcpy.JoinField_management(os.path.join(path, "Depressions"), "Depress_ID", HydricTable, "Depress_ID", ["MEAN"])
    arcpy.AddField_management(os.path.join(path, "Depressions"), "PctHydric", "FLOAT")
    arcpy.CalculateField_management(os.path.join(path, "Depressions"), "PctHydric", "!MEAN!", "PYTHON")
    arcpy.DeleteField_management(os.path.join(path, "Depressions"), ["MEAN"])
    DepthStats = ZonalStatisticsAsTable(os.path.join(path, "Depressions"), "Depress_ID", FillReg, "in_memory" + "\\DepthStats", "", "MAXIMUM")
    arcpy.JoinField_management(os.path.join(path, "Depressions"), "Depress_ID", DepthStats, "Depress_ID", ["MAX"])
    arcpy.AddField_management(os.path.join(path, "Depressions"), "MaxDepthCM")
    arcpy.CalculateField_management(os.path.join(path, "Depressions"), "MaxDepthCM", "!MAX!", "PYTHON")
    arcpy.DeleteField_management(os.path.join(path, "Depressions"), ["gridcode", "MAX", "id"])
    VolStats = ZonalStatisticsAsTable(os.path.join(path, "Depressions"), "Depress_ID", FillReg, "in_memory" + "\\VolStats", "", "SUM")
    arcpy.JoinField_management(os.path.join(path, "Depressions"), "Depress_ID", VolStats, "Depress_ID", ["SUM"])
    arcpy.AddField_management(os.path.join(path, "Depressions"), "VolAcreFt", "FLOAT")
    SumDepthtoAcFt = 1233.48 / 0.01 / 2
    arcpy.CalculateField_management(Depressions, "VolAcreFt", "(!SUM! / %s)" % SumDepthtoAcFt, "PYTHON")
    arcpy.DeleteField_management(Depressions, ["SUM"])
    m.removeLayer(m.listLayers("gSSURGO_copy")[0])
    m.removeLayer(m.listLayers("small_sinks")[0])
    m.removeTable(m.listTables("HydricTable")[0])
    m.removeTable(m.listTables("DepthStats")[0])
    m.removeTable(m.listTables("VolStats")[0])
    
# Select connected flow paths
def findConnected():
    i = 0
    j = 1
    while True:
        arcpy.SelectLayerByLocation_management("flowPaths", "INTERSECT", "flowPaths", "", "ADD_TO_SELECTION")
        j = int(arcpy.GetCount_management("flowPaths")[0])
        if j == i:
            break
        else:
            i = j

# Use flow path watershed to refine HUC12 boundary
def refineHUC12():
    arcpy.SelectLayerByAttribute_management("flowPaths", "SWITCH_SELECTION")
    arcpy.DeleteFeatures_management("flowPaths")
    arcpy.conversion.FeatureToRaster("flowPaths", "grid_code", os.path.join(path, "pathSeeds"), 2)
    arcpy.env.parallelProcessingFactor = "100%"
    pathWatersheds = arcpy.sa.Watershed("D8FlowDir", "pathSeeds", "Value")
    arcpy.env.parallelProcessingFactor = ""
    arcpy.RasterToPolygon_conversion(pathWatersheds, os.path.join(path, "pathWatershedsPoly"), "NO_SIMPLIFY")
    del pathWatersheds
    arcpy.Delete_management(os.path.join(path,"pathSeeds"))

# Delete depressions that are outside of pathWatersheds extent
def pruneDepressions():
    arcpy.SelectLayerByAttribute_management("Depressions", "CLEAR_SELECTION")
    allDepCount = int(arcpy.GetCount_management("Depressions")[0])
    arcpy.management.SelectLayerByLocation("Depressions", "HAVE_THEIR_CENTER_IN", "pathWatershedsPoly", None, "NEW_SELECTION", "INVERT")
    selDepCount = int(arcpy.GetCount_management("Depressions")[0])
    if allDepCount > selDepCount and selDepCount > 0:
        arcpy.DeleteFeatures_management("Depressions")
    arcpy.Delete_management(os.path.join(path,"pathWatershedsPoly"))

# Erase sections of flow paths that intersect depressions
def pruneFlowPaths():
    arcpy.management.CopyFeatures("flowPaths", os.path.join(path, "flowPathsRaw"), None, None, None, None)    
    arcpy.analysis.Erase("flowPathsRaw", "Depressions", os.path.join(path, "flowPathsMP"), None)
    arcpy.management.MultipartToSinglepart("flowPathsMP", os.path.join(path, "flowPaths"))
    arcpy.Delete_management(os.path.join(path,"flowPathsMP"))
    arcpy.management.FeatureVerticesToPoints("flowPaths", os.path.join(path, "flowPathNodes"), "BOTH_ENDS")
    arcpy.analysis.TabulateIntersection("Depressions", "Depress_ID", "flowPathNodes", os.path.join(path, "pathDepIntersect"), "ORIG_FID", None, None, "UNKNOWN")
    arcpy.JoinField_management("flowPaths", "OBJECTID", "pathDepIntersect", "ORIG_FID", ["Depress_ID","PNT_COUNT"])
    arcpy.management.SelectLayerByAttribute("flowPaths", "NEW_SELECTION", "PNT_COUNT = 2", None)
    arcpy.DeleteFeatures_management("flowPaths")
    arcpy.DeleteField_management("flowPaths", ["Depress_ID","PNT_COUNT"])
    arcpy.management.FeatureVerticesToPoints("flowPaths", os.path.join(path, "flowPathNodes"), "START")
    arcpy.management.SelectLayerByLocation("flowPathNodes", "INTERSECT", "Depressions", None, "NEW_SELECTION", "INVERT")
    arcpy.DeleteFeatures_management("flowPathNodes")
    arcpy.analysis.TabulateIntersection("Depressions", "Depress_ID", "flowPathNodes", os.path.join(path, "pathDepIntersect"), "to_node", None, None, "UNKNOWN")
    arcpy.management.SelectLayerByAttribute("pathDepIntersect", "NEW_SELECTION", "PNT_COUNT < 2", None)
    arcpy.management.DeleteRows("pathDepIntersect")
    arcpy.management.JoinField("flowPathNodes", "to_node", "pathDepIntersect", "to_node", "Depress_ID;PNT_COUNT")
    arcpy.management.SelectLayerByAttribute("flowPathNodes", "NEW_SELECTION", "PNT_COUNT IS NULL", None)
    arcpy.DeleteFeatures_management("flowPathNodes")
    arcpy.analysis.SpatialJoin("flowPathNodes", "Depressions", os.path.join(path, "flowPathNodesDepJoin"), "JOIN_ONE_TO_ONE", "KEEP_ALL", 'ORIG_FID "ORIG_FID" true true false 4 Long 0 0,First,#,flowPathsNoDepNodes,ORIG_FID,-1,-1;Depress_ID "Depress_ID" true true false 4 Long 0 0,First,#,flowPathsNoDepNodes,Depress_ID,-1,-1;Depress_ID_1 "Depress_ID_1" true true false 4 Long 0 0,First,#,Depressions,Depress_ID,-1,-1', "INTERSECT", None, None)
    arcpy.management.SelectLayerByAttribute("flowPathNodesDepJoin", "NEW_SELECTION", "Depress_ID <> Depress_ID_1", None)
    arcpy.DeleteFeatures_management("flowPathNodesDepJoin")
    arcpy.JoinField_management("flowPaths", "OBJECTID", "flowPathNodesDepJoin", "ORIG_FID", ["Depress_ID"])
    arcpy.management.SelectLayerByAttribute("flowPaths", "NEW_SELECTION", "Depress_ID IS NOT NULL", None)
    arcpy.sa.ZonalStatisticsAsTable("flowPaths", "grid_code", "D8FlowAcc", arcpy.env.scratchGDB + "\\flowPathAcc", "DATA", "MEAN")
    arcpy.JoinField_management("flowPaths", "grid_code", arcpy.env.scratchGDB + "\\flowPathAcc", "grid_code", ["MEAN"])

# Define topology for all flow paths and depressions that connect to flow paths
def defineTopology():
    arcpy.SelectLayerByAttribute_management("flowPaths", "CLEAR_SELECTION")
    arcpy.DeleteField_management("flowPaths", ["arcid","grid_code","from_node","to_node","grid_code_to","ORIG_FID","Depress_ID","MEAN"])
    arcpy.management.FeatureVerticesToPoints("flowPaths", os.path.join(path, "flowPathNodesUS"), "START")
    arcpy.management.FeatureVerticesToPoints("flowPaths", os.path.join(path, "flowPathNodesDS"), "END")
    arcpy.management.AddXY(os.path.join(path, "flowPathNodesUS"))
    arcpy.management.AddXY(os.path.join(path, "flowPathNodesDS"))
    arcpy.AddField_management("flowPathNodesUS", "LONLAT", "DOUBLE")
    arcpy.AddField_management("flowPathNodesDS", "LONLAT", "DOUBLE")
    arcpy.management.CalculateField("flowPathNodesUS", "LONLAT", "int(!POINT_X!)*10000000 + int(!POINT_Y!)", "PYTHON3", None)
    arcpy.management.CalculateField("flowPathNodesDS", "LONLAT", "int(!POINT_X!)*10000000 + int(!POINT_Y!)", "PYTHON3", None)
    arcpy.JoinField_management("flowPathNodesDS", "LONLAT", "flowPathNodesUS", "LONLAT", ["ORIG_FID"])
    arcpy.AlterField_management("flowPathNodesDS", "ORIG_FID", "FROM_ID", "FROM_ID")
    arcpy.AlterField_management("flowPathNodesDS", "ORIG_FID_1", "TO_ID", "TO_ID")
    arcpy.analysis.SpatialJoin("flowPathNodesDS", "Depressions", os.path.join(path, "flowPathNodesDSDep"), "JOIN_ONE_TO_ONE", "KEEP_ALL", 'FROM_ID "FROM_ID" true true false 4 Long 0 0,First,#,flowPathNodesDS,FROM_ID,-1,-1;TO_ID "TO_ID" true true false 4 Long 0 0,First,#,flowPathNodesDS,TO_ID,-1,-1;Depress_ID "Depress_ID" true true false 4 Long 0 0,First,#,Depressions,Depress_ID,-1,-1', "INTERSECT", None, None)
    arcpy.management.SelectLayerByAttribute("flowPathNodesDSDep", "NEW_SELECTION", "TO_ID IS NULL", None)
    arcpy.management.CalculateField("flowPathNodesDSDep", "TO_ID", "!Depress_ID!", "PYTHON3", None)
    arcpy.SelectLayerByAttribute_management("flowPathNodesDSDep", "CLEAR_SELECTION")
    arcpy.AlterField_management("flowPathNodesUS", "ORIG_FID", "TO_ID", "TO_ID")
    arcpy.analysis.SpatialJoin("flowPathNodesUS", "Depressions", os.path.join(path, "flowPathNodesUSDep"), "JOIN_ONE_TO_ONE", "KEEP_ALL", 'TO_ID "TO_ID" true true false 4 Long 0 0,First,#,flowPathNodesUS,TO_ID,-1,-1;Depress_ID "Depress_ID" true true false 4 Long 0 0,First,#,Depressions,Depress_ID,-1,-1', "INTERSECT", None, None)
    # Join node topology to flow paths and depressions
    arcpy.AddField_management("flowPaths", "FROM_ID", "LONG")
    arcpy.management.CalculateField("flowPaths", "FROM_ID", "!OBJECTID!", "PYTHON3", None)
    arcpy.JoinField_management("flowPaths", "FROM_ID", "flowPathNodesDSDep", "FROM_ID", ["TO_ID"])
    arcpy.JoinField_management("Depressions", "Depress_ID", "flowPathNodesUSDep", "Depress_ID", ["TO_ID"])
    # Delete node layers
    arcpy.Delete_management(os.path.join(path,"flowPathNodesUS"))
    arcpy.Delete_management(os.path.join(path,"flowPathNodesDS"))
    arcpy.Delete_management(os.path.join(path,"flowPathNodesUSDep"))
    arcpy.Delete_management(os.path.join(path,"flowPathNodesDSDep"))
    arcpy.Delete_management(os.path.join(path,"flowPathNodes"))
    arcpy.Delete_management(os.path.join(path,"flowPathNodesDepJoin"))
    arcpy.Delete_management(os.path.join(path,"pathDepIntersect"))
    arcpy.Delete_management(arcpy.env.scratchGDB + "\\flowPathAcc")

# Convert flow paths and depressions to rasters and combine to make seeds
def watershedSeeds():
    arcpy.conversion.FeatureToRaster("flowPaths", "FROM_ID", os.path.join(path, "pathSeeds"), 2)
    arcpy.conversion.FeatureToRaster("Depressions", "Depress_ID", os.path.join(path, "depSeeds"), 2)
    seeds = Con(IsNull("depSeeds"), "pathSeeds", "depSeeds")
    seeds.save(os.path.join(path, "seeds"))
    arcpy.Delete_management(os.path.join(path,"depSeeds"))
    arcpy.Delete_management(os.path.join(path,"pathSeeds"))

# Delineate watersheds for seeds
def watershed():
    arcpy.env.parallelProcessingFactor = "100%"
    watersheds = arcpy.sa.Watershed("D8FlowDir", os.path.join(path, "seeds"), "Value")
    arcpy.env.parallelProcessingFactor = ""
    watersheds.save(os.path.join(path,"watersheds"))
    arcpy.Delete_management(os.path.join(path,"seeds"))
    
# Convert watershed raster to polygons
def watershedPolygons():
    arcpy.RasterToPolygon_conversion(os.path.join(path, "watersheds"), os.path.join(path, "watershedsPoly"), "NO_SIMPLIFY", "Value", "MULTIPLE_OUTER_PART")
    
# LakeCat functions
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
    ymax = r.extent.YMax
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
    temp_z = arcpy.CopyRaster_management(zone_file, os.path.join(arcpy.env.scratchFolder, "z.tif"))
    temp_f = arcpy.CopyRaster_management(fdr_file, os.path.join(arcpy.env.scratchFolder, "fdr.tif"))
    z = arcpy.Raster(temp_z)
    f = arcpy.Raster(temp_f)
    for _, w in chunk_windows(z):  # currently defaults to 250MB
        nd = int(z.noDataValue)
        new_w = check_window(expand(w, 2), z.width, z.height)
        ll = lower_left_coord(z, window=new_w)
        nrows = new_w[0][1] - new_w[0][0]
        ncols = new_w[1][1] - new_w[1][0]
        data = arcpy.RasterToNumPyArray(z, lower_left_corner=ll, ncols=ncols, nrows=nrows)
        data = data.reshape((1, nrows, ncols))
        f_r = arcpy.RasterToNumPyArray(f, lower_left_corner=ll, ncols=ncols, nrows=nrows)
        f_r = f_r.reshape((1, nrows, ncols))
        flows = pd.concat([flows, compAll(data, f_r, moves, flows, nd)])
    del z, f
    return flows.drop_duplicates(['FROMCOMID', 'TOCOMID'])

# Run watershed topology script from LakeCat
def LakeCat():
    df = findFlows(os.path.join(path, "watersheds"), os.path.join(path, "D8FlowDir"))
    x = np.array(np.rec.fromrecords(df.values))
    names = df.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, arcpy.env.scratchGDB + "\\Topology")

# Compile and clean up watershed attributes
def watershedAttributes():
    # Join depression attributes to watersheds
    arcpy.AddField_management("watershedsPoly", "FROM_ID", "LONG")
    arcpy.CalculateField_management("watershedsPoly", "FROM_ID", "!gridcode!") 
    arcpy.JoinField_management("watershedsPoly", "gridcode", "Depressions", "Depress_ID", ["TO_ID"])
    arcpy.JoinField_management("watershedsPoly", "gridcode", "Depressions", "Depress_ID", ["PctHydric","MaxDepthCM","VolAcreFt"])
    arcpy.AddField_management("watershedsPoly", "AreaAcres", "FLOAT")
    arcpy.CalculateField_management("watershedsPoly", "AreaAcres", "!shape.area@acres!") 
    # Join node-based topology to watersheds
    arcpy.JoinField_management("watershedsPoly", "FROM_ID", "flowPaths", "FROM_ID", ["TO_ID"])
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "NEW_SELECTION", "TO_ID IS NULL And TO_ID_1 IS NOT NULL", None)
    arcpy.CalculateField_management("watershedsPoly", "TO_ID", "!TO_ID_1!") 
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "CLEAR_SELECTION")
    arcpy.DeleteField_management("watershedsPoly", ["TO_ID_1"])
    # Join LakeCat topology to watersheds
    arcpy.JoinField_management("watershedsPoly", "gridcode", arcpy.env.scratchGDB + "\\Topology", "FROMCOMID", ["TOCOMID"])
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "NEW_SELECTION", "TO_ID IS NULL And TOCOMID IS NOT NULL", None)
    arcpy.CalculateField_management("watershedsPoly", "TO_ID", "!TOCOMID!") 
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "CLEAR_SELECTION")
    arcpy.DeleteField_management("watershedsPoly", ["Id","gridcode","TOCOMID"])
    # Replace path-based TO_IDs that have no watershed with next downstream TO_ID
    arcpy.management.CopyRows("watershedsPoly", arcpy.env.scratchGDB + "\\watersheds", None)
    arcpy.JoinField_management("watershedsPoly", "TO_ID", arcpy.env.scratchGDB + "\\watersheds", "FROM_ID", ["FROM_ID"])
    arcpy.JoinField_management("watershedsPoly", "TO_ID", "flowPaths", "FROM_ID", ["FROM_ID","TO_ID"])
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "NEW_SELECTION", "FROM_ID_1 IS NULL And TO_ID IS NOT NULL", None)
    arcpy.CalculateField_management("watershedsPoly", "TO_ID", "!TO_ID_1!") 
    arcpy.SelectLayerByAttribute_management("watershedsPoly", "CLEAR_SELECTION")
    # Clean up
    arcpy.DeleteField_management("watershedsPoly", ["FROM_ID_1","FROM_ID_12","TO_ID_1"])
    arcpy.Delete_management(arcpy.env.scratchGDB + "\\Topology")
    arcpy.Delete_management(arcpy.env.scratchGDB + "\\watersheds")
    
# Calculate runoff volume from a single rainfall event
def runoffEvent(P): 
    # Convert watershedsPoly table to pandas dataframe
    arr = arcpy.da.FeatureClassToNumPyArray(os.path.join(path,"watershedsPoly"), ("FROM_ID", "TO_ID", "VolAcreFt", "AreaAcres", "CNlow", "DSorder"), null_value = 0)
    df2 = pd.DataFrame(arr)
    df2 = df2.sort_values(by = 'DSorder', ascending = False)
    df2 = df2.reset_index(drop=True)

    df2['S'] = 0 # 6
    df2['Qraw'] = 0 # 7
    df2['runDep'] = 0  # 8
    df2['runVolS'] = 0 # 9
    df2['runInS'] = 0 # 10
    df2['runOffS'] = 0 # 11
    df2['runInNDS'] = 0 # 12
    df2['RORUS'] = 0 # 13
    df2['RORInc'] = 0 # 14
    df2['RORDS'] = 0 # 15

    df2.S = 1000 / df2.CNlow - 10
    df2.Qraw = (P - 0.2 * df2.S)**2 / (P + 0.8 * df2.S)
    df2.loc[P < 0.2 * df2.S, 'runDep'] = 0
    df2.loc[P >= 0.2 * df2.S, 'runDep'] = df2.Qraw / 12

    # Loop through watersheds to calculate runIn and runOff
    i = 0
    while i < len(df2):
        us = df2[df2.TO_ID == df2.FROM_ID[i]]
        df2.iloc[i,9] = df2.AreaAcres[i] * df2.runDep[i]
        df2.iloc[i,10] = df2.AreaAcres[i] * df2.runDep[i] + sum(us.runOffS)
        df2.iloc[i,11] = max(df2.runInS[i] - df2.VolAcreFt[i], 0)
        df2.iloc[i,12] = df2.runVolS[i] + sum(us.runInNDS)
        i += 1

    df2.RORUS = df2.runOffS / df2.runInNDS
    df2.loc[df2.RORUS.isnull(), 'RORUS'] = 0
        
    # Loop through watersheds to calculate downstream runoff ratio
    df2.loc[df2.runInS == 0, 'RORInc'] = 0
    df2.loc[df2.runInS > 0, 'RORInc'] = df2.runOffS / df2.runInS
    i = 0
    while i < len(df2):
        if df2.runOffS[i] == 0:
            df2.iloc[i,15] = 0
            i += 1
        else:
            ds = [df2.RORInc[i]]
            to = int(df2.TO_ID[i])
            while True:
                RORto = df2.RORInc[df2.FROM_ID == to].values.tolist()
                ds.extend(RORto)
                if to == 0:
                    break
                to = int(df2.TO_ID[df2.FROM_ID == to])
            df2.iloc[i,15] = np.prod(ds)
            i += 1

    # Convert pandas dataframe to table
    x = np.array(np.rec.fromrecords(df2.values))
    names = df2.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, arcpy.env.scratchGDB + "\\runoff")

    # Join runoff table to watershedsPoly
    arcpy.JoinField_management("watershedsPoly", "FROM_ID", arcpy.env.scratchGDB + "\\runoff", "FROM_ID", ["runVolS","runInS","runOffS","RORUS","RORInc","RORDS"])
    arcpy.Delete_management(arcpy.env.scratchGDB + "\\runoff")

# Make transects across flowPaths for stream power and travel time calculations
def makeTransects():
    # Smoothing doesn't seem to fix the transect orientation issue
    # Identify transects that intersect cutLines for later exclusion?
    # arcpy.management.CopyFeatures("flowPaths", os.path.join(path,"flowPathsSmooth"), None, None, None, None)
    # arcpy.edit.Generalize("flowPathsSmooth", "6 Meters")
    arcpy.DeleteField_management(os.path.join(path,"flowPaths"), ["FlowAcc","RORUS","LengthM","n","USP","TT"])
    arcpy.DeleteField_management(os.path.join(path,"watershedsPoly"), ["TT","TTDS"])
    arcpy.management.GenerateTransectsAlongLines(os.path.join(path,"flowPaths"), os.path.join(path,"transects"), "100 Meters", "32 Meters", "END_POINTS")
    arcpy.AlterField_management(os.path.join(path,"transects"), "ORIG_FID", "pathID", "pathID")
    arcpy.management.GeneratePointsAlongLines(os.path.join(path,"transects"), os.path.join(path,"transectPoints2"), "DISTANCE", "2 Meters", None, "END_POINTS")
    arcpy.AlterField_management(os.path.join(path,"transectPoints2"), "ORIG_FID", "transectID", "transectID")
    arcpy.AddField_management(os.path.join(path,"transectPoints2"), "pointID", "LONG")
    arcpy.CalculateField_management(os.path.join(path,"transectPoints2"), "pointID", "!OBJECTID!")
    arcpy.DeleteField_management(os.path.join(path,"transectPoints2"), ["Shape_Length"]) # Delete Shape_Length field so extract values to points will work
    arcpy.sa.ExtractValuesToPoints(os.path.join(path,"transectPoints2"), os.path.join(path,"demFill"), os.path.join(path,"transectPoints"), "NONE", "VALUE_ONLY")
    arcpy.AlterField_management(os.path.join(path,"transectPoints"), "RASTERVALU", "Elevation", "Elevation")
    arcpy.conversion.FeatureToRaster(os.path.join(path,"flowPaths"), "OBJECTID", os.path.join(path,"flowPathsR"), 2)
    arcpy.sa.ZonalStatisticsAsTable(os.path.join(path,"flowPathsR"), "VALUE", os.path.join(path,"D8FlowAcc"), os.path.join(path,"flowPathsFA"), "DATA", "MEDIAN")
    arcpy.JoinField_management(os.path.join(path,"flowPaths"), "OBJECTID", os.path.join(path,"flowPathsFA"), "Value", ["MEDIAN"])
    arcpy.AlterField_management(os.path.join(path,"flowPaths"), "MEDIAN", "FlowAcc", "FlowAcc")
    arcpy.JoinField_management(os.path.join(path,"flowPaths"), "FROM_ID", os.path.join(path,"watershedsPoly"), "FROM_ID", ["RORUS"])
    arcpy.AddField_management(os.path.join(path,"flowPaths"), "LengthM", "DOUBLE")
    arcpy.CalculateGeometryAttributes_management(os.path.join(path,"flowPaths"), [["LengthM","LENGTH"]], "METERS")
    arcpy.Delete_management(os.path.join(path,"flowPathsFA"))
    arcpy.Delete_management(os.path.join(path,"flowPathsR"))
    arcpy.Delete_management(os.path.join(path,"transectPoints2"))
    
# Join results of R Manning test script to flowPaths, transects, and watershedsPoly
def USPTravel(): 
    arcpy.JoinField_management(os.path.join(path,"flowPaths"), "OBJECTID", os.path.join(path,"flowPathsUSP"), "OBJECTID_1", ["n","USP","TT"])
    arcpy.JoinField_management(os.path.join(path,"transects"), "OBJECTID", os.path.join(path,"transectsUSP"), "transectID", ["Qd","S","W","Dmax","A","Qprop","Wusp","V","L","USP","TT"])
    arcpy.JoinField_management(os.path.join(path,"watershedsPoly"), "FROM_ID", os.path.join(path,"watershedsTTDS"), "FROM_ID", ["TT","TTDS"])
    arcpy.Delete_management(os.path.join(path,"flowPathsUSP"))
    arcpy.Delete_management(os.path.join(path,"transectsUSP"))
    arcpy.Delete_management(os.path.join(path,"watershedsTTDS"))
    
# Flow path point elevation drop and flow accumulation
def flowPathPoints():
    arcpy.management.GeneratePointsAlongLines(os.path.join(path,"flowPathsRaw"), os.path.join(path,"flowPathPoints2"), "DISTANCE", "10 Meters", None, "END_POINTS")
    arcpy.DeleteField_management("flowPathPoints2", ["Shape_Length","FROM_ID","TO_ID","FlowAcc","RORUS","LengthM","n","USP","TT"])
    arcpy.AlterField_management("flowPathPoints2", "ORIG_FID", "pathID", "pathID")
    arcpy.AddField_management("flowPathPoints2", "pointID", "LONG")
    arcpy.CalculateField_management("flowPathPoints2", "pointID", "!OBJECTID!")
    arcpy.sa.ExtractValuesToPoints("flowPathPoints2", os.path.join(path,"demFill"), os.path.join(path,"flowPathPoints1"), "NONE", "VALUE_ONLY")
    arcpy.AlterField_management("flowPathPoints1", "RASTERVALU", "Elevation", "Elevation")
    arcpy.sa.ExtractValuesToPoints("flowPathPoints1", os.path.join(path,"D8FlowAcc"), os.path.join(path,"flowPathPoints"), "NONE", "VALUE_ONLY")
    arcpy.AlterField_management("flowPathPoints", "RASTERVALU", "FlowAcc", "FlowAcc")
    arcpy.AddField_management("flowPathPoints", "combID", "TEXT")
    arcpy.CalculateField_management("flowPathPoints", "combID", 'str(!pathID!) + "-" + str(!pointID!)', "PYTHON3", '', "TEXT")
    arcpy.AddField_management("flowPathPoints", "downID", "TEXT")
    arcpy.CalculateField_management("flowPathPoints", "downID", 'str(!pathID!) + "-" + str(!pointID!+1)', "PYTHON3", '', "TEXT")
    arcpy.management.CopyFeatures("flowPathPoints",os.path.join(path,"flowPathPointsCopy"))
    arcpy.AlterField_management("flowPathPointsCopy", "Elevation", "DownElevation", "DownElevation")
    arcpy.JoinField_management("flowPathPoints", "downID", "flowPathPointsCopy", "combID", ["DownElevation"])
    arcpy.management.SelectLayerByAttribute("flowPathPoints", "NEW_SELECTION", "DownElevation IS NULL", None)
    arcpy.management.DeleteFeatures("flowPathPoints")
    arcpy.AddField_management("flowPathPoints", "ElevDrop", "LONG")
    arcpy.management.CalculateField("flowPathPoints", "ElevDrop", "!Elevation! - !DownElevation!", "PYTHON3", '', "TEXT")
    arcpy.Delete_management(os.path.join(path,"flowPathPointsCopy"))
    arcpy.Delete_management(os.path.join(path,"flowPathPoints2"))
    arcpy.Delete_management(os.path.join(path,"flowPathPoints1"))
    
# Create a raster of cut depths around a stream to reduce slope to specified value            
def BankSlope(slope):
    arcpy.conversion.PolylineToRaster(os.path.join(path,"flowPathsMain"), "OBJECTID", os.path.join(path,"strm_ras"), "", "", 2)
    out_raster = arcpy.sa.ExtractByMask(os.path.join(path,"demCut"), os.path.join(path,"strm_ras"))
    out_raster.save(os.path.join(path,"strm_elev"))
    arcpy.env.parallelProcessingFactor = "100%"
    out_raster = arcpy.sa.Watershed(os.path.join(path,"D8FlowDir"), os.path.join(path, "strm_elev"), "Value")
    arcpy.env.parallelProcessingFactor = ""
    out_raster.save(os.path.join(path,"strm_shed"))
    out_raster = arcpy.sa.Minus(os.path.join(path,"demCut"), os.path.join(path,"strm_shed"))
    out_raster.save(os.path.join(path,"RelElev"))
    arcpy.analysis.Buffer(os.path.join(path,"flowPathsMain"), os.path.join(path,"buf50m"), "50 Meters", "FULL", "ROUND", "NONE", None, "PLANAR")
    arcpy.management.FeatureToLine("buf50m", os.path.join(path,"buf50mLine"), None, "ATTRIBUTES")
    arcpy.AddField_management("buf50mLine", "RelElev", "DOUBLE")
    arcpy.CalculateField_management("buf50mLine", "RelElev", float(50/slope))
    arcpy.management.Clip(os.path.join(path,"RelElev"), "", os.path.join(path,"RelElev50m"), "buf50m", 32767, "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    out_raster = arcpy.sa.Con(os.path.join(path,"RelElev50m"), 1, None, "Value <= 20")
    out_raster.save(os.path.join(path,"channel"))
    arcpy.conversion.RasterToPolygon(os.path.join(path,"channel"), os.path.join(path,"banksPoly"), "SIMPLIFY", "OBJECTID", "SINGLE_OUTER_PART", None)
    arcpy.management.SelectLayerByAttribute("banksPoly", "NEW_SELECTION", "Shape_Length = (SELECT MAX(Shape_Length) FROM banksPoly)", "INVERT")
    arcpy.DeleteFeatures_management("banksPoly")
    arcpy.management.FeatureToLine("banksPoly", os.path.join(path,"banksLine"), None, "ATTRIBUTES")
    arcpy.management.SelectLayerByAttribute("banksLine", "NEW_SELECTION", "Shape_Length = (SELECT MAX(Shape_Length) FROM banksLine)", "INVERT")
    arcpy.DeleteFeatures_management("banksLine")
    arcpy.AddField_management("banksLine", "RelElev", "DOUBLE")
    arcpy.CalculateField_management("banksLine", "RelElev", 0)
    arcpy.env.parallelProcessingFactor = "100%"
    slope = arcpy.sa.TopoToRaster("banksLine RelElev Contour;buf50mLine RelElev Contour;buf50m # Boundary", 2, "", 20, None, None, "ENFORCE", "CONTOUR", 20, None, 1, 0, 2.5, 100, None, None, None, None, None, None, None, None)
    arcpy.env.parallelProcessingFactor = ""
    slopecm = arcpy.sa.Times("slope", 100)
    cutRaw = arcpy.sa.Minus("RelElev50m", "slopecm")
    cut = arcpy.sa.Con("cutRaw", 0, "cutRaw", "VALUE <= 0")
    cut.save(os.path.join(path,"cut"))
 
# Select watersheds that are upstream of a selected watershed
def selectUpstream():
    i = 0
    j = 1
    FROMS = sorted(set([r[0] for r in arcpy.da.SearchCursor("watershedsPoly", "FROM_ID")]))
    while True:
        query = "TO_ID IN (" + str(FROMS)[1:-1] + ")"
        arcpy.SelectLayerByAttribute_management("watershedsPoly", "ADD_TO_SELECTION", query, None)
        FROMS = sorted(set([r[0] for r in arcpy.da.SearchCursor("watershedsPoly", "FROM_ID")]))
        j = int(arcpy.GetCount_management("watershedsPoly")[0])
        if j == i:
            break
        else:
            i = j      

# Select flow paths that are downstream of a selected flow path
def selectDownstream():
    i = 0
    j = 1
    TOS = sorted(set([r[0] for r in arcpy.da.SearchCursor("flowPaths", "TO_ID")]))
    while True:
        query = "FROM_ID IN (" + str(TOS)[1:-1] + ")"
        arcpy.SelectLayerByAttribute_management("flowPaths", "ADD_TO_SELECTION", query, None)
        TOS = sorted(set([r[0] for r in arcpy.da.SearchCursor("flowPaths", "TO_ID")]))
        j = int(arcpy.GetCount_management("flowPaths")[0])
        if j == i:
            break
        else:
            i = j                  
            
# Select watersheds that are upstream of a selected watershed in the same HUC12
def selectUpstreamInHUC12(HUC12):
    i = 0
    j = 1
    FROMS = sorted(set([r[0] for r in arcpy.da.SearchCursor("ACPFDaneWatersheds", "FROM_ID")]))
    while True:
        query = f"TO_ID IN (" + str(FROMS)[1:-1] + ") And HUC12 = " + "\'" + HUC12 + "\'"
        arcpy.SelectLayerByAttribute_management("ACPFDaneWatersheds", "ADD_TO_SELECTION", query, None)
        FROMS = sorted(set([r[0] for r in arcpy.da.SearchCursor("ACPFDaneWatersheds", "FROM_ID")]))
        j = int(arcpy.GetCount_management("ACPFDaneWatersheds")[0])
        if j == i:
            break
        else:
            i = j       
            
