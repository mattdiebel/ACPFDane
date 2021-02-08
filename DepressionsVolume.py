# Surface Depressions
#
# Original coding: S. Porter - NLAE 04/2014
# Modified by M. Diebel, 12/2019

# Import System Modules
import arcpy
from arcpy import env
from arcpy.sa import *
import sys, string, os, os.path, time
arcpy.CheckOutExtension("Spatial")

def DepressionIdentification(demCut, demFill, gSSURGO, ZFactor, MinimumSize, OutDepressions):

    # Fill DEM and find sinks
    arcpy.AddMessage("Finding depressions")
    FillReg = Minus(demFill, demCut)

    # Convert values < 0 to NODATA (non-depressional) else 1 
    AllSinks = Con(FillReg, 1, "", "VALUE > .01")
    
    # Convert depression raster to polygons
    Depressions = arcpy.RasterToPolygon_conversion(AllSinks, OutDepressions, "NO_SIMPLIFY")

    # Convert minimum surface area from acres to meters squared
    arcpy.AddMessage("Deleting depressions less than %s acres" % MinimumSize)
    MinSize = (float(MinimumSize) * float(4046.00))

    # Delete depressions smaller than minimum size
    arcpy.MakeFeatureLayer_management(Depressions, "small_sinks", "\"Shape_Area\" < %s" % MinSize)
    arcpy.DeleteFeatures_management("small_sinks")

    # Add and populate unique "Depress_id" field, populate with 10,000 + OID to avoid duplicates with stream ids
    desc = arcpy.Describe(Depressions)
    OID = desc.OIDFieldName
    arcpy.AddField_management(Depressions, "Depress_ID")
    arcpy.CalculateField_management(Depressions, "Depress_ID", "10000 + !%s!" % OID, "PYTHON") 

    # Calculate mean percent hydric soil within each depression - zonal mean
    arcpy.AddMessage("Calculating percent hydric soil")
    gSSURGO_copy = arcpy.CopyRaster_management(gSSURGO, "gSSURGO_copy", "", 0, 0)
    arcpy.AddField_management(gSSURGO_copy,"HydFlt","FLOAT","8","0")
    arcpy.CalculateField_management(gSSURGO_copy,"HydFlt", '!Hydric!', "PYTHON")
    HydricPct = Lookup(gSSURGO_copy, "HydFlt")
    HydricTable = ZonalStatisticsAsTable(Depressions, "Depress_ID", HydricPct, "in_memory" + "\\HydricTable", "#", "MEAN")

    # Join and populate "PctHydric" field
    arcpy.JoinField_management(Depressions, "Depress_ID", HydricTable, "Depress_ID", ["MEAN"])
    arcpy.AddField_management(Depressions, "PctHydric", "FLOAT")
    arcpy.CalculateField_management(Depressions, "PctHydric", "!MEAN!", "PYTHON")
    arcpy.DeleteField_management(Depressions, ["MEAN"])

    # Calculate max depth in cm 
    arcpy.AddMessage("Calculating maximum depth in cm")
    DepthStats = ZonalStatisticsAsTable(Depressions, "Depress_ID", FillReg, "in_memory" + "\\DepthStats", "", "MAXIMUM")

    # Add and calculate MaxDepthCM field 
    arcpy.JoinField_management(Depressions, "Depress_ID", DepthStats, "Depress_ID", ["MAX"])
    arcpy.AddField_management(Depressions, "MaxDepthCM")
    arcpy.CalculateField_management(Depressions, "MaxDepthCM", "(!MAX! * (%s * 100))" % ZFactor, "PYTHON")
    arcpy.DeleteField_management(Depressions, ["gridcode", "MAX", "id"])
    
    # Calculate volume in acre-ft
    arcpy.AddMessage("Calculating volume in acre-ft")
    VolStats = ZonalStatisticsAsTable(Depressions, "Depress_ID", FillReg, "in_memory" + "\\VolStats", "", "SUM")

    # Add and calculate VolAcreFt field 
    arcpy.JoinField_management(Depressions, "Depress_ID", VolStats, "Depress_ID", ["SUM"])
    arcpy.AddField_management(Depressions, "VolAcreFt", "FLOAT")
    
    SumDepthtoAcFt = 1233.48 / float(ZFactor) / resolution
    arcpy.CalculateField_management(Depressions, "VolAcreFt", "(!SUM! / %s)" % SumDepthtoAcFt, "PYTHON")
    arcpy.DeleteField_management(Depressions, ["SUM"])

    # Cleanup
    arcpy.Delete_management(gSSURGO_copy)
    arcpy.Delete_management(HydricPct)
    arcpy.Delete_management(HydricTable)
    arcpy.Delete_management(DepthStats)
    arcpy.Delete_management(VolStats)
    arcpy.Delete_management(demFill)
    arcpy.Delete_management(FillReg)
    arcpy.Delete_management(AllSinks)
    del[demFill, FillReg, AllSinks, DepthStats, VolStats, gSSURGO_copy, HydricPct, HydricTable]
    arcpy.Delete_management("in_memory")

##############################################################################
    
if __name__ == "__main__":

    # Define Parameters
    demCut = arcpy.GetParameterAsText(0)
    demFill = arcpy.GetParameterAsText(1)
    gSSURGO = arcpy.GetParameterAsText(2)
    ZFactor = arcpy.GetParameterAsText(3)
    MinimumSize = arcpy.GetParameterAsText(4)
    OutDepressions = arcpy.GetParameterAsText(5)

    # Set environments
    arcpy.env.snapRaster = demCut
    rObj = Raster(demCut)
    arcpy.env.cellSize = rObj.meanCellHeight
    arcpy.env.overwriteOutput = True
    arcpy.env.outputCoordinateSystem = demCut
    
    # Determine the cellsize and resolution of DEM
    cellsize = float(arcpy.GetRasterProperties_management(demCut, "CELLSIZEX").getOutput(0))
    resolution = float(cellsize * cellsize)

    if not arcpy.Exists(env.workspace):
        arcpy.AddError("workspace does not exist!! Please set your workspace to a valid path directory in Arcmap --> Geoprocessing --> Environments --> Workspace")
        sys.exit(0)

    # Run modules
    DepressionIdentification(demCut, demFill, gSSURGO, ZFactor, MinimumSize, OutDepressions)

