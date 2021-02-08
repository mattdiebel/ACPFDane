# Flow Paths ----------------------------------------------------------------------------------------------
#
# Defines a flow path network using an area threshold in acres
#
##--------------------------------------------------------------------------------
## Original coding: S. Porter - NLAE 11/2017
## Revised by M. Diebel 12/2019
##--------------------------------------------------------------------------------

# Import System Modules
import arcpy
from arcpy import env
from arcpy.sa import *
import sys
arcpy.CheckOutExtension("Spatial")

def FlowPaths(D8FlowAcc, D8FlowDir, AreaThreshold):

    # Convert area threshold from acres (input) to meters: 1 acre = 4046 meters
    thresh_meters = float(AreaThreshold) * 4046
    number_cells = float(thresh_meters / resolution)
    arcpy.AddMessage("Area threshold of %s acres...." % (AreaThreshold))
   
    # Threshold the flow accumulation raster using a CON statement and convert 0 background values to null
    FlowNetRas = Con(D8FlowAcc, 1, "", "VALUE >= %s" % number_cells)
    
    # Assigns unique values to sections of a raster linear network between intersections
    Link = StreamLink(FlowNetRas, D8FlowDir)
    
    # Converts a raster representing a linear network to features representing the linear network
    StreamToFeature(Link, D8FlowDir, flowPaths, "NO_SIMPLIFY")
    
    # Add downstream grid_code to flowPaths
    flowPathsCopy = arcpy.management.CopyFeatures(flowPaths)
    arcpy.AddField_management(flowPathsCopy, "grid_code_to", "LONG")
    arcpy.CalculateField_management(flowPathsCopy, "grid_code_to", "!grid_code!") 
    arcpy.JoinField_management(flowPaths, "to_node", flowPathsCopy, "from_node", ["grid_code_to"])
        
##############################################################################
    
if __name__ == "__main__":

    # Define Parameters
    D8FlowAcc = arcpy.GetParameterAsText(0)  
    D8FlowDir = arcpy.GetParameterAsText(1)
    AreaThreshold = arcpy.GetParameterAsText(2)
    flowPaths = arcpy.GetParameterAsText(3)

    # Set environments
    arcpy.env.extent = D8FlowAcc
    arcpy.env.snapRaster = D8FlowAcc

    # Determine the cellsize and resolution of input flow accumulation raster
    cellsize = float(arcpy.GetRasterProperties_management(D8FlowAcc, "CELLSIZEX").getOutput(0))
    resolution = float(cellsize * cellsize)

    if not arcpy.Exists(env.workspace):
        arcpy.AddError("workspace does not exist!! Please set your workspace to a valid path directory in Arcmap --> Geoprocessing --> Environments --> Workspace")
        sys.exit(0)

    # Run Modules
    FlowPaths(D8FlowAcc, D8FlowDir, AreaThreshold)
    

