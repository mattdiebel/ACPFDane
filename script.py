# This script was developed for ArcGIS 2.3
# Some of the functions in this script need memory management improvements - they tend to cause ArcGIS to crash. Run one function at a time if this happens.
# Comments between function calls are manual steps.

# Add the following layers (with the specified names) to an ArcGIS Pro project:
    # "dem" (uncut, unfilled) for ACPF buffered watershed
    # "cutLines" (culvert lines)
    # "depOutlets" (depression outlet points from preliminary analysis, if applicable)
    # "CNlow" (runoff curve number raster)
    # "gSSURGO" (version from ACPF database)
    
# Specify path to geodatabase:
path = "acpf.gdb"

exec(open("functions.py").read())

arcpy.HydroConditioning("cutLines", None, "dem", os.path.join(path, "demCut"), os.path.join(path, "demFill"), os.path.join(path, "D8FlowDir"), os.path.join(path, "D8FlowAcc"), os.path.join(path, "Hshd"), 0.01)

arcpy.FlowPaths("D8FlowAcc", "D8FlowDir", 10, os.path.join(path, "flowPaths"))

arcpy.DepressionVolume("demCut", "demFill", "gSSURGO", 0.01, 0.459, os.path.join(path, "Depressions"))

deleteFalseDepressions() # Only run this function if depOutlets is present

# Leave "stub" flowLines on all inlets and outlets, and delete flowLines that connect to those stubs; Save Edits
# Select any reach in watershed from flowPaths

findConnected()

# Check selection for accuracy

refineHUC12()

# Delete pathWatersheds that are upstream of watershed inlet or downstream of watershed outlet

pruneDepressions()

pruneFlowPaths()

# Sort flowPaths by DepressID, then for each DepressID:
    # Delete all but one flowPath that exits the depression, then
    # Merge remaining path with next downstream path (i.e., remove pseudonodes)

defineTopology()

watershedSeeds()

watershed()

watershedPolygons()

LakeCat()

watershedAttributes()

# Delete watershedsPoly and flowPaths that are upstream or downstream of watershed boundary.
# Edit TO_ID of outlet flowPath and watershedsPoly to Null
# Check that all watersheds except outlet have a TO_ID

runoff()

# Map RORDS to check topology

makeTransects()

# R Manning Test Script

joinUSP()
