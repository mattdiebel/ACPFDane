# Follow this script to run the ACPFDane functions for a single HUC12 watershed.
# The functions were written for ArcGIS Pro 2.8.

exec(open("functions.py").read())

# Add the following layers (with the specified names) to an ArcGIS Pro project:
    # "dem" (uncut, unfilled) for ACPF buffered watershed
    # "cutLines" (culvert lines)
    # "stormsewers" (storm sewer lines)
    # "CNlow" (runoff curve number raster)
    # "gSSURGO" (version from ACPF database)
    
# Specify path to geodatabase:
path = "acpf.gdb"

# Run the following functions in order.    
# Comments between function calls are manual steps.
# Functions without arguments use input layers and outputs from previous functions

conditionDEM()

flowPaths(AreaThreshold) # AreaThreshold is flow accumulation (acres) threshold for flow paths

Depressions(MinimumSize) # MinimumSize is minimum depression area (acres)

# Leave "stub" flowLines on all inlets and outlets, and delete flowLines that connect to those stubs
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

runoffEvent()

# Map RORDS to check topology

makeTransects()

# R Manning Test Script

USPTravel()
