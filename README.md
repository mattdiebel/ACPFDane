# ACPFDane
Modified version of the Agricultural Conservation Planning Framework for Dane County, Wisconsin

This script takes a HUC12 LiDAR DEM as input and conducts the following analyses:
1. Run ACPF tools:
    - Hydro conditioning: Cut DEM with culvert lines, fill depressions, calculate flow direction and accumulation.
    - Flow paths: Define flow paths with 10 acre flow accumulation.
    - Depressions: Define depressions in DEM above 20,000 sq ft in surface area.
2. Refine HUC12 boundary
3. Delineate flow path and depression watersheds
4. Define watershed topology through node intersections and LakeCat tool
5. Loop through design storms
    - Loop through watersheds to calculate incremental and cumulative runoff with curve number hydrology
    - Loop through watersheds to calculate downstream runoff ratio
6. Calculate unit stream power of flow paths
