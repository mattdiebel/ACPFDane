# ACPFDane
Modified version of the Agricultural Conservation Planning Framework [ACPF](https://acpf4watersheds.org/) for Dane County, Wisconsin

This script and set of functions takes a HUC12 LiDAR DEM as input and conducts the following analyses:
1. Condition DEM with culverts and storm sewers, fill depressions, calculate flow direction and accumulation.
2. Define flow paths and topographic depressions
3. Delineate flow path and depression watersheds
4. Define watershed topology from node intersections and LakeCat tool
5. Calculate runoff from a design storm with curve number hydrology
6. Calculate unit stream power of flow paths

This presentation provides more conceptual background and an example workflow: Concepts and Tools for Hydrologic Modeling with LiDAR.pdf
