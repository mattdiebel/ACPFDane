setwd("acpf.gdb")
source("ManningWidth.R")

library(arcgisbinding)
arc.check_product()

in_param_flowPaths = 'flowPaths'
in_param_transectPoints = 'transectPoints2'
out_param_flowPaths = 'flowPathsUSP'
out_param_transects = 'transectsUSP'

flowPathsFC = arc.open(in_param_flowPaths)
flowPaths = arc.select(flowPathsFC)
flowPaths$n = NA
flowPaths$USP = NA
flowPaths$RORUS[is.na(flowPaths$RORUS)] = 0

transectPoints = arc.open(in_param_transectPoints)
data = arc.select(transectPoints)

for (p in 1:nrow(flowPaths)) {
  flowPaths$n[p] = length(unique(data$transectID[data$pathID==flowPaths$OBJECTID[p]]))
  if (is.na(flowPaths$FlowAcc[p])) {
    flowPaths$FlowAcc[p] = sum(flowPaths$FlowAcc[flowPaths$TO_ID==flowPaths$FROM_ID[p]], na.rm = TRUE)
  }
}

all_transects = data.frame()

for (p in 1:nrow(flowPaths)) {
  remainder = flowPaths$LengthM[p] - (flowPaths$n[p]-2)*100
  pdata = data[data$pathID==flowPaths$OBJECTID[p],]
  transects = data.frame(pathID = flowPaths$OBJECTID[p],
                         transectID = unique(pdata$transectID),
                         order = NA,
                         Elevation = NA,
                         FlowAcc = flowPaths$FlowAcc[p],
                         RORUS = flowPaths$RORUS[p],
                         Qd = NA,
                         S = NA,
                         W = NA,
                         L = 100)

  # Sort transects downstream (Arc tool does the end transects last (e.g., 2,3,4,1,5))
  transects$order = rank(transects$transectID, ties.method = "first")
  transects = transects[order(transects$order),]
  transects$order[nrow(transects)-1] = 0
  transects = transects[order(transects$order),]

  # Copy channel elevation to transects
  for (t in 1:nrow(transects)) {
    tdata = data[data$transectID==transects$transectID[t],]
    tdata = tdata[order(tdata$pointID),]
    transects$Elevation[t] = tdata$Elevation[9]
  }

  # Calculate wetted width for each transect
  for (t in 1:(nrow(transects))) {
    if(flowPaths$n[p]==2) {
      A = transects$FlowAcc[1] * 4 / 2590000 * transects$RORUS[1]
      if (t==1) {
        transects$S[t] = (transects$Elevation[t] - transects$Elevation[t+1]) / (100 * remainder)
        transects$L[t] = remainder/2
      } else {
        transects$S[t] = (transects$Elevation[t-1] - transects$Elevation[t]) / (100 * remainder)
        transects$L[t] = remainder/2
      }
    }

    if(flowPaths$n[p]>2) {
      if (t==1) {
        A = transects$FlowAcc[1] * 4 / 2590000 * transects$RORUS[1]
        transects$S[t] = (transects$Elevation[t] - transects$Elevation[t+1]) / 10000
      }
      if (t>1 & t<flowPaths$n[p]) {
        A = transects$FlowAcc[t] * 4 / 2590000 * transects$RORUS[t]
        transects$S[t] = (transects$Elevation[t-1] - transects$Elevation[t]) / 10000
      }
      if (t==flowPaths$n[p]) {
        A = transects$FlowAcc[t-1] * 4 / 2590000 * transects$RORUS[t-1]
        transects$S[t] = (transects$Elevation[t-1] - transects$Elevation[t]) / (100 * remainder)
        transects$L[t] = remainder
      }
    }

    if (transects$S[t]<=0) {transects$S[t] = 0.0001}
    transects$Qd[t] = 35 * A^0.649 
    tdata = data[data$transectID==transects$transectID[t],]
    tdata = tdata[order(tdata$pointID),]
    tdata = tdata[1:17,]
    tdata$X = c(0:16*2)
    tdata$Y = tdata$Elevation/100
    tdata = tdata[,c("X","Y")]
    transects$W[t] = ManningWidth(tdata, transects$Qd[t], transects$S[t], 0.04)
  }

  # Calculate unit stream power
  transects = transects[!is.na(transects$W),]
  transects$USP = (1000 * 9.8 * transects$Qd/3.28084^3 * transects$S) / (transects$W/3.28084)
  all_transects = rbind(all_transects, transects)

  pct = round(p/nrow(flowPaths),3)
  print(pct)
  flush.console()

}

all_transects$USPL = all_transects$USP * all_transects$L
for (p in 1:nrow(flowPaths)) {
  pt = all_transects[all_transects$pathID==flowPaths$OBJECTID[p],]
  flowPaths$USP[p] = sum(pt$USPL, na.rm = TRUE) / sum(pt$L, na.rm = TRUE)
}

arc.write(out_param_transects, all_transects, overwrite = TRUE)
arc.write(out_param_flowPaths, flowPaths, overwrite = TRUE)
  

# Test arguments for ManningWidth function
# xsec = tdata
# Qd = transects$Qd[t]
# S = transects$S[t]
# n = 0.04


# USGS SIR 2016-5140 Area 8 Q50p equation
# A = 10/640 # (square miles)
# Ksat = 10.1 # (inches per hour)
# LUw = 1.59 # (percent open water)
# Q50p = 150 * A^0.649 * Ksat^-0.895 * ((LUw+0.01)/100)^-0.147

