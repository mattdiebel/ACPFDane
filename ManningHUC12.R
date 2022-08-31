# Calculates hydraulic properties for transects in an ACPF HUC12 database

ManningHUC12 <- function(gdbpath) {
  
  # Read in datasets
  flowPathsFC = arc.open(paste0(gdbpath,'/flowPaths'))
  flowPaths = arc.select(flowPathsFC)
  flowPaths$n = NA
  
  transectPoints = arc.open(paste0(gdbpath,'/transectPoints'))
  data = arc.select(transectPoints)
  
  watershedsFC = arc.open(paste0(gdbpath,'/watershedsPoly'))
  watersheds = arc.select(watershedsFC)
  
  # Count transects per flow path and fill in null flow accumulation and RORUS values from sum of upstream reaches
  for (p in 1:nrow(flowPaths)) {
    flowPaths$n[p] = length(unique(data$transectID[data$pathID==flowPaths$OBJECTID[p]]))
    if (is.na(flowPaths$FlowAcc[p])) {
      flowPaths$FlowAcc[p] = sum(flowPaths$FlowAcc[flowPaths$TO_ID==flowPaths$FROM_ID[p]], na.rm = TRUE)
    }
    if (is.na(flowPaths$RORUS[p])) {
      flowPaths$RORUS[p] = mean(flowPaths$RORUS[flowPaths$TO_ID==flowPaths$FROM_ID[p]], na.rm = TRUE)
    }
  }
  
  flowPaths$RORUS[flowPaths$RORUS<0.01] = 0.01 # Avoid very long travel times from very low Qd
  
  # USGS 2-year recurrence interval discharge
  flowPaths$Qd = 35 * (flowPaths$FlowAcc * 4 / 2590000 * flowPaths$RORUS) ^ 0.649 
  
  all_transects = data.frame()
  pct_old = 0
  print(gdbpath)
  print(paste0(nrow(flowPaths)," flow paths"))
  flush.console()
  
  for (p in 1:nrow(flowPaths)) {
    n = flowPaths$n[p]
    remainder = flowPaths$LengthM[p] - (flowPaths$n[p]-2)*100
    remainderCM = remainder * 100
    pdata = data[data$pathID==flowPaths$OBJECTID[p],]
    transects = data.frame(pathID = flowPaths$OBJECTID[p],
                           transectID = unique(pdata$transectID),
                           order = NA,
                           Elevation = NA,
                           Qd = flowPaths$Qd[p],
                           S = NA,
                           W = NA,
                           Dmax = NA,
                           A = NA,
                           Qprop = NA,
                           Wusp = NA,
                           V = NA,
                           L = 100)
    
    # Sort transects downstream (Arc tool does the end transects last (e.g., 2,3,4,1,5))
    transects$order = rank(transects$transectID, ties.method = "first")
    transects = transects[order(transects$order),]
    transects$order[nrow(transects)-1] = 0
    transects = transects[order(transects$order),]
    
    
    # Copy channel elevation to transects
    for (t in 1:(nrow(transects))) {
      tdata = data[data$transectID==transects$transectID[t],]
      tdata = tdata[order(tdata$pointID),]
      transects$Elevation[t] = tdata$Elevation[9]
    }
    
    # Calculate cross section dimensions at Qd
    for (t in 1:(nrow(transects))) {
      
      # Calculate slope
      if(n==2) {
        transects$L[t] = remainder/2
        transects$S[t] = (transects$Elevation[1] - transects$Elevation[2]) / remainderCM
      }
      
      if(n>2) {
        if (t==1) {
          transects$S[t] = (transects$Elevation[t] - transects$Elevation[t+1]) / 10000
        }
        if (t>1 & t<n) {
          transects$S[t] = (transects$Elevation[t-1] - transects$Elevation[t]) / 10000
        }
        if (t==n) {
          transects$S[t] = (transects$Elevation[t-1] - transects$Elevation[t]) / remainderCM
        }
        if (t>=(n-1)) {
          transects$L[t] = remainder/2
        }
      }
      
      if (transects$S[t]<=0) {transects$S[t] = 0.0001}
      
      # Format cross section geometry for ManningWidth function
      tdata = data[data$transectID==transects$transectID[t],]
      tdata = tdata[order(tdata$pointID),]
      tdata = tdata[1:17,]
      tdata$X = c(0:16*2)
      tdata$Y = tdata$Elevation/100
      tdata = tdata[,c("X","Y")]
      
      # Run Manning Width function
      MW = ManningWidth(tdata, transects$Qd[t], transects$S[t], 0.04)
      transects$W[t] = MW[1]
      transects$A[t] = MW[2]
      transects$Dmax[t] = MW[3]
      transects$Wusp[t] = MW[4]
      transects$Qprop[t] = MW[5]
      transects$V[t] = transects$Qd[t]/transects$A[t]
      transects$TT[t] = transects$L[t]*3.28084/transects$V[t]/60/60
    }
    
    # Calculate unit stream power
    transects = transects[!is.na(transects$W),]
    transects$USP = (1000 * 9.8 * transects$Qd/3.28084^3 * transects$Qprop * transects$S) / (transects$Wusp / 3.28084)
    all_transects = rbind(all_transects, transects)
    
    pct = p/nrow(flowPaths)*100
    pct = round(pct,0)
    if (pct >= pct_old+10) {
      pct_old = pct
      print(paste0(pct,"%"))
      flush.console()
    }

    
  }
  
  all_transects$USPL = all_transects$USP * all_transects$L
  flowPaths$USP = NA
  for (p in 1:nrow(flowPaths)) {
    pt = all_transects[all_transects$pathID==flowPaths$OBJECTID[p],]
    flowPaths$USP[p] = sum(pt$USPL, na.rm = TRUE) / sum(pt$L, na.rm = TRUE)
    flowPaths$TT[p] = sum(pt$TT, na.rm = TRUE)
  }
  
  if(!("DSorder" %in% colnames(watersheds))) {
    ds = c()
    froms = watersheds$FROM_ID[is.na(watersheds$TO_ID)]
    ds = append(ds,froms)
    len_start = length(ds)
    repeat {
      len_start = length(ds)
      froms = watersheds$FROM_ID[watersheds$TO_ID %in% froms]
      ds = append(ds,froms)
      len_end = length(ds)
      if(len_end==len_start) {break}
    }
    DSorder = data.frame("FROM_ID" = ds, "DSorder" = NA)
    DSorder$DSorder = as.numeric(row.names(DSorder))
    watersheds = merge(watersheds, DSorder, by="FROM_ID", all.x = TRUE)
  }
  
  watersheds = merge(watersheds[,c("FROM_ID","TO_ID","DSorder")], flowPaths[,c("FROM_ID","TT")], by="FROM_ID", all.x = TRUE)
  watersheds = watersheds[order(watersheds$DSorder),]
  watersheds$TT[is.na(watersheds$TT)] = 0
  watersheds$TTDS = 0
  watersheds$TO_ID[is.na(watersheds$TO_ID)] = 0
  
  for (w in 1:nrow(watersheds)) {
    DS = watersheds[watersheds$FROM_ID==watersheds$TO_ID[w],]
    if (nrow(DS)==0) watersheds$TTDS[w] = watersheds$TT[w] else
      watersheds$TTDS[w] = watersheds$TT[w] + DS$TTDS[1]
  }
  
  arc.write(paste0(gdbpath,'/transectsUSP'), all_transects, overwrite = TRUE)
  arc.write(paste0(gdbpath,'/flowPathsUSP'), flowPaths, overwrite = TRUE)
  arc.write(paste0(gdbpath,'/watershedsTTDS'), watersheds, overwrite = TRUE)
  
}
