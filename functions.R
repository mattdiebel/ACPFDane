# This file defines two R functions for calculating hydraulic properties of flow paths in an ACPF HUC12 database
# ManningHUC12 reads in transects, flow paths, and watersheds from the ACPF database, runs the ManningWidth function, and writes the results back to the ACPF database.
# ManningWidth solves the Manning equation for width, depth, and velocity with known discharge, slope, roughness, and irregular channel geometry

ManningHUC12 <- function(gdbpath) {
  
  # gdbpath is the file path to the ACPF geodatabase
  
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

#######################################################################################

ManningWidth <- function(xsec, Qd, S, n = 0.04) {
  
  # xsec is a data frame of X and Y coordinates of the channel cross section in meters
  # Qd is the design discharge in cfs (typically bankfull)
  # S is the channel slope
  # n is Manning's n (roughness coefficient)
  
  data = xsec
  dX = data$X[2] - data$X[1]
  ends = data.frame(X=c(min(data$X)-dX,max(data$X)+dX), Y=rep(max(data$Y)+1,2))
  data = rbind(data,ends)
  data$Yr = rank(data$Y, ties.method = "first")
  data = data[order(data$X),]
  data$Q = NA
  data$Q[data$Yr==1] = 0
  
  # Calculate Q at elevation of each transect point
  for (i in 2:nrow(data)) {
    
    if(data$Y[data$Yr==i] == data$Y[data$Yr==i-1]) {
      data$Q[data$Yr==i] = data$Q[data$Yr==i-1]
      next
    }
    
    dati = data
    Yi = dati$Y[dati$Yr==i]
    dati$D = Yi - dati$Y
    dati = dati[,c(1,2,5)]
    datij = data.frame()
  
    for (j in 1:(nrow(dati)-1)) {
      if(dati$D[j]>0 & dati$D[j+1]<0) {
        datj = dati[1,]
        datj[1,] = NA
        datj$X = dati$X[j] + dati$D[j]/(dati$D[j]-dati$D[j+1])*dX
        datj$Y = Yi
        datj$D = 0
        datij = rbind(datij, datj)
      }
      if(dati$D[j]<0 & dati$D[j+1]>0) {
        datj = dati[1,]
        datj[1,] = NA
        datj$X = dati$X[j+1] - dati$D[j+1]/(dati$D[j+1]-dati$D[j])*dX
        datj$Y = Yi
        datj$D = 0
        datij = rbind(datij, datj)
      }
    }
    
    dati = rbind(dati, datij)
    dati = dati[order(dati$X),]
    
    dati$W = 0
    dati$H = 0
    dati$A = 0
    dati$P = 0
    
    for (j in 1:(nrow(dati)-1)) {
      if(dati$D[j]<0 | dati$D[j+1]<0) {next} else {
        dati$W[j] = dati$X[j+1] - dati$X[j]
        dati$H[j] = (dati$D[j+1] + dati$D[j]) / 2
        dati$A[j] = dati$W[j] * dati$H[j]
        dati$P[j] = sqrt((dati$X[j+1] - dati$X[j])^2 + (dati$Y[j+1] - dati$Y[j])^2)
      }
    }
    
    W = sum(dati$W) * 3.28084
    Dmax = max(dati$D) * 3.28084
    A = sum(dati$A) * 3.28084^2
    Amax = max(dati$A) * 3.28084^2
    Qprop = Amax / A
    Wusp = dati$W[dati$A==max(dati$A)] * 3.28084
    Wusp = Wusp[1]
    P = sum(dati$P) * 3.28084
    R = A/P
    Q = (1.49/n)*A*R^(2/3)*sqrt(S)
    
    data$Q[data$Yr==i] = Q
    if(Q>Qd) {break}
  }
  
  Qdiff = min(abs(Qd - data$Q), na.rm = TRUE)
  Ymax = data$Y[data$Yr==max(data$Yr[!is.na(data$Q)])]
  Ymin = data$Y[data$Yr==(max(data$Yr[!is.na(data$Q)])-1)]
  
  
  # If Qd fills the transect, return full transect dimensions
  if(max(data$Q, na.rm = TRUE) < Qd) {
    return(c(W,A,Dmax,Wusp,Qprop))
    break
  }
  
  # If Qd does not fill the transect, adjust the flow elevation until Q is within 1% of Qd
  while (Qdiff > 0.01 * Qd) {
    
    Ystep = 0.1 * (Ymax - Ymin)
    Ys = data.frame(Y = Ymin + c(0:10*Ystep), Q = NA)
    
    for (i in 1:nrow(Ys)) {
      dati = data
      dati$D = Ys$Y[i] - dati$Y
      dati = dati[,c(1,2,5)]
      datij = data.frame()
      
      for (j in 1:(nrow(dati)-1)) {
        if(dati$D[j]>0 & dati$D[j+1]<0) {
          datj = dati[1,]
          datj[1,] = NA
          datj$X = dati$X[j] + dati$D[j]/(dati$D[j]-dati$D[j+1])*dX
          datj$Y = Ys$Y[i]
          datj$D = 0
          datij = rbind(datij, datj)
        }
        if(dati$D[j]<0 & dati$D[j+1]>0) {
          datj = dati[1,]
          datj[1,] = NA
          datj$X = dati$X[j+1] - dati$D[j+1]/(dati$D[j+1]-dati$D[j])*dX
          datj$Y = Ys$Y[i]
          datj$D = 0
          datij = rbind(datij, datj)
        }
      }
      
      dati = rbind(dati, datij)
      dati = dati[order(dati$X),]
      
      dati$W = 0
      dati$H = 0
      dati$A = 0
      dati$P = 0
      
      for (j in 1:(nrow(dati)-1)) {
        if(dati$D[j]<0 | dati$D[j+1]<0) {
          next
          } else {
          dati$W[j] = dati$X[j+1] - dati$X[j]
          dati$H[j] = (dati$D[j+1] + dati$D[j]) / 2
          dati$A[j] = dati$W[j] * dati$H[j]
          dati$P[j] = sqrt((dati$X[j+1] - dati$X[j])^2 + (dati$Y[j+1] - dati$Y[j])^2)
        }
      }
      
      W = sum(dati$W) * 3.28084
      Dmax = max(dati$D) * 3.28084
      A = sum(dati$A) * 3.28084^2
      Amax = max(dati$A) * 3.28084^2
      Qprop = Amax / A
      Wusp = dati$W[dati$A==max(dati$A)] * 3.28084
      Wusp = Wusp[1]
      P = sum(dati$P) * 3.28084
      R = A/P
      if (A==0) {
        Q = 0
        } else {
          Q = (1.49/n)*A*R^(2/3)*sqrt(S)
        }
      
      Ys$W[i] = W
      Ys$Dmax[i] = Dmax
      Ys$A[i] = A
      Ys$Amax[i] = Amax
      Ys$Qprop[i] = Qprop
      Ys$Wusp[i] = Wusp
      Ys$Q[i] = Q
      if(Q>Qd) {break}
    }
    
    Ys = Ys[!is.na(Ys$Q),]
    Ys$Qdiff = abs(Qd - Ys$Q)
    W = Ys$W[Ys$Qdiff==min(Ys$Qdiff)]
    Dmax = Ys$Dmax[Ys$Qdiff==min(Ys$Qdiff)]
    A = Ys$A[Ys$Qdiff==min(Ys$Qdiff)]
    Wusp = Ys$Wusp[Ys$Qdiff==min(Ys$Qdiff)]
    Qprop = Ys$Qprop[Ys$Qdiff==min(Ys$Qdiff)]
    
    Qdiff = min(Ys$Qdiff)
    Ymax = Ys$Y[nrow(Ys)]
    Ymin = Ys$Y[nrow(Ys)-1]
    
  }
  
  return(c(W,A,Dmax,Wusp,Qprop))
}
