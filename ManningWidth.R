# Solve Manning equation for width with known discharge and irregular channel geometry
# xsec is a data frame of X and Y coordinates of the channel cross section in meters
# Qd is the design discharge in cfs (typically bankfull)
# S is the channel slope
# n is Manning's n (roughness coefficient)

ManningWidth <- function(xsec, Qd, S, n = 0.04) {
  
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
