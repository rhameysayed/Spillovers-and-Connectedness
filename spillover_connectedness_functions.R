# install.packages("vars")
# install.packages("igraph")
# install.packages("PerformanceAnalytics")
# install.packages("quantmod")
# install.packages("plyr")
# install.packages("xts")
# install.packages("reshape2")
# install.packages("dplyr")

library(vars)
library(igraph)
library(PerformanceAnalytics)
library(quantmod)
library(plyr)
library(dplyr)
library(reshape2)
library(xts)
library(ggplot2)
library(stringr)

######################################################################
###############             Purpose                ###################
######################################################################
### This code implements the models developed by Diebold and Yilmaz, and Pesaran and Shin
### to estimate correlated impulse responses which are used to estimate forecast error variance
### decomposition. Using these estimates Diebold and Yilmaz (2011) develop a network visualization
### of directional variance between financial firms. Also, Diebold and Yilmaz (2009,2011) develop
### two indices of connectedness using two kinds of measures. The code below implements Generalized
### Variance Decomposition developed by Pesaran and Shin (1997), creates network visualizations
### and indices for connectedness
###


######################################################################
###############             Sources                ###################
######################################################################
### (Diebold & Yilmaz 2011)
### "On the Network Topology of Variance Decompositions: Measuring
### the Connectedness of Financial Firms". The URL for the white paper:
### http://www.ssc.upenn.edu/~fdiebold/papers/paper106/dieboldandyilmaz2011.pdf
###
### (Diebold & Yilmaz 2010)
### "Better to Give than to Receive: Predictive Directional Measurement of Volatility 
### Spillovers", 
### URL: http://www.ssc.upenn.edu/~fdiebold/papers/paper99/DirectionofSpillovers_Mar10.pdf
### 
### (Diebold and Yilmaz 2009)
### "Measuring Financial Asset Return and Volatility Spillovers, with Application to Global Equity Markets"
### URL: http://www.ssc.upenn.edu/~fdiebold/papers/paper75/DY2final.pdf
### 
### (Pesaran and Shin 1997)
### The model uses the Generalized Variance Decomposition method to calculate
### Forecast Error Variance Decomposition as described by Pesaran and Shin
### in "Generalized Impulse Response Analysis in Linear Multivariate Models"
### The URL for the whitepaper:
### http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.153.5596&rep=rep1&type=pdf
################################################################################

######################################################################
###############             Functions                #################
######################################################################

### the getSymbols function in quantmod retuns a list of xts objects.
### merge.all is used within in the Reduce function to combine all elements
### in the list into one xts object for time series analysis
merge.all <- function(x,y) { merge.xts(x,y,join="outer")}

### See Diebold and Yilmaz (2010) for the reasoning behind vol.fn
vol.fn <- function(history){log(100*sqrt(252*0.361*log(history[,3]/history[,2])^2))}

### the simple.ret function computes the simple returns on interday prices.
### The ClCl function in quantmod does this too but it
### eliminates the names of the columns which have been reordered by 
### the getSymbols function. Hence to track the column names
### and assure data integrity, the custom function is used instead.
simple.ret <- function(closing) {100*((closing[,4]-stats::lag(closing[,4]))/stats::lag(closing[,4]))}

### the directional.matrix function estimates the forecast error variance
### decomposition matrix using generalized impulse responses developed by
### Pesaran and Shin 1997. Diebold and Yilmaz use this matrix to develop
### a "Connectedness Table", see Table 2 in Diebold and Yilmaz 2011 and
### Table 3 in Diebold and Yilmaz 2010.

directional.matrix <- function(resid.var,ma.coef){
  ### ei and ej are identity matrices, or
  ### each column in ei and ej is a vector with unity at the i/jth element
  ### see equation 2.7 in Pesaran and Shin (1997)
  cov.error <- cov(resid.var)
  nobs <- dim(cov.error)[1]
  ei <- diag(nobs)
  ej <- diag(nobs)
  
  ### Initialize matrix of network connections
  Dg <- matrix(NA,ncol=nobs,nrow=nobs)
  
  ### Derive Forecast Error Variance Decomposition (FEVD) with 
  ### Generalized Variance Decomposition method
  for(j in 1:(ncol(ej))){
    ### the following steps calculate the denominator
    ### to estimate FEVD on page 3 of Pesaran and Shin (1997)
    ejj <- ej[,j]
    little.sigma <- 1/sqrt(cov.error[j,j])
    
    denom.step1 <- llply(ma.coef,function(x) t(ejj)%*%x%*%cov.error%*%t(x)%*%ejj)
    denom = Reduce('+',denom.step1)
    for(i in 1:(ncol(ej))){
      ### the following steps calculates equation 2.10 in Pesaran and Shin (1997)
      eii <- ei[,i]
      num.step1 <- llply(ma.coef,function(x) t(ejj)%*%x%*%cov.error%*%eii)
      num.squared <- llply(num.step1,function(x) x^2)
      num.sum <- Reduce('+',num.squared)
      num <- little.sigma*num.sum
      
      Dg[j,i] <- (num/denom)
      
    }
  }
  Dg
}

fevd.matrix <- function(resid.var,ma.coef){
  ### ei and ej are identity matrices, or
  ### each column in ei and ej is a vector with unity at the i/jth element
  ### see equation 2.7 in Pesaran and Shin (1997)
  Qinv <- t(chol(resid.var))
  cov.error <- cov(resid.var)
  
  nobs <- dim(Qinv)[1]
  ei <- diag(nobs)
  ej <- diag(nobs)
  
  ### Initialize matrix of network connections
  Dg <- matrix(NA,ncol=nobs,nrow=nobs)
  
  ### Derive Forecast Error Variance Decomposition (FEVD) with 
  ### Generalized Variance Decomposition method
  for(j in 1:(ncol(ej))){
    ### the following steps calculate the denominator
    ### to estimate FEVD on page 3 of Pesaran and Shin (1997)
    ejj <- ej[,j]
    
    denom.step1 <- llply(ma.coef,function(x) t(ejj)%*%x%*%cov.error%*%t(x)%*%ejj)
    denom = Reduce('+',denom.step1)
    for(i in 1:(ncol(ej))){
      ### the following steps calculates equation 2.10 in Pesaran and Shin (1997)
      eii <- ej[,i]
      num.step1 <- llply(ma.coef,function(x) t(ejj)%*%x%*%Qinv%*%eii)
      num.squared <- llply(num.step1,function(x) x^2)
      num <- Reduce('+',num.squared)
      
      Dg[j,i] <- (num/denom)
      
    }
  }
  Dg
}

### Normalize the columns forecast error variance decomposition specified
### on page 6 of Diebold and Yilmaz 2011
normalize <- function(amatrix){
  
  norm <- apply(amatrix,1,sum)
  
  amatrix/norm
}

### Create Connectedness Index as in Diebold and Yilmaz 2009 -
###
### Diebold/Yilmaz's notation is a bit sloppy when it comes to matrix operations.
### Their definition of the trace of (A0*A0), where A0 is matrix, is the sum
### of all the squared elements in A0 not just the diagonal. The functions get.trace
### and cond.sum perform this kind of calculation. Also, note that 
### the multiplication is elementwise and not a matrix product.

### get.trace calculates variance share for a particular variable
### i.e. it calculates the row sum for the directional network and includes
### the variable's impulse response on itself
get.trace <- function(x){ 
  A0 <- x%*%t(x)
  sum(A0)
}

### get.trace calculates the cross-variance share for a particular variable
### i.e. it calculates the row sum for the directional network and subtracts
### the variable's impulse response on itself

cond.sum <- function(mat){
  mat.sq <- mat%*%t(mat)
  low.tri <- mat.sq[lower.tri(mat.sq,diag=F)]
  up.tri <- mat.sq[upper.tri(mat.sq,diag=F)]
  new.mat <- sum(low.tri) + sum(up.tri)
  new.mat
}


sum.table <- function(mat,names,to,from){
  mat.df <- as.data.frame(mat)
  colnames(mat.df) <- names
  mat.df$From <- from
  
  to2 <- c(to,mean(to))
  
  mat.df <- rbind(mat.df,to2)
  rownames(mat.df) <- c(names,"To")
  mat.df
}

make.table <- function(mat,names,date.col){
  
  temp.mat <- mat
  
  diag(temp.mat) <- 0
  to <- apply(temp.mat,2,sum)
  from <- apply(temp.mat,1,sum)
  
  mat.df <- sum.table(mat,names,to,from)
  
  net <- to - from
  
  from <- as.data.frame(from,row.names = names)
  to <- as.data.frame(to,row.names = names)
  net <- as.data.frame(net,row.names = names)
  
  colnames(from) <- date.col
  colnames(to) <- date.col
  colnames(net) <- date.col
  
  list(table=mat.df,To=to,From=from,Net=net)
  
}

datesFn <- function(theData.xts,window=100,iter=10,lag=3){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  
  mySeq <- seq(1,nrow(theData.xts)-window-1,by=iter)
  min_iter = min(diff(index(theData.xts[mySeq])),na.rm=T)
  date_names <- seq.Date(from=min(index(theData.xts)),by=min_iter,length.out = length(mySeq))
  
  list(dates = date_names,step = min_iter)
}

vector_autoreg <- function(theData.xts,window=100,iter=10,ar_lag=3,ma_lag=10){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  
  mySeq <- seq(1,nrow(theData.xts)-window-1,by=iter)
  min_iter = min(diff(index(theData.xts[mySeq])),na.rm=T)
  date_names <- as.character(seq.Date(from=min(index(theData.xts)),by=min_iter,length.out = length(mySeq)))
  
  model_list <- vector("list",length=length(date_names))
  names(model_list) <- date_names
  ### initialize a list to hold index values
  fevd.list <- vector("list",length=length(date_names))
  ### names of the list come from the index of .xts data object
  names(fevd.list) <- date_names
  
  spillover.table.list <- vector("list",length=length(date_names))
  names(spillover.table.list) <- date_names
  
  vol.conn.index <- vector("list",length=length(date_names))
  ### names of the list come from the index of .xts data object
  names(vol.conn.index) <- date_names
  
  network.list <- vector("list",length=length(date_names))
  names(network.list) = date_names
  net.df.list <- network.list
  nvars <- ncol(theData.xts)
  
  for(d in 1:length(mySeq)){
    ### estimate VAR(p) and SVAR(p)
    var_vol_temp <- VAR(theData.xts[mySeq[d]:(mySeq[d]+window-1)],p=ar_lag,type="none")
    amat <- diag(nvars)
    amat[lower.tri(amat)] <- NA
    svar_vol_temp <- SVAR(var_vol_temp,estmethod = "direct",Amat=amat)
    
    model_list[[date_names[d]]] <- list(orthogonal=svar_vol_temp,
                                        non_orthogonal = var_vol_temp)
    ### get MA coefficients
    theta.temp <- Phi(var_vol_temp,nstep=ma_lag)
    svar.theta <- Phi(svar_vol_temp,nstep=ma_lag)
    ### convert array to listn.ahead
    theta.list <- alply(theta.temp,3)
    svar.theta.list <- alply(svar.theta,3)
    ################## SPILLOVER #########################
    svar_ecov <-svar_vol_temp$Sigma.U
    error.variance.total <- fevd.matrix(svar_ecov,svar.theta.list)*100
    spillover.table.list[[date_names[d]]] <- make.table(error.variance.total,colnames(theData.xts),date_names[d])
    
    
    Qinv <- t(chol(svar_ecov))
    irf <- llply(theta.list,function(x) x%*%Qinv)
    num = llply(irf,cond.sum)
    ### calculate own variance shares (i.e. do not subtract the diagonal of the matrix)
    den = llply(irf,get.trace)
    sum.num = sum(unlist(num))
    sum.den = sum(unlist(den))
    S = sum.num*100/sum.den
    fevd.list[[date_names[d]]] <- S
    ########################################################
    
    ################## CONNECTEDNESS #######################
    e.cov <- cov(residuals(var_vol_temp))
    ### create directional network table
    Sconn = directional.matrix(residuals(var_vol_temp),theta.list)
    ### to normalize or not to normalize
    D = normalize(Sconn)*100
    ### calculate total connectedness index level
    vol.conn.index[[date_names[d]]] = (sum(apply(D,1,sum)  - diag(D)))/ncol(D)
    ########################################################
    
    ################## Networks #######################
    rownames(D) <- colnames(theData.xts)
    colnames(D) <- rownames(D)
    
    df.list <- make.table(D,rownames(D),date_names[d])
    
    net.df.list[[date_names[d]]] <- df.list
    ### Create network Graph
    
    ### calculate net pairwise directional connectedness, Diebold pg 4.
    ### the calculation t(D) - D is the net pairwise equivalent of 
    ### the gross net = to - from, the calculation just above.
    abs.D = abs(D)
    net.mat = t(abs.D) - abs.D
    net.melt1 <- melt(net.mat)
    net.melt <- net.melt1 %>% filter(value != 0 )

    ### calculate each percentile of the net pairwise connectedness values
    ### and choose only the top 10%
    net.quant <- quantile(net.melt$value,prob=seq(0,1,by=0.01))
    new.net1 <- net.melt[net.melt$value > net.quant[90],]
    new.net <- new.net1[new.net1[,1] != new.net1[,2],]
    
    ### create igraph graph
    net.network <- graph.data.frame(new.net,direct=T)
    ### set graph nodes
    V(net.network)
    ### set edge colors
    E(net.network)$color <- ifelse(E(net.network)$value >= net.quant[99],"black",
                                   ifelse((E(net.network)$value < net.quant[99] & E(net.network)$value >= net.quant[95]),"red",
                                          ifelse((E(net.network)$value < net.quant[95] & E(net.network)$value >= net.quant[90]),"orange","grey")))
    
    ### set node size
    V(net.network)$size <- degree(net.network)/.5
    
    network.list[[date_names[d]]] <- net.network
    #####################################################################################
    remove(var_vol_temp,theta.temp,theta.list,e.cov,Sconn,
           num,den,sum.num,sum.den,S)
  }
  
  list(models = model_list,
       spillover_table = spillover.table.list,
       spill = fevd.list,
       conn=vol.conn.index,
       nets = network.list,
       net_df = net.df.list,
       dates = date_names)
  
}
