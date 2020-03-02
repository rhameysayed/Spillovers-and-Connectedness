library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinycssloaders)
library(shinyFiles)
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
vol.fn <- function(history){asinh(sqrt(252*0.361*log(history[,3]/lag(history[,2]))^2))}

### the simple.ret function computes the simple returns on interday prices.
### The ClCl function in quantmod does this too but it
### eliminates the names of the columns which have been reordered by 
### the getSymbols function. Hence to track the column names
### and assure data integrity, the custom function is used instead.
simple.ret <- function(closing) {100*((closing[,4]-lag(closing[,4]))/lag(closing[,4]))}
cumulative_return <- function(closing){100*((closing[,4]/as.numeric(closing[1,4]))-1)}
log_prices <- function(closing){log(closing[,4])}
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
    
    denom.step1 <- llply(ma.coef,function(x) t(ejj)%*%x%*%resid.var%*%t(x)%*%ejj)
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

datesFn <- function(theData.xts,window=100,iter=10){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  
  thedates <- index(theData.xts)
  start <- c()
  end <- c()
  maxdate <- max(thedates)
  t <- (length(thedates) - window)/iter
  i <- 1
  bigT <- t*iter
  while(i < bigT){
    
    start <- c(start,thedates[i])
    end <- c(end,thedates[i+window])
    i <- i + iter
  }
  
  date.df <- na.omit(data.frame(start=as.Date(start),end=as.Date(end)))
  
  end <- na.omit(end)
  end <- as.character(as.Date(end))
  
  
  list(dates = end,
       date.df = date.df)
}

vector_autoreg <- function(theData.xts,window=100,iter=10,ar_lag=3,ma_lag=10){
  ### create a sequence to use as cut offs for the rolling window
  ### estimation. iter determines how far to move the window forward
  dates <- datesFn(theData.xts,window=100,iter=10)
  date_names <- dates[["dates"]]
  date.df <- dates[["date.df"]]
  ndates <- length(date_names)
  alldates <- index(theData.xts)
  
  model_list <- vector("list",length=length(date_names))
  names(model_list) <- date_names
  centrality <- model_list
  ### initialize a list to hold index values
  fevd.list <- model_list
  ### names of the list come from the index of .xts data object
  spillover.table.list <- model_list
  vol.conn.index <- model_list
  ### names of the list come from the index of .xts data object
  network.list <- model_list
  net.df.list <- model_list
  nvars <- ncol(theData.xts)
  
  
  for(d in 1:ndates){
    ### estimate VAR(p) and SVAR(p)
    var_dates <- alldates[alldates >= date.df[d,"start"] & alldates <= date.df[d,"end"]]
    
    var.vol.temp <- VAR(theData.xts[var_dates,],p=ar_lag,type="none")
    amat <- diag(nvars)
    amat[lower.tri(amat)] <- NA
    svar.vol.temp <- SVAR(var.vol.temp,estmethod = "direct",Amat=amat)
    
    model_list[[date_names[d]]] <- list(orthogonal=svar.vol.temp,
                                        non_orthogonal = var.vol.temp)
    ### get MA coefficients
    theta.temp <- Phi(var.vol.temp,nstep=ma_lag)
    svar.theta <- Phi(svar.vol.temp,nstep=ma_lag)
    ### convert array to listn.ahead
    theta.list <- alply(theta.temp,3)
    svar.theta.list <- alply(svar.theta,3)
    ################## SPILLOVER #########################
    svar_ecov <-svar.vol.temp$Sigma.U
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
    e.cov <- cov(residuals(var.vol.temp))
    ### create directional network table
    Sconn = directional.matrix(residuals(var.vol.temp),theta.list)
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
    
    net.mat = t(D) - D
    net.melt1 <- melt(net.mat)
    net.melt <- net.melt1 %>% filter(value != 0 )
    colnames(net.melt)[3] <- "weight"
    
    ### calculate each percentile of the net pairwise connectedness values
    ### and choose only the top 10%
    net.quant <- quantile(net.melt$weight,prob=seq(0,1,by=0.01))
    new.net1 <- net.melt[net.melt$weight > net.quant[75],]
    new.net <- new.net1[new.net1[,1] != new.net1[,2],]
    
    ### create igraph graph
    
    net.network <- graph.data.frame(new.net,direct=T)
    E(net.network)$weight <- as.numeric(new.net$weight)
    ### set graph nodes
    #V(net.network)
    ### set edge colors
    E(net.network)$color <- ifelse(E(net.network)$weight >= net.quant[99],"black",
                                   ifelse(E(net.network)$weight < net.quant[99] & E(net.network)$weight >= net.quant[95],"red",
                                          ifelse(E(net.network)$weight < net.quant[95] & E(net.network)$weight >= net.quant[90],"orange","blue")))
    
    ### set node size
    V(net.network)$size <- degree(net.network)/.5
    
    network.list[[date_names[d]]] <- net.network
    #####################################################################################
    remove(var.vol.temp,theta.temp,theta.list,e.cov,Sconn,
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


ui <-dashboardPage(
   
   dashboardHeader(title = "Spillovers"),
   
   dashboardSidebar(
     sidebarMenu(
       
       selectInput("choice","Choose One:",choices=c("Returns","Volatility"),selected = "Volatility"),
       dateRangeInput("dateRange","Select Dates",
                      start = "2006-01-02",
                      end = "2016-11-11"),
       textInput("stocks","Ticker Input",value = "BK,BNS,TD,WFC,GS,JPM,MS,BAC,CS,DB,HSBC,C,AIG"),
       numericInput("window","Rolling Window (Days)",value=100),
       numericInput("windowDivider","Window Increment",value=10),
       numericInput("ar_lag","AR Lag (Days)",value=3),
       numericInput("ma_lag","MA Lag",value=10),
       actionButton("getIt","Submit")
       
     )
   ),
   dashboardBody(
    navbarPage("",
       tabPanel(title = "Index",
                fluidRow(
                  plotOutput("spilloverIndex",height="750px",hover="index_hover") %>% withSpinner(color="#0dc5c1")
                ),
                fluidRow(verbatimTextOutput("index_hover_data"))
                ),
       tabPanel(title = "Network Visual",
                fluidRow(uiOutput("net_dates")),
                fluidRow(
                  plotOutput("connectedness",height="800px")
                )),
       tabPanel(title = "Give/Receive Graphs",
                fluidRow(
                  selectInput("direction","Direction",choices = c("To","From","Net"),selected="To"),
                  plotOutput("grid_graph",height = "800px")  %>% withSpinner(color="#0dc5c1")
                )),
       tabPanel(title = "Give/Receive Table",
                fluidRow(uiOutput("table_dates")),
                
                fluidRow(column(5,tableOutput("connected_table")),
                         column(1,downloadButton("downloadConnectedness",label="Download Connectedness Table"),offset=3)),
                       
                fluidRow(column(5,tableOutput("spillover_table")),
                column(1,downloadButton("downloadSpillover",label="Download Spillover Table"),offset=3))
                ),
       tabPanel(title = "Pairwise Plot",
                fluidRow(uiOutput("sender_tickers"),uiOutput("receiver_tickers")),
                fluidRow(plotOutput("net_pairwise_plot",hover = "net_pairwise_hover")  %>% withSpinner(color="#0dc5c1")),
                fluidRow(verbatimTextOutput("net_pairwise_hover_data"))),
       tabPanel(title = "Data Viz",
                fluidRow(uiOutput("data_viz_checkbox")),
                fluidRow(plotOutput("log_price_plot")  %>% withSpinner(color="#0dc5c1")),
                fluidRow(plotOutput("cumulative_returns")  %>% withSpinner(color="#0dc5c1")),
                fluidRow(plotOutput("return_data_plot")  %>% withSpinner(color="#0dc5c1")),
                fluidRow(plotOutput("vol_data_plot")  %>% withSpinner(color="#0dc5c1"))
                ),
       tabPanel(title = "Returns Data",
                fluidRow(
                  column(5,tableOutput("return_table")),
                  column(1,downloadButton("downloadReturns",label="Download Returns Data"),offset=2))),
       tabPanel(title = "Volatility Data",
                fluidRow(
                  column(5,tableOutput("vol_table")),
                  column(1,downloadButton("downloadVolatility",label="Download Volatility Data"),offset=2)))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$sender_tickers <- renderUI({
    ticker_symbols <- str_split(input$stocks,",")[[1]]
    selectInput("sender","Sender",choices = ticker_symbols,selected = ticker_symbols[1])
  })
  output$receiver_tickers <- renderUI({
    ticker_symbols <- str_split(input$stocks,",")[[1]]
    selectInput("receiver","Receiver",choices = ticker_symbols,selected = ticker_symbols[2])
  })
  output$data_viz_checkbox <- renderUI({
    ticker_symbols <- str_split(input$stocks,",")[[1]]
    checkboxGroupInput("data_viz_check","Select Ticker to Show:",choices = ticker_symbols,selected = ticker_symbols[1],inline=TRUE)
  })
  
  ### read symbols from yahoo using getSymbols from quantmod package
  returns.data <- eventReactive(input$getIt,{
    tickers <- input$stocks
    tickers <- unlist(strsplit(tickers,split=","))
    tickers <- tickers[order(tickers)]
    
    
    dat <- new.env()
    getSymbols(tickers,
               src="yahoo",
               from=input$dateRange[1],
               to = input$dateRange[2],
               env=dat)
    
    
    logprices <- lapply(dat,log_prices)
    logprices <- Reduce(merge.all,logprices)
    
    
    cum_returns <- lapply(dat,cumulative_return)
    cum_returns <- Reduce(merge.all,cum_returns)
    
    rets <- lapply(dat,simple.ret)
    vol <- lapply(dat,vol.fn)
    ## combine each ticker's xts object into one xts object
    vol.data <- Reduce(merge.all,vol)
    vol.data <- na.omit(vol.data)
    ret.data <- Reduce(merge.all,rets)
    ret.data <- na.omit(ret.data)
    
    ret1.names <- colnames(vol.data)
    ## replace names with actual ticker symbols
    newNames <- sapply(as.character(ret1.names),function(x) substr(x,1,nchar(x)-4))
    colnames(ret.data) <- newNames
    colnames(vol.data) <- newNames
    colnames(cum_returns) <- newNames
    colnames(logprices) <- newNames
    
    list(log_price = logprices,
         returns = ret.data,
         cum_returns = cum_returns,
         return_vol = vol.data)
  })
  
  returns_data <- reactive({
    dat <- returns.data()
    dat[["returns"]]
  })
  
  volatility_data <- reactive({
    dat <- returns.data()
    dat[["return_vol"]]
  })
  
  output$return_table <- renderTable({
    returns_data()
  },rownames = TRUE)
  
  output$return_data_plot <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["returns"]][,input$data_viz_check],lwd=2,auto.grid=F,ylab="Return (%)",xlab="Time",
                     main="Daily Return (%)",lty=1,
                     legend.loc="topright")
  })
  
  output$vol_table <- renderTable({
    volatility_data()
  },rownames = TRUE)
  
  output$vol_data_plot <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["return_vol"]][,input$data_viz_check],lwd=2,auto.grid=F,ylab="Annualized Log Volatility",xlab="Time",
                     main="Log Volatility",lty=1,
                     legend.loc="topright")
  })
  
  output$log_price_plot <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["log_price"]][,input$data_viz_check],lwd=2,auto.grid=F,ylab="Log Prices",xlab="Time",
                     main="Log Prices",lty=1,
                     legend.loc="topright")
  })
  
  output$cumulative_returns <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["cum_returns"]][,input$data_viz_check],lwd=2,auto.grid=F,ylab="Cumulative Returns",xlab="Time",
                     main="Cumulative Returns (%)",lty=1,
                     legend.loc="topright")
  })
  
  theDates <- eventReactive(input$getIt,{
    dat <- returns.data()
    if(input$choice == "Returns"){
      dat = dat[["returns"]]
    }else{
      dat = dat[["return_vol"]]
    }
    datesFn(dat,window=input$window,iter=input$windowDivider)
  })
  
  output$net_dates <- renderUI({
    dates <- theDates()
    days <- dates[["dates"]]
    shinyWidgets::sliderTextInput(inputId = "net_dates", 
                                  label = "Date", 
                                  choices = days)
  })
  output$shock_dates <- renderUI({
    
    dates <- theDates()
    days <- dates[["dates"]]
    shinyWidgets::sliderTextInput(inputId = "shock_dates", 
                                  label = "Date", 
                                  choices = days)
  })
  output$table_dates <- renderUI({
    
    dates <- theDates()
    days <- dates[["dates"]]
    shinyWidgets::sliderTextInput(inputId = "table_dates", 
                                  label = "Date", 
                                  choices = days)
  })
  
  
  model_run <- eventReactive(input$getIt,{
    
    dat <- returns.data()
    if(input$choice == "Returns"){
      dat = dat[["returns"]]
    }else{
      dat = dat[["return_vol"]]
    }
    
    vector_autoreg(dat,window=input$window,iter=input$windowDivider,ar_lag=input$ar_lag,ma_lag=input$ma_lag)
  })
  
  connected_grid_graphs <- reactive({
    net.df.list <- model_run()[["net_df"]]
    spill.df.list <- model_run()[["spillover_table"]]
    
    net.cols <- names(net.df.list)
    
    selected.net.df.list <- llply(net.df.list,function(df){
      df[[input$direction]]
      })
    
    selected.net.df <- as.data.frame(do.call(rbind,lapply(selected.net.df.list, function(x) t(data.frame(x)))))
    selected.net.df <- selected.net.df %>% mutate(Date = as.POSIXct(net.cols,format="%Y-%m-%d"))
    
    spill.cols <- names(spill.df.list)
    
    selected.spill.df.list <- llply(spill.df.list,function(df){
      df[[input$direction]]
    })
    
    selected.spill.df <- as.data.frame(do.call(rbind,lapply(selected.spill.df.list, function(x) t(data.frame(x)))))
    selected.spill.df <- selected.spill.df %>% mutate(Date = as.POSIXct(spill.cols,format="%Y-%m-%d"))
    
    net.melt <- melt(selected.net.df,id="Date",value.name = "Connectedness")
    spill.melt <- melt(selected.spill.df,id="Date",value.name = "Spillover")
    merge(net.melt,spill.melt,by=c("Date","variable"))
  })
  
  output$grid_graph <- renderPlot({
    tbl <- connected_grid_graphs()
    ggplot(tbl) + 
      geom_line(aes(x=Date,y=Connectedness),colour="black") +
      geom_line(aes(x=Date,y=Spillover),colour="red") +
      facet_wrap(~variable,scales="free_y") +
      theme(legend.position = "bottom")
    
  })
  
  conn_table <- reactive({
    dat <- model_run()
    tables <- dat[["net_df"]]
    tables[[as.character(input$table_dates)]][["table"]]
  })
  output$connected_table <- renderTable({
    conn_table()
  },rownames = TRUE)
  
  spill_table <- reactive({
    dat <- model_run()
    tables <- dat[["spillover_table"]]
    tables[[as.character(input$table_dates)]][["table"]]
  })
  
  output$spillover_table <- renderTable({
    spill_table()
  },rownames = TRUE)
  
  
  
  output$net_pairwise_plot <- renderPlot({
    sender <- input$sender
    receiver <- input$receiver
    
    dat <- model_run()
    net_stuff <- dat[["net_df"]]
    tables <- llply(net_stuff,function(stuff) stuff[["table"]])
    
    pairwise_data <- ldply(tables,function(table){
      sender_index <- which(colnames(table) == sender)
      receiver_index <- which(colnames(table) == receiver)
      
      sender_shock <- table[receiver_index,sender_index]
      receiver_shock <- table[sender_index,receiver_index]
      
      net_shock <- sender_shock - receiver_shock
    })
    colnames(pairwise_data) <- c("Time","Connectedness")
    pairwise_data <- pairwise_data %>% mutate(Time = as.POSIXct(Time,format="%Y-%m-%d"))
    
    ggplot(pairwise_data,aes(x=Time,y=Connectedness,group=1)) + 
      geom_line() +
      geom_hline(yintercept=0,color="red",linetype="dashed")+
      ggtitle(paste(sender,"to", receiver)) + 
      xlab("Time") + 
      ylab("Connectedness")
  })
  
  output$net_pairwise_hover_data <- renderText({
    
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("Date=", as.Date(as.POSIXct(e$x,origin="1970-01-01",format="%Y-%m-%d",tz="UTC")), " Index=", round(e$y,2), "\n")
    }
    xy_str(input$net_pairwise_hover)
  })
  
  net.net <- reactive({
    nets <- model_run()
    nets[["nets"]]
  })
  
  output$connectedness <- renderPlot({
    netnets <- net.net()
    nets <- netnets[[as.character(input$net_dates)]]
    plot.igraph(nets,layout=layout.circle(nets),
         main="Firm Connectedness \n (intra day return volatility)",
         xlab = "black: 1st percentile \n red: 5th percentile \n orange: 10th percentile \n blue: 25th percentile \n node size: number of edeges connected to the node")
    
  })
   
   indexes <- reactive({
     
     index <- model_run()
     
     ############## CONNECTEDNESS ################
     vol.conn.index.out <- index[["conn"]]
     ### unlist output data
     vol.conn.index.df <- data.frame(unlist(vol.conn.index.out))
     colnames(vol.conn.index.df) = c("Index")
     ### .xts object
     vol.conn.index.df$Date <- as.POSIXct(rownames(vol.conn.index.df),format="%Y-%m-%d")
     
     ############## SPILLOVER ###################
     fevd.list.out <- index[["spill"]]
     
     fevd.df <- data.frame(unlist(fevd.list.out))
     colnames(fevd.df) = c("Index")
     fevd.df$Date <- as.POSIXct(rownames(fevd.df),format="%Y-%m-%d")
     
     ### compare the two index measures
     indice = merge(vol.conn.index.df,fevd.df,by="Date")
     colnames(indice) = c("Date","Connectedness","Spillover")
     indice
   })
   
   output$spilloverIndex <- renderPlot({
     dat <- indexes()
     dat <- xts(dat[,-1],order.by=dat[,1])
     
     ggplot(dat) + 
       geom_line(aes(x=Index,y=Connectedness,colour="Connectedness"),size=1.2) +
       geom_line(aes(x=Index,y=Spillover,colour="Spillover"),size=1.2) + theme_bw() +
       theme(legend.position = "bottom",legend.title = element_blank(),text = element_text(size=20)) + 
       xlab("Time") + ylab("Connectedness/Spillover Index") +
       scale_color_manual(values = c("black","red"))
   })
   
   output$index_hover_data <- renderText({
    
       xy_str <- function(e) {
         if(is.null(e)) return("NULL\n")
         paste0("Date=", as.Date(as.POSIXct(e$x,origin="1970-01-01",format="%Y-%m-%d",tz="UTC")), " Index=", round(e$y,2), "\n")
       }
       xy_str(input$index_hover)
   })
   
   output$downloadConnectedness <- downloadHandler(
     filename = function() {
       paste("Connectedness-",as.character(input$table_dates), ".csv", sep = "")
     },
     content = function(file){
       write.csv(conn_table(), file,row.names = TRUE)
     }
   )
   
   output$downloadSpillover <- downloadHandler(
     filename = function() {
       paste("Spillover-",as.character(input$table_dates), ".csv", sep = "")
     },
     content = function(file) {
       
       write.csv(spill_table(), file,row.names = TRUE)
     }
   )
   
   output$downloadReturns <- downloadHandler(
     filename = function() {
       "Returns.csv"
     },
     content = function(file){
       write.csv(returns_data(), file,row.names = TRUE)
     }
   )
   
   output$downloadVolatility <- downloadHandler(
     filename = function() {
       "Volatility.csv"
     },
     content = function(file) {
       
       write.csv(volatility_data(), file,row.names = TRUE)
     }
   )
}

# Run the application 
shinyApp(ui = ui, server = server)
