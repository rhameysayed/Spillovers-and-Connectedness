library(shiny)
library(shinydashboard)

#source("C:\\Users\\Rhamey\\Desktop\\financeR\\spillover_connectedness_functions.R")


ui <-dashboardPage(
   
   dashboardHeader(title = "Spillovers"),
   
   dashboardSidebar(
     sidebarMenu(
       checkboxGroupInput("choice","Choose One:",choices=c("Returns","Volatility"),selected = "Volatility"),
       dateRangeInput("dateRange","Select Dates",
                      start = "2001-01-02",
                      end = Sys.Date()),
       textInput("stocks","Ticker Input",value = "BK,BNS,TD,WFC,GS,JPM,MS,BAC,CS,DB,HSBC,C,AIG"),
       numericInput("window","Rolling Window (Days)",value=100),
       numericInput("windowDivider","Forecast Horizion",value=10),
       numericInput("ar_lag","AR Lag (Days)",value=3),
       numericInput("ma_lag","MA Lag",value=10),
       actionButton("getIt","Done")
     )
   ),
   dashboardBody(
    navbarPage("",
       tabPanel(title = "Index",
                fluidRow(
                  plotOutput("spilloverIndex",height="800px")
                )),
       tabPanel(title = "Network Visual",
                fluidRow(uiOutput("net_dates")),
                fluidRow(
                  plotOutput("connectedness",height="800px")
                )),
       tabPanel(title = "Give/Receive Graphs",
                fluidRow(
                  selectInput("direction","Direction",choices = c("To","From","Net"),selected="To"),
                  plotOutput("grid_graph",height = "800px")
                )),
       tabPanel(title = "Give/Receive Table",
                fluidRow(uiOutput("table_dates")),
                fluidRow(tableOutput("connected_table")),
                fluidRow(tableOutput("spillover_table"))),
       tabPanel(title = "Pairwise Plot",
                fluidRow(uiOutput("sender_tickers"),uiOutput("receiver_tickers"))),
       tabPanel(title = "Returns Data",
                fluidRow(plotOutput("return_data_plot")),
                fluidRow(tableOutput("return_table"))),
       tabPanel(title = "Volatility Data",
                fluidRow(plotOutput("vol_data_plot")),
                fluidRow(tableOutput("vol_table")))
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
  
  ### read symbols from yahoo using getSymbols from quantmod package
  returns.data <- eventReactive(input$getIt,{
    tickers <- input$stocks
    tickers <- unlist(strsplit(tickers,split=","))
    tickers <- tickers[order(tickers)]
    
    withProgress(message="Getting Data",value=0,
                 {
                   for(i in 1:length(tickers)){
                    incProgress(1/length(tickers))
                    Sys.sleep(0.75)
                  }
                 })
    
    dat <- new.env()
    getSymbols(tickers,
               src="yahoo",
               from=input$dateRange[1],
               to = input$dateRange[2],
               env=dat)
    
    rets <- lapply(dat,simple.ret)
    vol <- lapply(dat,vol.fn)
    ## combine each ticker's xts object into one xts object
    vol.data <- Reduce(merge.all,vol)
    ret.data <- Reduce(merge.all,rets)
    ret.data <- na.omit(ret.data)
    
    ret1.names <- colnames(vol.data)
    ## replace names with actual ticker symbols
    newNames <- sapply(as.character(ret1.names),function(x) substr(x,1,nchar(x)-4))
    colnames(ret.data) <- newNames
    colnames(vol.data) <- newNames
    
    list(returns = ret.data,
         return_vol = vol.data)
  })
  
  output$return_table <- renderTable({
    dat <- returns.data()
    dat[["returns"]]
  },rownames = TRUE)
  
  output$return_data_plot <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["returns"]],lwd=2,auto.grid=F,ylab="Return (%)",xlab="Time",
                     main="Daily Return Percentage",lty=1,
                     legend.loc="topright")
  })
  
  output$vol_table <- renderTable({
    dat <- returns.data()
    dat[["return_vol"]]
  },rownames = TRUE)
  
  output$vol_data_plot <- renderPlot({
    dat <- returns.data()
    chart.TimeSeries(dat[["return_vol"]],lwd=2,auto.grid=F,ylab="Annualized Volatility Percentage",xlab="Time",
                     main="Daily Return Percentage",lty=1,
                     legend.loc="topright")
  })
  
  theDates <- eventReactive(input$getIt,{
    dat <- returns.data()
    if(input$choice == "Returns"){
      dat = dat[["returns"]]
    }else{
      dat = dat[["return_vol"]]
    }
    datesFn(dat,window=input$window,iter=input$windowDivider,lag=input$ar_lag)
  })
  
  output$net_dates <- renderUI({
    dates <- theDates()
    sliderInput("net_dates","Date",
                min=min(dates[["dates"]]),
                max=max(dates[["dates"]]),
                step=dates[["step"]],
                value=min(dates[["dates"]]),
                width = "800px")
  })
  output$shock_dates <- renderUI({
    dates <- theDates()
    sliderInput("shock_dates","Date",
                min=min(dates[["dates"]]),
                max=max(dates[["dates"]]),
                step=dates[["step"]],
                value=min(dates[["dates"]]),
                width = "800px")
  })
  output$table_dates <- renderUI({
    dates <- theDates()
    sliderInput("table_dates","Date",
                min=min(dates[["dates"]]),
                max=max(dates[["dates"]]),
                step=dates[["step"]],
                value=min(dates[["dates"]]),
                width = "800px")
  })
  
  
  model_run <- eventReactive(input$getIt,{
    # req(input$window)
    # req(input$windowDivider)
    # req(input$ar_lag)
    # req(input$ma_lag)
    
    dat <- returns.data()
    if(input$choice == "Returns"){
      dat = dat[["returns"]]
    }else{
      dat = dat[["return_vol"]]
    }
    
    withProgress(message='Running VAR',detail="this will take a minute..",value=0,{
      steps <- round((nrow(dat)-input$window)/input$windowDivider,0)
      for(i in 1:steps){
        incProgress(1/steps)
        Sys.sleep(0.1)
      }
    })
    
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
  
  output$connected_table <- renderTable({
    dat <- model_run()
    tables <- dat[["net_df"]]
    tables[[as.character(input$table_dates)]][["table"]]
  },rownames = TRUE)
  
  output$spillover_table <- renderTable({
    dat <- model_run()
    tables <- dat[["spillover_table"]]
    tables[[as.character(input$table_dates)]][["table"]]
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
  
  net.net <- reactive({
    nets <- model_run()
    nets[["nets"]]
  })
  
  output$connectedness <- renderPlot({
    netnets <- net.net()
    nets <- netnets[[as.character(input$net_dates)]]
    plot.igraph(nets,layout=layout.circle(nets),
         main="Firm Connectedness \n (intra day return volatility)",
         xlab = "black: 1st percentile \n red: 5th percentile \n orange: 10th percentile \n node size: number of edeges connected to the node")
    
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
     vol.conn.index.xts <- xts(vol.conn.index.df[,-2],order.by=vol.conn.index.df[,2])
     
     ############## SPILLOVER ###################
     fevd.list.out <- index[["spill"]]
     
     fevd.df <- data.frame(unlist(fevd.list.out))
     colnames(fevd.df) = c("Index")
     fevd.df$Date <- as.POSIXct(rownames(fevd.df),format="%Y-%m-%d")
     fevd.xts <- xts(fevd.df[,-2],order.by=fevd.df[,2])
     
     ### compare the two index measures
     indice = merge.all(vol.conn.index.xts,fevd.xts)
     colnames(indice) = c("Connectedness","Spillover")
     indice
   })
   
   output$spilloverIndex <- renderPlot({
     chart.TimeSeries(indexes(),lwd=2,auto.grid=F,ylab="Index",xlab="Time",
                      main="Comparing Spillover Index levels",lty=1,
                      legend.loc="topright")
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
