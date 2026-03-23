#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(plotly)
library(magrittr)
library(markdown)
source("k0_estimation.R")
source("simulate_heavy_tails.R")

models = c('Pareto','Burr','Frechet','T');

D=read.table('calcium.csv')
Xcal=D[[1]]
res <- reactiveValues();
#data<-data.frame("x"=Xcal,"i"=c(1:length(Xcal)));

# Turn off warnings 

options(warn=-1)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("The Trimmed Hill Estimator: Robust and adaptive tail inference "),
  
  fluidRow(
    
    column(12,
           includeMarkdown("Description.md")
    )
  ),
  
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      radioButtons("type","Category",c("Simulated Data","Uploaded Data","Calcium Data"),selected = 'Calcium Data'),
      numericInput("alpha","Type I error (for Adaptive Trimmed Hill)",min = 0, max = 1,value =0.05,step = 0.01),
      numericInput("k0","Value of K0 (for Trimmed Hill plot)",min = 0, max = 500,value =7),
      sliderInput("k","Value of k:",min = 1,max = 500,value = 50),
      # Only show this panel if the plot type is a histogram
      conditionalPanel(
        condition = "input.type == 'Uploaded Data'",
        fileInput("file1", "Choose CSV File (single column)",multiple = FALSE,
                  accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
        numericInput("d","Dithering Parameter",min = 0, max = 1000,value =0.01),
       # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        # Input: Select separator ----
        radioButtons("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ","),
        # Input: Select quotes ----
       
        radioButtons("quote", "Quote",choices = c(None = "","Double Quote" = '"',"Single Quote" = "'"),selected = '"')
      ),
    
      conditionalPanel(
        condition = "input.type == 'Simulated Data'",
        radioButtons("mod","Data",choiceValues =c("Pareto","Burr","Frechet","T"),choiceNames =
                       c("Pareto","Burr","Frechet","T")),
        sliderInput("n","Sample size",min=2,value=2000,max=10000),
        numericInput("k0_outliers","Number of top outliers",min = 0, max = 500,value =7),
        numericInput("L","Outlier exponent parameter L",min = -10, max = 50,value =2,step = 0.2),
        numericInput("xi","Tail index xi ",min = 0.0001, max = 50,value =0.5,step = 0.01)
      )
   ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data", plotlyOutput("dataPlot")),
        tabPanel("Trim Diagnostics", plotlyOutput("TrimDiagnostic")),
        tabPanel("Trimmed Hill plot", plotlyOutput("trimmedHillPlot"))
      )
    )))
# Define server logic required to draw a histogram
server <- function(session,input, output) {
  
  options(shiny.maxRequestSize=10*1024^2)
  

  dataInput <- reactive({
    if (input$type=='Calcium Data'){x <-Xcal+runif(length(Xcal), min =-0.01,max=0.01 ); res$name='Calcium Data'}
    if (input$type=='Simulated Data'){
      if (input$mod=='Pareto'){x <- simulate_pareto(input$n,1/input$xi,1); res$name='Pareto Data'};
      if (input$mod=='Burr'){x <- simulate_burr(input$n,1,0.5,1/input$xi);res$name= 'Burr Data'};
      if (input$mod=='Frechet'){x <- simulate_frec(input$n,1/input$xi); res$name='Frechet Data'};
      if (input$mod=='T'){x <- simulate_T(input$n,1/input$xi); res$name='|T| Data'}
      out =sort.int(x,index.return = T,decreasing = T);
      idx <- out$ix;
      if (input$k0_outliers>0){
        x[idx[1:(input$k0_outliers+1)]]  <- (1+x[idx[1:(input$k0_outliers+1)]] - x[idx[input$k0_outliers+1]])^input$L-1 + x[idx[input$k0_outliers+1]];
      }
    }
    if (input$type=='Uploaded Data'){
      req(input$file1)
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          x <- c(as.matrix(read.csv(input$file1$datapath,header = input$header,sep = input$sep,quote = input$quote)))
          x=x/min(x)+runif(length(x),min =-input$d,max=input$d)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      res$name='Uploaded Data'
    }
     res$x <-x;
    res$i <- c(1:length(x));
    res$data <- data.frame("x"=c(1:length(x)),"y"=x);
  })
  
 
  observe({
    dataInput();
    val <- length(res$x)
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session,"k", min = 0, max = val-1,value = 50)
  })
  
  observe({
    val <- input$k
    # Control the value, min, max, and step.
    # Step size is 2 when input value is even; 1 when value is odd.
    updateSliderInput(session,"k0", min = 0, max = val-2,value = 10)
    updateSliderInput(session,"k0_outliers", min = 0, max = val-2,value = 10)
  })
  
  
  
  
  output$dataPlot <- renderPlotly({
    dataInput();
    data = res$data;
    
    sorted_data = sort(data$y,index.return = TRUE, decreasing = TRUE);
    est = compute_ADAP(data$y,input$k,level = input$alpha);
    
    if (est$k0>0){
      idx_outliers = sorted_data$ix[c(1:est$k0)];
    } else {
      idx_outliers = c(); 
    } 
    
    plot_ly(data,x=~data$x,y=~data$y,type="scatter",mode="points",name="Data")%>%
      add_trace(x=idx_outliers,y=data$y[idx_outliers],name="Estimated outliers",mode="markers",marker=list("size"=12,color="red"))%>%
      layout(title=paste0( "Data Set:",res$name),
             xaxis=list("title"="i"),
             yaxis=list("title"="x(i)"));
  })
  
  output$trimmedHillPlot <- renderPlotly({
    dataInput();
    
    x <- res$data$y;
    
    #
    # Change the name of the functions:
    #  (1) "censor_hill_plot" to "trimmed_hill_plot"
    #  (2) "trimmed_hill_plot" to "trim_diagnostic_plot"
    #
    censor_hill = censor_hill_plot(x,k0=input$k0,k=input$k);
    classic_hill = censor_hill_plot(x,k0=0,k=input$k);
    est = compute_ADAP(x,input$k,level = input$alpha);
    ##adap_hill = censor_hill_plot(x,k0=est$k0,k=input$k);
    
    data  = data.frame("x"=censor_hill$k0+c(1:length(censor_hill$xi)),
                       "y"=censor_hill$xi,
                       "y_wrong" = censor_hill$xi_wrong,
                       "sigma"=censor_hill$sigma);
    
    data_hill = data.frame("x"=classic_hill$k0+c(1:length(classic_hill$xi)),
                           "y"=classic_hill$xi,"sigma"=classic_hill$sigma);
    ##data_adap_hill = data.frame("x"=adap_hill$k0+c(1:length(adap_hill$xi)),
    ##                            "y"=adap_hill$xi);
    
    #
    # Plotting
    #
    
    name_censored = paste("Trimmed Hill plot k0=",input$k0,sep="");
    name_Hill = paste("Classic Hill plot (k0=0) based on all the data");
    ##name_adap_Hill = paste("Adaptive Trimmed Hill plot k0=",est$k0,sep="");
    name_Hill_wrong = paste("Biased Hill plot (Classic Hill based on missing top-k0=",input$k0,")",sep="");
    
    if (input$type=='Simulated Data'){
      plot_ly(data, x=~data$x,y=~data$y,type="scatter",mode="lines+markers", name=name_censored,
              error_y = ~list(type = "data", array = data$sigma, opacity=0.2, color="red"),
              line=list("color"="red")) %>%
        add_trace(data_hill, x=~data_hill$x,y=~data_hill$y, name=name_Hill, error_y = NULL, mode="lines+markers",
                  line=list("color"="black","dash"="dot")) %>%
        ##add_trace(data_adap_hill, x=~data_adap_hill$x,y=~data_adap_hill$y, name=name_adap_Hill, error_y = NULL, mode="lines+markers",
        ##         line=list("color"="black","dash"="dash")) %>%
        add_trace(data,x=~data$x,y=~data$y_wrong,name=name_Hill_wrong,mode="lines",error_y=NULL,
                  line=list("color"="black","dash"="dash")) %>%
        add_trace(x=as.numeric(c(min(data_hill$x),max(data_hill$x))),y=as.numeric(c(input$xi,input$xi)),mode="lines",error_y=NULL, 
                  showlegend = T, name="True xi",
                  type="scatter",line= list("color"="black","width"=1,"dash"="dash"),opacity=0.5)%>%
        layout(title = paste("The Trimmed Hill plot vs its biased alternatives. Data Set: ",input$mod,sep=""),
               xaxis=list("title"="k"),
               yaxis=list("title"="xi"));
    } else{
      plot_ly(data, x=~data$x,y=~data$y,type="scatter",mode="lines+markers", name=name_censored,
              error_y = ~list(type = "data", array = data$sigma, opacity=0.2, color="red"),
              line=list("color"="red")) %>%
        add_trace(data_hill, x=~data_hill$x,y=~data_hill$y, name=name_Hill, error_y = NULL, mode="lines+markers",
                  line=list("color"="black","dash"="dot")) %>%
        ##add_trace(data_adap_hill, x=~data_adap_hill$x,y=~data_adap_hill$y, name=name_adap_Hill, error_y = NULL, mode="lines+markers",
        ##          line=list("color"="black","dash"="dash")) %>%
        add_trace(data,x=~data$x,y=~data$y_wrong,name=name_Hill_wrong,mode="lines",error_y=NULL,
                  line=list("color"="black","dash"="dash")) %>%
        layout(title = paste("The Trimmed Hill plot vs its biased alternatives. Data Set: ",input$mod,sep=""),
               xaxis=list("title"="k"),
               yaxis=list("title"="xi"));
    }
  })
  
  output$TrimDiagnostic <- renderPlotly({
    
    dataInput();
    
    x <- res$data$y;
    #
    # Computations: trimmed Hill and ADAP
    #
    trim_hill = trim_hill_plot(x,input$k);
    est = compute_ADAP(x,input$k,level = input$alpha);
    
    data = data.frame("x"=trim_hill$k0,"y"=trim_hill$xi,"sigma"=trim_hill$sigma);
    
    
    if (est$k0>0){
      idx_k0 = c(0:(est$k0-1));
    } else {
      idx_k0 = c(); 
    } 
    
    if ( input$type=='Simulated Data' ){
      plot_ly(data, x=~data$x,y=~data$y,type="scatter",mode="lines+markers", 
              error_y = ~list(type = "data", array = data$sigma, opacity=0.2, color="red"), showlegend = F)%>%
        add_trace(x = as.numeric(est$k0,est$k0),y=c(0,max(data$y)),line = list("color"="red",dash="dash"),name="k0 (estimated)")%>%
        add_trace(x=as.numeric(c(min(data$x),max(data$x))),y=as.numeric(c(input$xi,input$xi)),error_y=NULL, showlegend = T, name="True xi",
                  type="scatter",line= list("color"="black","width"=1,"dash"="dash"),mode="lines",opacity=0.5) %>% 
        add_trace(x=idx_k0,y=data$y[1+idx_k0],name="Estimated outliers",mode="markers",marker=list("size"=12,color="red"))%>%
        layout(title = paste("Trim Diagnostic plot: k0(estimated) = ",est$k0,",  xi=", format(est$xi,digits=4),sep=""),
               xaxis=list("title"="k0"),
               yaxis=list("title"="xi"));
    } else {
      plot_ly(data, x=~data$x,y=~data$y,type="scatter",mode="lines+markers", 
              error_y = ~list(type = "data", array = data$sigma, opacity=0.2, color="red"), showlegend = F)%>%
        add_trace(x = as.numeric(est$k0,est$k0),y=c(0,max(data$y)),line = list("color"="red",dash="dash"),name="k0 (estimeted)")%>%
        add_trace(x=idx_k0,y=data$y[1+idx_k0],name="Estimated outliers",mode="markers",marker=list("size"=12,color="red"))%>%
        layout(title = paste("Trim Diagnostic plot: k0(estimated) = ",est$k0,",  xi=", format(est$xi,digits=4),sep=""),
               xaxis=list("title"="k0"),
               yaxis=list("title"="xi"));
    }
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

