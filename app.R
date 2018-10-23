#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(mclust)
library(stats)
library(DT)
library(MASS)
library(MethComp)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        div(style=paste0("display: inline-block;vertical-align:top; width: 200px;"),sliderInput("cls",
                     "Number of clusters:",
                     min = 3,
                     max = 15,
                     value = 8)),
        div(style=paste0("display: inline-block;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: inline-block;vertical-align:top; width: 200px;"),sliderInput("clsid",
                     "Display cluster:",
                     min = 1,
                     max = 3,
                     value = 1)),
        fileInput(
         "lasfile1",
         "Загрузить LAS файл:",
          accept = '.las',
          buttonLabel = "Открыть..."
        ),
        dataTableOutput('datasetDT'),
        div(style=paste0("display: table-cell;vertical-align:top; width: 120px;"),selectInput("selectInp", "Исходный",choices = NULL)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),selectInput("selectOut", "Прогнозный",choices = NULL)),
        div(style=paste0("display: vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),checkboxInput('logInp','log(1+x)',value = T)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),checkboxInput('logOut','log(1+x)', value = T)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),actionButton("runCalc",label = 'Рассчитать',icon = icon("play-circle"))),
        #div(style=paste0("display: table-cell;vertical-align:top; width: 10px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 120px;"),downloadButton('downloadDataset', 'Save dataset'))
        #div(style=paste0("display: vertical-align:top; width: 20px;"),HTML("<br>")),
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot",height = "800px")
      )
   )
)
myReactives <- reactiveValues()

plotError <- function (message = "Error!!!") {
  par(mfrow=c(1,1),pch=16)
  plot(0,0,t="l");text(0,0,paste("ОШИБКА:",message),col = "red")
}

getFormStr <- function(dlm) {
  if(is.null(dlm) || is.null(dlm$coefficients))
    return('')
  form = '' 
  for(term in names(dlm$coefficients)) { 
    if(term!='(Intercept)') form = paste(form,term,'*')
    form = paste(form,prettyNum(dlm$coefficients[[term]]),'+')
    }
  return(form)
}

extract_wellname <- function(oneline) {
  
  temp <- stringr::str_trim(oneline) # remove initial and final spaces
  temp <- gsub("\\s+", " ", temp) # replace multiple spaces with 1 space
  temp <- sub("WELL *\\.", "", temp) # remove initial "WELL ." with any number of spaces between WELL and .
  temp <- sub("\\: *\\w* *\\w* *\\w*", "", temp) # remove final ": WELL" with any number of spaces between : and WELL
  temp <- stringr::str_trim(temp) #remove initial and final spaces
  
  return(temp)
}

read_las <- function(filename,
                     lasnull="-999.250000") {
  
  if(!file.exists(filename)) stop("File '", filename, "' does not exist!")
  #browser()
  lasfile <- file(filename, open="r")
  
  ### 1. read header of the las file (LAS files have all the headers within first 100 lines)
  headerlines <- readLines(con = lasfile, n = 100L, ok = TRUE, skipNul = TRUE)
  
  # 1.a get the well name
  #pattern <- paste("WELL *\\.", "WELL\\.", sep = "|")
  oneline <- headerlines[grep("WELL *\\.", headerlines)]
  #wellname <- stringr::word(gsub("\\s+"," ", stringr::str_trim(oneline)),2)[1]
  wellname <- extract_wellname(oneline)
  
  # 1.b get the names of the logs
  pattern_c <- paste("~C" ,"~Curve", "~CURVE", sep = "|")
  aa <- grep(pattern_c, headerlines)
  bb <- grep("~", headerlines)[which(grep("~", headerlines) %in% aa) + 1] # look for the next line with ~ after ~C
  
  # count number of log curves, discounting comment lines starting with "#"
  ncurves <- bb - aa - 1 - length(grep("#", headerlines[aa:bb]))
  logname <- seq(0, 0, ncurves)
  curveslines <- subset(headerlines[(aa):(bb-1)], !grepl("#", headerlines[(aa):(bb-1)])) # remove comment lines
  for (i in 1:ncurves) {
    oneline <- curveslines[1 + i]
    logname[i] <- stringr::word(gsub("\\s+", " ", stringr::str_trim(oneline)), 1)[1]
  }
  
  # 1.c get the first line with data
  pattern_a <- paste("~A" ,"~Ascii", "~ASCII", sep = "|")
  dataline <- as.numeric(grep(pattern_a, headerlines) + 1)
  
  ### 2. read the log data from the line after ~Ascii
  temp <- utils::read.table(lasfile,
                            header = FALSE,
                            na.strings = lasnull,
                            skip = dataline-1,
                            col.names = logname)
  
  close(lasfile)
  
  #temp
  #temp <- temp %>%
  #  dplyr::mutate(WELL = wellname) %>%
  #  dplyr::select(WELL, 1:ncurves)
  
  return(temp)
}

tls <- function(A, b){
  
  n <- ncol(A)
  C <- cbind(A, b)
  
  V <- svd(C)$v
  VAB <- V[1:n, (n+1):ncol(V)]
  VBB <- V[(n+1):nrow(V), (n+1):ncol(V)]
  return(-VAB/VBB)
}

tls_pca <- function (x, y) {
  v <- prcomp(cbind(x,y))$rotation
  beta <- v[2,1]/v[1,1]
  b0=mean(y)-beta*mean(x)
  return(c(Intercept=b0,Gradient=beta))
}

getTranspRamp <- function(color = NULL ,grades = 0.5,max_alpha=0.5) {
  if(is.null(colors)) return(NULL)
  
  #browser()
  colors=rep(color,grades)
  grades=seq(0,max_alpha,length.out = grades)
  colors=sapply(colors, col2rgb)/255
  colors = sapply(seq_along(colors[1,]),  
                  function(x) rgb(colors[1,x], colors[2,x], colors[3,x], alpha=grades[x]))  
  #dbgmes("transpal=",colors)
  return(colors)
}
options(shiny.maxRequestSize = 500 * 1024 ^ 2)
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  showModDial <- function(message = "Ожидайте...") {
    showModal(modalDialog( list(imageOutput("cycle", 
                                            width = "50px",
                                            height = "50px",
                                            inline = T),
                                message),
                           title = "Ожидайте...", footer = modalButton("Закрыть")))
  }
  
  output$downloadDataset <- downloadHandler(
     filename = function() {
       paste0('Mclust_regr',Sys.Date(), '.txt')
     },
     contentType = '.txt',
     content = function (fname) {
       #browser()
       if(is.null(myReactives$mcres) || length(myReactives$mcres)<2) return(NULL)
       save = try( expr = write.table(x = myReactives$mcres,
                                      file = fname,
                                      sep = '    \t',
                                      dec = '.',
                                      row.names = F,
                                      col.names = T) 
                   , TRUE)
       if(class(save)=="try-error")
         showNotification(ui = "Ошибка при сохранении файла",
                          type = "error")      
     }
   )
   
   #### Precess data ####
   observe({
     if(is.null(myReactives$data) || all(!(dim(myReactives$data)>2))) {
       #plotError("No Data")
       print("no data")
       return(NULL)
     }
     if(is.null(input$datasetDT_rows_selected) || length(input$datasetDT_rows_selected)<2) {
       updateSelectInput(session = session,"selectInp",choices = character(0))
       updateSelectInput(session = session,"selectOut",choices = character(0))
       print("selected < 2")
       return(NULL)
     } 
     #browser() 
     parList = names(myReactives$data)
     parList = parList[parList %in% parList[input$datasetDT_rows_selected]]
     if(length(parList)<1) {
       print("no pars")
       return(NULL)
     }
     selInp=input$selectInp
     selOut=input$selectOut
     if(!(selInp %in% parList)) selInp=NULL
     if(!(selOut %in% parList)) selOut=NULL
     updateSelectInput(session = session,"selectInp",choices = parList,selected = selInp)
     updateSelectInput(session = session,"selectOut",choices = parList,selected = selOut)
     #q = myReactives$data[,c('DEPT','DT','dCALI','SP','SN')]
     #print(parList)
     q = myReactives$data[,c('DEPT',parList)]
     q = q[complete.cases(q),]
     #browser()
     #q$dCALI = log10(q$dCALI)
     qlog = log(q)
     qlog$DEPT=q$DEPT

     qlog = qlog[complete.cases(qlog),]
     #plot(qlog[-1])
     #cls = 4
     #dat = mcres[mcres$VVV12 == cls,]
     #plot(dat$DT~dat$SN)
     #dlm = lm(dat$DT~dat$SN)
     #browser()
     #par(new=T,mfrow=c(2,2),pch=16)
     if(dim(qlog)[1]<3) return(NULL)
     myReactives$qlog = qlog
     
     #myReactives$mbic = Mclust((myReactives$qlog[-1]),G=seq(2,3,3),modelNames = c('EEI','VVV','VEE'))
   })
  
   #### Calc clusters ####
   #observe({
   observeEvent(input$runCalc,{
     if(is.null(myReactives$qlog) || length(myReactives$qlog)<3
        || is.null(input$selectInp) || is.null(input$selectOut)
        || input$selectInp=="" || input$selectOut=="") {
       #plotError("No Data")
       return(NULL)
     }
     #browser()
     showModDial()
     qlog = myReactives$qlog
     if(input$logInp) qlog[,input$selectInp] = log1p(qlog[,input$selectInp])
     if(input$logOut) qlog[,input$selectOut] = log1p(qlog[,input$selectOut])
     # qlog$DT = log1p(qlog$DT)
     # qlog$SN = log1p(qlog$SN)
     
     qlog = qlog[complete.cases(qlog),]
     if(dim(qlog)[1]<3) return(NULL)
     print(dim(qlog))
     #browser();
     myReactives$mc = Mclust((qlog[-1]),G=input$cls,modelNames = 'VVV')
     if(is.null(myReactives$mc)) {
       removeModal()
       return(NULL)
     }
     myReactives$mcres = cbind(qlog,VVV12 = predict(myReactives$mc)$classification)
     updateSliderInput(session = session,inputId = 'clsid',min = 1,max=input$cls,value = min(input$clsid,input$cls))
     removeModal()
   })
   #### Plot data ####
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
     if(is.null(myReactives$mc) ||
        length(myReactives$data)<2
        || is.null(input$selectInp) || is.null(input$selectOut)
        || input$selectInp=="" || input$selectOut=="") {
       plotError("No Data")
       return(NULL)
     }
     #split.screen(c(2,2))
     #screen(1);par(pch=16)
     par(mfrow = c(2,2))
     
     mbic=myReactives$mbic
     #par(mfg = c(2,1))
     #plot(mbic,what = "BIC")
     mbic=myReactives$mc
     #par(mfg = c(2,2),pch=16)
     #par(pch=16)
     #plot(mbic,what = "classification")
     mcres = myReactives$mcres
     #mcres = cbind(qlog,VVV12 = predict(mbic)$classification)
     dlms=list()
     pal=mclust.options("classPlotColors")
     #screen(1);par(pch=16)
     par(mfg = c(1,1),pch=16)
     #plot all classes ####
     plot(mcres[,input$selectOut]~mcres[,input$selectInp],col=pal[mcres$VVV12]);
     classes= sort(unique(mcres$VVV12))
     legend("topright",inset = c(0.01,0.01),legend = classes ,col = pal[classes],pch = 16)
     #browser()
     #for(cls in unique(mcres$VVV12)) 
     cls = input$clsid
     {
       dat = mcres[mcres$VVV12 == cls,]
       lm=glm(dat[,input$selectOut]~dat[,input$selectInp])#,family = Gamma)
       ab=lm$coefficients
       predicted=ab[1]+ab[2]*dat[,input$selectInp]
       lm$cc=ccf(dat[,input$selectOut],predicted,lag.max = 0,plot = F)$acf[[1]]
       
       lm.inv=glm(dat[,input$selectInp]~dat[,input$selectOut])#,family = Gamma)
       ab.inv=c(-lm.inv$coefficients[1]/lm.inv$coefficients[2],b=1/lm.inv$coefficients[2])
       predicted.inv=ab.inv[1]+ab.inv[2]*dat[,input$selectInp]
       lm.inv$cc=ccf(dat[,input$selectOut],predicted.inv,lag.max = 0,plot = F)$acf[[1]]
       
       #lm.odr=Deming(dat[,input$selectInp],dat[,input$selectOut])
       #ab.odr=c(lm.odr[1],lm.odr[2])
       #predicted.odr=ab.odr[1]+ab.odr[2]*dat[,input$selectInp]
       #odr.cc=ccf(dat[,input$selectOut],predicted.odr,lag.max = 0,plot = F)$acf[[1]]
       #browser()
       
       ab.pca=tls_pca(dat[,input$selectInp],dat[,input$selectOut])
       predicted.pca=ab.pca[1]+ab.pca[2]*dat[,input$selectInp]
       pca.cc=ccf(dat[,input$selectOut],predicted.pca,lag.max = 0,plot = F)$acf[[1]]
       
       
       
       if(sd(lm$residuals) < 1.5*sd(lm.inv$residuals)) {
         #predicted=predicted
         cc=lm$cc
       } else {
         predicted=predicted.inv
         cc=lm.inv$cc
       }
       cc=pca.cc
       predicted=predicted.pca
       
       print(paste("LM    :",cls,"="));print(ab);print(mad(lm$residuals))
       #print(capture.output(lm))
       print(paste("LM.INV:",cls,"="));print(ab.inv);print(mad(lm.inv$residuals))
       #print(capture.output(lm.inv))
       print(paste("LM.ODR:",cls,'='));print(ab.pca);print(mad(predicted.pca-dat[,input$selectOut]))
     }
     #names(dlms) <- unique(mcres$VVV12)

     #plot single class ####
     mcres_2=mcres[mcres$VVV12 == input$clsid,]
     k <- kde2d(y = mcres_2[,input$selectOut],x = mcres_2[,input$selectInp])
     
     #screen(2);par(pch=16)
     par(mfg = c(1,2),pch=16)
     #lm=dlms[[paste(input$clsid)]]
     title = paste("CC for class",input$clsid,'=',prettyNum(cc))
     plot(mcres[,input$selectOut]~mcres[,input$selectInp],col='grey80',main = getFormStr(lm));
     #browser()
     colors=getTranspRamp(pal[input$clsid],grades = 5,max_alpha = 0.7)
     .filled.contour(x = k$x,y = k$y,z = k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=colors)

     #points(y=mcres[,input$selectOut],x=mcres[,input$selectInp],col='grey80');
     abline(ab,col=pal[input$clsid])
     
     par(lwd=2,lty="dashed")
     abline(ab.pca,col=pal[input$clsid])
     
     par(lwd=4,lty="dotted")
     abline(ab.inv,col=pal[input$clsid])
     
     par(lwd=1,lty="solid")
     #points(y=mcres_2[,input$selectOut],x=mcres_2[,input$selectInp],col=pal[mcres_2$VVV12]);
     #contour(k,nlevels=5,col='gray50',lwd = 2,add=T,drawlabels = F)
     #print(title)
     
     #plot LM residuals ####
     #screen(3);par(pch=16)
     par(mfg = c(2,1),pch=16)
     plot(lm.inv,which=1)
     #browser()
     
     #plot LM xplot ####
     #screen(4);par(pch=16)
     par(mfg = c(2,2),pch=16)
     if(input$logInp) measured=exp(exp(predicted)-1)
     else measured=exp(predicted)
     if(input$logOut) predicted=exp(exp(mcres_2[,input$selectOut])-1)
     else predicted = exp(mcres_2[,input$selectOut])
     plot(predicted~measured,col='grey30',main=title)
     #abline(lm(mcres_2$SN~predict(lm)),col='black')
     #split.screen(c(2,2))
     #screen(1)
     #plot(dat$DT,predict(dlm),main=ccf(dat$DT,predict(dlm),lag.max = 0,plot = F)$acf[[1]])
     
   })
  
   #### Load data ####
   observeEvent(input$lasfile1, {
     #data = read_table2("D:/Libya - Pyroclastic Study/to_send/LAS/SIRT_GF44_Subonsh_Geo100_A7-100_Logs.txt", 
     #             col_types = cols(X11 = col_skip()))
     data=read_las(input$lasfile1$datapath,"-999.25")
     #browser()
     #data = read_table2(input$lasfile1$datapath, 
     #             col_types = cols(X11 = col_skip()))
     #data = data[complete.cases(data.frame(apply(data,c(1,2),FUN = function(x) {if(x==-999.25) NA else x}))),]
     #mbic = Mclust(q,G = 9);plot(mbic)
     #mcres = cbind(q,predict(mbic)$classification)
     #browser()
     data = data.frame(apply(data,c(1,2),FUN = function(x) {if(x==-999.25) return(NA) else return(x)}))
     data$dBH = data$CALI
     data$dBH[data$DEPT<3427] = 13+3/8
     data$dBH[data$DEPT>=3427 & data$DEPT < 10877] = 9+5/8
     data$dBH[data$DEPT>=10877 & data$DEPT < 13567] = 7
     data$dBH[data$DEPT>=13567 & data$DEPT < 14663] = 5
     data$dCALI = data$CALI - data$dBH
     data$dCALI[data$dCALI<0]=NA
     #browser()
     myReactives$data=data
   })
   #### Plot Table ####
   output$datasetDT <- renderDataTable({
     data = myReactives$data
     if(is.null(data) || length(data)<2) return(NULL)
     #browser()
     logList = data.frame(names = names(data), 
                          min = sapply(names(data), FUN = function(x) min(data[, x], na.rm = T)), 
                          max = sapply(names(data), FUN =function(x) max(data[, x], na.rm = T)),
                          cover = sapply(names(data), FUN = function(x) 
                            round(100*length(data$DEPT[!is.na(data[,x])])/length(data$DEPT)))
                          )
     #cover = sapply(names(data), FUN = function(x) 
     #round(100*(max(data$DEPT[!is.na(data[,x])])-min(data$DEPT[!is.na(data[,x])]))/(max(data$DEPT)-min(data$DEPT))))
     myReactives$logList = logList
     datatable(myReactives$logList,
               rownames = F,
               class = 'compact',
               options = list(
                 pagingType = "simple",
                 paging = FALSE,
                 searching = FALSE,
                 ColumnRender = prettyNum,
                 scrollY = "400px",
                 scrollCollapse = TRUE,
                 pageLength = 15)
               ) %>%
       formatStyle(
         'cover',
        color = styleInterval(seq(0,100,34), c('white', 'red', 'blue','green'))
#         backgroundSize = '100% 10%',
#         backgroundRepeat = 'no-repeat',
#         backgroundPosition = 'center'
       )
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

