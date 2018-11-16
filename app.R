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
library(e1071)

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
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),actionButton("runCalc",label = 'Рассчитать',icon = icon("play-circle"))),
        #div(style=paste0("display: table-cell;vertical-align:top; width: 10px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 120px;"),downloadButton('downloadDataset', 'Save dataset')),
        div(style=paste0("display: vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 120px;"),selectInput("selectInp", "Исходный",choices = NULL)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),selectInput("selectOut", "Прогнозный",choices = NULL)),
        div(style=paste0("display: vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),checkboxInput('logInp','log(1+x)',value = T)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),checkboxInput('logOut','log(1+x)', value = T)),
        div(style=paste0("display: table-cell;vertical-align:top; width: 20px;"),HTML("<br>")),
        div(style=paste0("display: table-cell;vertical-align:top; width: 150px;"),checkboxInput('logAll','log(x)', value = T)),
        div(style=paste0("display: vertical-align:top; width: 20px;"),HTML("<br>"))
      ),
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot",height = "800px"),
         selectInput("selectXP", "Прогнозный",choices = c('glm','inv','avg','pca','lqs'))
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
  headerlines <- readLines(con = lasfile, n = 1000L, ok = TRUE, skipNul = TRUE)
  
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
  return(c(Intercept = b0,Gradient = beta))
}

descaleAB <- function (ab,dat_scaled,inout) {
  K=ab[2]
  B=ab[1]
  Kx=attributes(dat_scaled)$`scaled:scale`[inout[1]]
  Bx=attributes(dat_scaled)$`scaled:center`[inout[1]]
  Ky=attributes(dat_scaled)$`scaled:scale`[inout[2]]
  By=attributes(dat_scaled)$`scaled:center`[inout[2]]
  #descaled coefficient
  KK=Ky*K/Kx
  BB=By+Ky*B-K*Bx*Ky/Kx
  ab=c(Intercept=BB,Gradient=KK)
  return(ab)
}

getModel <- function(dat,inout=c(1,2),type='glm',scale=TRUE, logpOut=FALSE, logOut) {
  dat_scaled=scale(dat)
  measured=dat[,inout[2]]
  #dat_scaled[,inout[1]] = -sinh(dat_scaled[,inout[1]])
  #browser()
  if(type=='pca') ab=tls_pca(dat_scaled[,inout[1]],dat_scaled[,inout[2]])
  else if(type=='glm') ab=glm(dat_scaled[,inout[2]]~dat_scaled[,inout[1]])$coefficients
  else if(type=='lqs') ab=lqs(dat_scaled[,inout[2]]~dat_scaled[,inout[1]])$coefficients
  else if(type=='rlm') ab=rlm(dat_scaled[,inout[2]]~dat_scaled[,inout[1]])$coefficients
  else if(type=='inv') {
    ab=glm(dat_scaled[,inout[1]]~dat_scaled[,inout[2]])$coefficients
    ab=c(Intercept = -ab[1]/ab[2],Gradient = 1/ab[2])
  }
  if(type =='avg') {
    ab1=tls_pca(dat_scaled[,inout[1]],dat_scaled[,inout[2]])
    ab2=glm(dat_scaled[,inout[2]]~dat_scaled[,inout[1]])$coefficients
    ab3=glm(dat_scaled[,inout[1]]~dat_scaled[,inout[2]])$coefficients
    ab3=c(Intercept = -ab3[1]/ab3[2],Gradient = 1/ab3[2])
    ab=apply(rbind(ab2,ab3),2,FUN=mean)
  }
  #descale ab
  ab = descaleAB(ab = ab,dat_scaled = dat_scaled,inout = inout)
  
  predicted=ab[1]+ab[2]*dat[,inout[1]]  # no scaling
  #predicted=asinh(-predicted)
  resid=predicted-measured
  #estimate prob.densities
  denm=density(measured,from = min(measured),to = max(measured))
  denp=density(predicted,from = min(measured),to = max(measured))
  denr=density(resid,from = min(measured),to = max(measured))
  # estimating QC factors
  qcf0=sd(resid)
  qcf1=denr$bw
  #browser()
  #qcf2=sd(abs(denp$y-denm$y))#/max(denm$y))
  qcf2=abs(sum((denp$y-denm$y)**2))
  qcf3=abs(lm(resid~measured)$coefficients[2])
  qcf=qcf2
  #print(c(qcf0,qcf1,qcf))
  #predicted.pca=(predicted.pca-B)/K    #descale
  cc=ccf(measured,predicted,lag.max = 0,plot = F)$acf[[1]]
  
  #if(logpOut) predicted=exp(exp(predicted)-1)
  #if(logOut) predicted=exp(predicted)
  
  return(list(ab=ab,predicted=predicted,resid=resid,cc=cc,qcf=qcf))
}

setPaletteTransp <- function(colors = NULL ,alpha = 0.5) {
  if(is.null(colors)) return(NULL)
  
  colors = apply(sapply(colors, col2rgb)/255, 2, 
                 function(x) rgb(x[1], x[2], x[3], alpha=alpha)) 
  #dbgmes("transpal=",colors)
  return(colors)
}

getTranspRamp <- function(color = NULL ,grades = 0.5,alphaRange=c(0.1,0.5)) {
  if(is.null(colors)) return(NULL)
  
  #browser()
  colors=rep(color,grades)
  grades=seq(alphaRange[1],alphaRange[2],length.out = grades)
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
                                      row.names = T,
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
     print(input$datasetDT_rows_selected)
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
     #if(selOut==selInp) {}
     updateSelectInput(session = session,"selectInp",choices = parList,selected = selInp)
     updateSelectInput(session = session,"selectOut",choices = parList,selected = selOut)
     #q = myReactives$data[,c('DEPT','DT','dCALI','SP','SN')]
     #print(parList)
     q = myReactives$data[,c(parList)]
     rownames(q) <- as.character(myReactives$data$DEPT)
     #colnames(q)[1]='DPH0'
     qlog = q[complete.cases(q),]
     #browser()
     #q$dCALI = log10(q$dCALI)
     #qlog$DPH0=q$DPH0
     if(input$logAll) {
       qlog = log(q)
       qlog = qlog[complete.cases(qlog),]
       #print(capture.output(summary(qlog)))
     }
     # if(input$logInp) qlog[,input$selectInp] = log1p(q[,input$selectInp])
     # if(input$logOut) qlog[,input$selectOut] = log1p(q[,input$selectOut])
     # qlog = qlog[complete.cases(qlog),]
     
     print(c("qlog=",colnames(qlog)))
     #plot(qlog[-1])
     #cls = 4
     #dat = mcres[mcres$VVV12 == cls,]
     #plot(dat$DT~dat$SN)
     #dlm = lm(dat$DT~dat$SN)
     #browser()
     #par(new=T,mfrow=c(2,2),pch=16)
     if(dim(qlog)[1]<2) return(NULL)
     myReactives$qlog = qlog
     
     #myReactives$mbic = Mclust((myReactives$qlog[-1]),G=seq(2,3,3),modelNames = c('EEI','VVV','VEE'))
   })
  
   #### Calc clusters ####
   #observe({
   observeEvent(input$runCalc,{
     if(is.null(myReactives$qlog) || length(myReactives$qlog)<2
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
     qlog = qlog[complete.cases(qlog),]

     if(dim(qlog)[1]<2) return(NULL)
     print(dim(qlog))
     #browser();
     myReactives$mc = Mclust(qlog,G=input$cls,modelNames = 'VVV')
     if(is.null(myReactives$mc)) {
       removeModal()
       return(NULL)
     }
     #print(c("qlog=",colnames(qlog)))
     myReactives$mcres = cbind(qlog,VVV12 = predict(myReactives$mc)$classification)
     updateSliderInput(session = session,inputId = 'clsid',min = 1,max=input$cls,value = min(input$clsid,input$cls))
     removeModal()
   })
   #### Plot data ####
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
     #print(input$selectInp)
     #print(input$selectOut)
     #print(colnames(myReactives$data))
     #print(colnames(myReactives$mc$data))
     #print(myReactives$mc$parametrs)
     #print(all(colnames(myReactives$mc$data) %in% colnames(myReactives$data)))
     if(is.null(myReactives$mc) ||
        length(myReactives$data)<2
        || is.null(input$selectInp) || is.null(input$selectOut)
        || input$selectInp=="" || input$selectOut==""
        || !(input$selectInp %in% colnames(myReactives$mc$data) && input$selectOut %in% colnames(myReactives$mc$data)) ) {
       plotError("Выборка изменена, либо расчет не выполнен")
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
     # if(input$logInp) mcres[,input$selectInp] = log1p(mcres[,input$selectInp])
     # if(input$logOut) mcres[,input$selectOut] = log1p(mcres[,input$selectOut])
     # mcres = mcres[complete.cases(mcres),]
     #mcres = cbind(qlog,VVV12 = predict(mbic)$classification)
     dlms=list()
     pal=mclust.options("classPlotColors")
     #screen(1);par(pch=16)
     par(mfg = c(1,1),pch=16)
     #plot all classes ####
     #browser()
     plot(mcres[,input$selectOut]~mcres[,input$selectInp],col=pal[mcres$VVV12]);
     classes= sort(unique(mcres$VVV12))
     legend("topright",inset = c(0.01,0.01),legend = classes ,col = pal[classes],pch = 16)
     #browser( )
     #for(cls in unique(mcres$VVV12)) 
     cls = input$clsid
     {
       dat = mcres[mcres$VVV12 == cls,]
       mod.pca = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'pca',scale=T,logpOut = input$logOut,logOut = input$logAll)
       mod.glm = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'glm',scale=T,logpOut = input$logOut,logOut = input$logAll)
       mod.inv = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'inv',scale=T,logpOut = input$logOut,logOut = input$logAll)
       mod.avg = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'avg',scale=T,logpOut = input$logOut,logOut = input$logAll)
       mod.lqs = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'lqs',scale=T,logpOut = input$logOut,logOut = input$logAll)
       mod.rlm = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'rlm',scale=T,logpOut = input$logOut,logOut = input$logAll)
       #mod.svm = getModel(dat=dat,inout = c(input$selectInp,input$selectOut),type = 'svm',scale=T)

       measured=dat[,input$selectOut]
       #if(input$logInp) measured=exp(exp(measured)-1)
       #else measured=exp(measured)
       #if(input$logAll) measured=exp(measured)
       print('measured range')
       print(range(measured))
       
       lms<-data.frame(QCF=c(mod.glm$qcf,mod.inv$qcf,mod.avg$qcf,mod.pca$qcf,mod.lqs$qcf,mod.rlm$qcf),
                  CC=c(mod.glm$cc,mod.inv$cc,mod.avg$cc,mod.pca$cc,mod.lqs$cc,mod.rlm$cc),
                  rbind(mod.glm$ab,mod.inv$ab,mod.avg$ab,mod.pca$ab,mod.lqs$ab,mod.rlm$ab))
       rownames(lms) <- c('LM.glm','LM.inv','LM.avg','LM.pca','LM.lqs','LM.rlm')
     }
     #browser()
     lms <- lms[with(lms,order(QCF)),]
     print(lms)
     #plot single class ####
     mcres_2=mcres[mcres$VVV12 == input$clsid,]
     
     #screen(2);par(pch=16)
     par(mfg = c(1,2),pch=16)
     #lm=dlms[[paste(input$clsid)]]
     plot(mcres[,input$selectOut]~mcres[,input$selectInp],col='grey80',
          main = paste0(rownames(lms)[1],'\'s QCF: ',prettyNum(lms$QCF[1])));
     #browser()
     range = rbind(range(mcres[,input$selectOut]),range(mcres[,input$selectInp]))
     colors=getTranspRamp(pal[input$clsid],grades = 5,alphaRange = c(0.0,0.8))
     k <- kde2d(y = mcres_2[,input$selectOut],x = mcres_2[,input$selectInp])#,lims = c(range[1,],range[2,]))
     .filled.contour(x = k$x,y = k$y,z = k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=colors)
     #contour(k,col=pal[mcres_2$VVV12],add=TRUE)

     getTxtLoc <- function (ab,ranges, margin = 0.05) {
       # get the plot aspect ratio as 'asp'
       w <- diff(par("usr")[1:2]) / par("pin")[1] # plot units per inch horizontal axis
       h <- diff(par("usr")[3:4]) / par("pin")[2] # plot units per inch vertical axis
       asp <- w/h
       txt.y=ranges[1,2]-(ranges[1,2]-ranges[1,1])*margin
       txt.x=(txt.y-ab[1])/ab[2]
       txt.rot=180/pi*atan(ab[2]*asp)
       if(txt.x>ranges[2,2]) {
         txt.x=ranges[2,2]-(ranges[2,2]-ranges[2,1])*margin*asp
         txt.y=txt.x*ab[2]+ab[1]
       } else if(txt.x<ranges[2,1]) {
         txt.x=ranges[2,1]+(ranges[2,2]-ranges[2,1])*margin*asp
         txt.y=txt.x*ab[2]+ab[1]
       }
       return(list(x=txt.x,y=txt.y,srt=txt.rot))
     }
     
     
     #plot glm line
     abline(mod.glm$ab,col=pal[input$clsid])
     txt=getTxtLoc(mod.glm$ab,range)
     text(txt$x,txt$y,"glm",srt=txt$srt,pos=3,col=pal[input$clsid])
     #plot inverse glm line
     abline(mod.inv$ab,col=pal[input$clsid])
     txt=getTxtLoc(mod.inv$ab,range)
     text(txt$x,txt$y,"inv",srt=txt$srt,pos=3,col=pal[input$clsid])
     #plot PCA (ortho) line
     par(lwd=2,lty="dotted")
     abline(mod.pca$ab,col=pal[input$clsid])
     txt=getTxtLoc(mod.pca$ab,range)
     text(txt$x,txt$y,"pca",srt=txt$srt,pos=3,col=pal[input$clsid])
     #plot glm-inv averaged  line
     par(lwd=2,lty="longdash") 
     abline(mod.avg$ab,col=pal[input$clsid])
     txt=getTxtLoc(mod.avg$ab,range)
     text(txt$x,txt$y,"avg",srt=txt$srt,pos=3,col=pal[input$clsid])

     #plot glm-inv averaged  line
     par(lwd=2,lty="dashed") 
     abline(mod.rlm$ab,col=pal[input$clsid])
     txt=getTxtLoc(mod.rlm$ab,range)
     text(txt$x,txt$y,"rlm",srt=txt$srt,pos=3,col=pal[input$clsid])
     
     
     par(lwd=1,lty="solid")
     #points(y=mcres_2[,input$selectOut],x=mcres_2[,input$selectInp],col=pal[mcres_2$VVV12]);
     #points(y=mcres_2[,input$selectOut],x=mcres_2[,input$selectInp],col=colors[1]);
     #contour(k,nlevels=5,col=pal[mcres_2$VVV12],lwd = 1,add=T,drawlabels = F)
     #print(title)
     
     # select best fit 1) by sd(resid)?
     
     #plot LM residuals ####
     #screen(3);par(pch=16)
     par(mfg = c(2,1),pch=16,lwd = 1)
     #plot(lm.inv,which=3)
     cols <- c('black','blue','red','green')
     cols = setPaletteTransp(colors = cols,alpha = 0.1)
     names(cols) <- c('glm','inv','avg','pca')
       plot(mod.inv$resid ~ measured, col=cols['inv'])
     points(mod.avg$resid ~ measured, col=cols['avg'])
     points(mod.pca$resid ~ measured, col=cols['pca'])
     points(mod.glm$resid ~ measured, col=cols['glm'])
     points(mod.rlm$resid ~ measured, col=cols['rlm'])
     cols = setPaletteTransp(colors = cols,alpha = 1)
     resRange = c(range(measured),range(mod.inv$resid))
     k <- kde2d(x = measured,y = mod.glm$resid,lims = resRange);contour(k,nlevels=5,col=cols['glm'],add=T,drawlabels = F)
     k <- kde2d(x = measured,y = mod.inv$resid,lims = resRange);contour(k,nlevels=5,col=cols['inv'],add=T,drawlabels = F)
     k <- kde2d(x = measured,y = mod.avg$resid,lims = resRange);contour(k,nlevels=5,col=cols['avg'],add=T,drawlabels = F)
     k <- kde2d(x = measured,y = mod.pca$resid,lims = resRange);contour(k,nlevels=5,col=cols['pca'],add=T,drawlabels = F)
     k <- kde2d(x = measured,y = mod.rlm$resid,lims = resRange);contour(k,nlevels=5,col=cols['rlm'],add=T,drawlabels = F)
     abline(lm(mod.glm$resid~measured)$coefficients,col=cols['glm'])
     abline(lm(mod.inv$resid~measured)$coefficients,col=cols['inv'])
     abline(lm(mod.avg$resid~measured)$coefficients,col=cols['avg'])
     abline(lm(mod.pca$resid~measured)$coefficients,col=cols['pca'])
     abline(lm(mod.rlm$resid~measured)$coefficients,col=cols['rlm'])
     #abline(c(0,0),col='red')
     legend("topright",inset = c(0.01,0.01),
            legend = names(cols) ,
            col = cols,pch = 16)
     #browser()
     
     ab=mod.inv$ab
     cc=mod.inv$cc
     predicted=mod.inv$predicted
     
     mod=switch(input$selectXP,
               glm=mod.glm,
               inv=mod.inv,
               avg=mod.avg,
               pca=mod.pca,
               lqs=mod.lqs)
     
     ab=mod$ab
     cc=mod$cc
     predicted=mod$predicted
     #plot LM xplot ####
     #screen(4);par(pch=16)
     par(mfg = c(2,2),pch=16)
     #browser()
     ranges<-range(measured,predicted)
     title = paste("CC for class",input$clsid,'=',prettyNum(cc))
     cols = setPaletteTransp(colors = cols,alpha = 0.2)
     plot(data.frame(measured,measured),col='red',main=title,xlim=ranges,ylim=ranges)
     #plot(mod.glm$predicted~measured,col=cols['glm'],main=title,xlim=ranges,ylim=ranges)
     #points(mod.inv$predicted~measured,col=cols['inv'],main=title,xlim=ranges,ylim=ranges)
     #points(mod.avg$predicted~measured,col=cols['avg'],main=title,xlim=ranges,ylim=ranges)
     #points(mod.pca$predicted~measured,col=cols['pca'],main=title,xlim=ranges,ylim=ranges)
     #points(mod.rlm$predicted~measured,col=cols['rlm'],main=title,xlim=ranges,ylim=ranges)
     #plot(mod.lqs$predicted~measured,col=cols['lqs'],main=title,xlim=ranges,ylim=ranges)
     #plot(predicted~measured,col='grey30',main=title,xlim=ranges,ylim=ranges)
     
     abline(c(0,1),lty='dashed',lwd=2,col='red')
     # draw 2d densities with lines
     #k <- kde2d(x = measured,y = mod.glm$predicted,lims = rep(ranges,2));contour(k,nlevels=5,col=cols['glm'],add=T,drawlabels = F)
     #k <- kde2d(x = measured,y = mod.inv$predicted,lims = rep(ranges,2));contour(k,nlevels=5,col=cols['inv'],add=T,drawlabels = F)
     #k <- kde2d(x = measured,y = mod.avg$predicted,lims = rep(ranges,2));contour(k,nlevels=5,col=cols['avg'],add=T,drawlabels = F)
     #k <- kde2d(x = measured,y = mod.pca$predicted,lims = rep(ranges,2));contour(k,nlevels=5,col=cols['pca'],add=T,drawlabels = F)
     #k <- kde2d(x = measured,y = mod.rlm$predicted,lims = rep(ranges,2));contour(k,nlevels=5,col=cols['rlm'],add=T,drawlabels = F)
     # draw 2d densities with filled polys ####
     col.glm=getTranspRamp(cols['glm'],grades = 5,alphaRange = c(0.0,0.3))
     col.inv=getTranspRamp(cols['inv'],grades = 5,alphaRange = c(0.0,0.3))
     col.pca=getTranspRamp(cols['pca'],grades = 5,alphaRange = c(0.0,0.3))
     col.avg=getTranspRamp(cols['avg'],grades = 5,alphaRange = c(0.0,0.3))
     k <- kde2d(x = measured,y = mod.glm$predicted,lims = rep(ranges,2));.filled.contour(x=k$x,y=k$y,z=k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=col.glm)
     k <- kde2d(x = measured,y = mod.inv$predicted,lims = rep(ranges,2));.filled.contour(x=k$x,y=k$y,z=k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=col.inv)
     k <- kde2d(x = measured,y = mod.avg$predicted,lims = rep(ranges,2));.filled.contour(x=k$x,y=k$y,z=k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=col.avg)
     k <- kde2d(x = measured,y = mod.pca$predicted,lims = rep(ranges,2));.filled.contour(x=k$x,y=k$y,z=k$z,levels=seq(min(k$z),max(k$z),length.out = 5),col=col.pca)
#     k <- kde2d(x = measured,y = mod.rlm$predicted,lims = rep(ranges,2));.filled.contour(x=k$x,y=k$y,z=k$z,nlevels=5,col=col.glm,add=T,drawlabels = F)
     # draw correlation lines for each model ####
     cols = setPaletteTransp(colors = cols,alpha = 1)
     abline(lm(mod.glm$predicted~measured)$coefficients,col=cols['glm'])
     abline(lm(mod.inv$predicted~measured)$coefficients,col=cols['inv'])
     abline(lm(mod.avg$predicted~measured)$coefficients,col=cols['avg'])
     abline(lm(mod.pca$predicted~measured)$coefficients,col=cols['pca'])
     #abline(lm(mod.rlm$predicted~measured)$coefficients,col=cols['rlm'])
     #### draw 1d densities with lines ####
     scaler=(ranges[2]-ranges[1])*0.5
     km=density(measured);lines(x = km$x,y=(km$y-km$y)*scaler/max(km$y)+mean(ranges),col='magenta',lwd=3,lty='dashed')
     k<-density(mod.glm$predicted,from=ranges[1],to=ranges[2]);lines(x = k$x,y=(k$y-km$y)*scaler/max(km$y)+mean(ranges),col=cols['glm'],lwd=2)
     k<-density(mod.inv$predicted,from=ranges[1],to=ranges[2]);lines(x = k$x,y=(k$y-km$y)*scaler/max(km$y)+mean(ranges),col=cols['inv'],lwd=2)
     k<-density(mod.pca$predicted,from=ranges[1],to=ranges[2]);lines(x = k$x,y=(k$y-km$y)*scaler/max(km$y)+mean(ranges),col=cols['pca'],lwd=2)
     k<-density(mod.avg$predicted,from=ranges[1],to=ranges[2]);lines(x = k$x,y=(k$y-km$y)*scaler/max(km$y)+mean(ranges),col=cols['avg'],lwd=2)
     legend("topright",inset = c(0.01,0.01),
            legend = names(cols) ,
            col = cols,pch = 16)
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
     colnames(data)[1] <- 'DEPT'
     print(colnames(data))
     #browser()
     #data = read_table2(input$lasfile1$datapath, 
     #             col_types = cols(X11 = col_skip()))
     #data = data[complete.cases(data.frame(apply(data,c(1,2),FUN = function(x) {if(x==-999.25) NA else x}))),]
     #mbic = Mclust(q,G = 9);plot(mbic)
     #mcres = cbind(q,predict(mbic)$classification)
     #browser()
     data = data.frame(apply(data,c(1,2),FUN = function(x) {if(x==-999.25) return(NA) else return(x)}))
     # data$dBH = data$CALI
     # data$dBH[data$DEPT<3427] = 13+3/8
     # data$dBH[data$DEPT>=3427 & data$DEPT < 10877] = 9+5/8
     # data$dBH[data$DEPT>=10877 & data$DEPT < 13567] = 7
     # data$dBH[data$DEPT>=13567 & data$DEPT < 14663] = 5
     # data$dCALI = data$CALI - data$dBH
     # data$dCALI[data$dCALI<0]=NA
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

