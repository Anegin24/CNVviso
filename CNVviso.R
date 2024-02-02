library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)
library(DT)
library(ggtext)
library(tidyverse)
library(plotly)
library(smoother)
# UI
ui <- fluidPage(
  titlePanel("CNVviso: A tool for visualize CNV data from low-coverage WGS"),
  sidebarLayout(
    sidebarPanel(width=2,
      fileInput("file_dfcnr", "Choose LLR file(.cnr)", accept = c(".cnr")),
      fileInput("file_dfcns", "Choose Seqmentation file(.cns)", accept = c(".cns")),
      fileInput("file_dfcall", "Choose CNVcall file(.call.cns)", accept = c(".call.cns")),
      checkboxInput("displayTable", "Display Data", value = FALSE),
      downloadButton("downloadData", "Download Data"),
      sliderInput("xSlider", "Chromosome size", min = 0, max = 66, value = c(0, 66))
    ),
    mainPanel(width=10,
      plotOutput("myPlot", width = "100%", height = "400px"),
      DTOutput("dataTable",width= "50%")
    )
  ),
  fluidRow(
    column(2,plotOutput("myPlot1",width="100%",height="200px")),
    column(2,plotOutput("myPlot2",width="100%",height="200px")),
    column(2,plotOutput("myPlot3",width="100%",height="200px")),
    column(2,plotOutput("myPlot4",width="100%",height="200px")),
    column(2,plotOutput("myPlot5",width="100%",height="200px")),
    column(2,plotOutput("myPlot6",width="100%",height="200px"))),
  fluidRow(
    column(2,plotOutput("myPlot7",width="100%",height="200px")),
    column(2,plotOutput("myPlot8",width="100%",height="200px")),
    column(2,plotOutput("myPlot9",width="100%",height="200px")),
    column(2,plotOutput("myPlot10",width="100%",height="200px")),
    column(2,plotOutput("myPlot11",width="100%",height="200px")),
    column(2,plotOutput("myPlot12",width="100%",height="200px"))),
  fluidRow(
    column(2,plotOutput("myPlot13",width="100%",height="200px")),
    column(2,plotOutput("myPlot14",width="100%",height="200px")),
    column(2,plotOutput("myPlot15",width="100%",height="200px")),
    column(2,plotOutput("myPlot16",width="100%",height="200px")),
    column(2,plotOutput("myPlot17",width="100%",height="200px")),
    column(2,plotOutput("myPlot18",width="100%",height="200px"))),
  fluidRow(
    column(2,plotOutput("myPlot19",width="100%",height="200px")),
    column(2,plotOutput("myPlot20",width="100%",height="200px")),
    column(2,plotOutput("myPlot21",width="100%",height="200px")),
    column(2,plotOutput("myPlot22",width="100%",height="200px")),
    column(2,plotOutput("myPlot23",width="100%",height="200px")),
    column(2,plotOutput("myPlot24",width="100%",height="200px")))
)

server <- function(input, output) {
  # Reactive function to read uploaded dfcnr file
  dfcnr_data <- reactive({
    req(input$file_dfcnr)
    read.table(input$file_dfcnr$datapath, header = TRUE)
  })
  
  # Reactive function to read uploaded dfcns file
  dfcns_data <- reactive({
    req(input$file_dfcns)
    read.table(input$file_dfcns$datapath, header = TRUE)
  })
  
  dfcall_data <- reactive({
    req(input$file_dfcall)
    read.table(input$file_dfcall$datapath, header = TRUE)
  })
  
  observe({
    req(dfcnr_data(), dfcns_data(),dfcall_data())
    
    # Check if data is being loaded correctly
    print(dfcnr_data())
    print(dfcns_data())
    print(dfcall_data())
    
    dfcnr <- dfcnr_data()
    dfcns <- dfcns_data()
    dfcall <- dfcall_data()
    # Data processing
    dfcnr <- dfcnr %>%
      mutate(lrr = 2 * 2 ^ log2)
    for (i in 1:length(dfcnr$chromosome)) {
      dfcnr$order[i] <- i
    }
    
    dfcns <- dfcns %>%
      mutate(cncall = 2 * 2 ^ log2) %>%
      select(chromosome, start, end, cncall) %>%
      filter(chromosome != "chrY")
    dfcall <- dfcall %>% filter(chromosome == "chrY") %>% 
      select(chromosome,start,end,cn)
      names(dfcall)[names(dfcall) == "cn"] <- "cncall"
    dfcns<-merge(dfcall,dfcns,all=TRUE)
    df <- merge(dfcnr, dfcns, by = c("chromosome", "end"), all = TRUE) %>%
      select(-start.y)
      names(df)[names(df) == "start.x"] <- "start"
    df <- merge(df, dfcns, by.x = c("chromosome", "start"),
                by.y = c("chromosome", "start"), all = TRUE) %>%
      select(-end.y)
      names(df)[names(df) == "end.x"] <- "end"
    df[is.na(df)] <- 0
    df$cncall <- df$cncall.x + df$cncall.y
    df <- df %>% arrange(order)
    df <- df %>% select(chromosome, start, end, lrr, cncall)
    df <- df %>% mutate(lrr=ifelse(chromosome=="chrY",dfcall$cncall,lrr))
    
    CNvector <- dfcns$cncall
    CNcall <- as.numeric(df$cncall)
    
    start_index <- vector()
    end_index <- vector()
    
    for (i in 1:length(CNvector)) {
      start_index[i] <- which(CNcall == CNvector[i])[1]
    }
    for (i in 1:length(CNvector)) {
      end_index[i] <- which(CNcall == CNvector[i])[2]
    }
    for (i in 1:length(start_index)) {
      df$cncall[(start_index[i] + 1):(end_index[i] - 1)] <- df$cncall[start_index[i]]
    }
    hg38size<-data.frame(
      chromosome=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                   "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                   "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                   "chrX","chrY"),
      length=c(248956422,242193529,198295559,190214555,181538259,170805979,
               159345973,145138636,138394717,135086622,133797422,133275309,
               114364328,107043718,101991189,90338345,83257441,80373285,
               58617616,64444167,46709983,50818468,156040895,57227415)
    )
    striphg38 <- hg38size %>%
      rename_all(tolower) %>%
      mutate(pos = cumsum(length / min(length))) %>%
      mutate(ymin = 0, ymax = 4,
             xmin = c(0, pos[1:23]), xmax = pos,
             fill = rep(c("xam", "trang"), length.out = nrow(hg38size))) %>%
      pivot_longer(c(ymin, ymax), values_to = "y", names_to = "yy") %>% select(-yy)
    sizepros <- hg38size %>%
      rename_all(tolower) %>%
      mutate(propos = length / min(length)) %>% select(-length) %>%
      mutate(pos = cumsum(propos)) %>%
      mutate(c = c(0, pos[1:23])) %>%
      mutate(poschr = propos / 2 + c)
    
    data <- df %>%
      inner_join(., sizepros, by = "chromosome") %>%
      group_by(chromosome) %>%
      mutate(position = cumsum(propos / length(chromosome))) %>%
      mutate(point_x = c + position) %>%
      select(chromosome, lrr, cncall, point_x)
    
    
    data <- data %>%
      mutate(lrr = ifelse(lrr > cncall + 0.1, NA, lrr)) %>% na.omit()
    data <- data %>%
      mutate(lrr = ifelse(lrr < cncall - 0.1, NA, lrr)) %>% na.omit()
    cncallcount<-table(data$cncall)
    singletons<-names(cncallcount[cncallcount==1])
    data<-data[!data$cncall %in% singletons, ]
    data <- data %>%
      group_by(chromosome,cncall) %>%
      mutate(lrr = smth.gaussian(lrr)) %>% na.omit()
    
    table<-data%>%select(chromosome,lrr)
    names(table)[names(table) == "lrr"] <- "cn"  
    
    # Print the processed data
    print(data)
    
    # Render the plot
    output$myPlot <- renderPlot({
      data %>%
        ggplot(aes(x = point_x, y = lrr)) +
        geom_ribbon(data = striphg38,
                    aes(y = y, xmin = xmin, xmax = xmax, group = chromosome, fill = fill),
                    inherit.aes = FALSE,
                    show.legend = NULL) +
        geom_vline(size = 0.3, xintercept = c(unique(sizepros$pos)), linetype = 2, color = "grey90") +
        geom_line(linewidth = 0.3, color = "green") +
        geom_point(size = 0.9, shape = 21, color = "black", fill = "chartreuse3") +
        geom_hline(size = 0.3, yintercept = c(1, 3)) +
        geom_hline(size = 0.5, yintercept = 2, color = "green", linetype = "solid") +
        geom_vline(size = 0.3, xintercept = 61.550044) +
        scale_fill_manual(name = NULL,
                          breaks = c("xam", "trang"),
                          values = c("#F3F3F3", "#F3FAFE")) +
        scale_y_continuous(limits = c(0, 4),
                           breaks = seq(0, 4, by = 0.4),
                           expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0),
                           breaks = c(sizepros$poschr),
                           labels = c(sizepros$chromosome),
                           guide = guide_axis(angle = 45)) +
        labs(x = "Chromosomes", y = "Copy Number") +
        theme(
          axis.text.x = element_markdown(size=12),
          axis.text.y = element_markdown(size=12),
          plot.background = element_rect(fill = "#F3FAFE"),
          panel.background = element_blank()
        ) +
        theme_bw()+
        coord_cartesian(xlim = c(input$xSlider[1], input$xSlider[2]))
    })
    output$myPlot1 <- renderPlot({
      data%>%filter(chromosome == "chr1")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 1",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot2 <- renderPlot({
      data%>%filter(chromosome == "chr2")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 2",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot3 <- renderPlot({
      data%>%filter(chromosome == "chr3")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 3",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot4 <- renderPlot({
      data%>%filter(chromosome == "chr4")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 4",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot5 <- renderPlot({
      data%>%filter(chromosome == "chr5")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 5",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot6 <- renderPlot({
      data%>%filter(chromosome == "chr6")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 6",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot7 <- renderPlot({
      data%>%filter(chromosome == "chr7")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 7",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot8 <- renderPlot({
      data%>%filter(chromosome == "chr8")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 8",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot9 <- renderPlot({
      data%>%filter(chromosome == "chr9")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 9",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot10 <- renderPlot({
      data%>%filter(chromosome == "chr10")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 10",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot11 <- renderPlot({
      data%>%filter(chromosome == "chr11")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 11",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot12 <- renderPlot({
      data%>%filter(chromosome == "chr12")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 12",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot13 <- renderPlot({
      data%>%filter(chromosome == "chr13")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 13",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot14 <- renderPlot({
      data%>%filter(chromosome == "chr14")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 14",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot15 <- renderPlot({
      data%>%filter(chromosome == "chr15")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 15",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot16 <- renderPlot({
      data%>%filter(chromosome == "chr16")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 16",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot17 <- renderPlot({
      data%>%filter(chromosome == "chr17")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 17",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot18 <- renderPlot({
      data%>%filter(chromosome == "chr18")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 18",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot19 <- renderPlot({
      data%>%filter(chromosome == "chr19")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 19",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot20 <- renderPlot({
      data%>%filter(chromosome == "chr20")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 20",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot21 <- renderPlot({
      data%>%filter(chromosome == "chr21")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 21",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot22 <- renderPlot({
      data%>%filter(chromosome == "chr22")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes 22",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot23 <- renderPlot({
      data%>%filter(chromosome == "chrX")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes X",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$myPlot24 <- renderPlot({
      data%>%filter(chromosome == "chrY")%>%
        ggplot(aes(x=point_x,y=lrr))+
        geom_line(linewidth=0.3,color="green")+
        geom_point(size=0.9,shape=21,color="black",fill="chartreuse3")+
        geom_hline(size=0.3,yintercept = c(1,3))+
        geom_hline(size=0.5,yintercept=2,color="green",linetype="solid")+
        scale_fill_manual(name=NULL,
                          breaks=c("xam","trang"),
                          values=c("#F3F3F3","#F3FAFE"))+
        scale_y_continuous(limits=c(0,4),
                           breaks=seq(0,4,by=0.4),
                           expand=c(0,0))+
        scale_x_continuous(expand=c(0,0),
                           breaks=c(sizepros$chromopos),
                           labels=c(sizepros$chromosome),
                           guide = guide_axis(angle=45))+
        labs(x="Chromosomes Y",y="Copy Number" )+
        theme(
          plot.background = element_rect(fill="#F3FAFE"),
          panel.background = element_blank())+
        theme_bw()
    })
    output$dataTable <- renderDT({
      table
    })
    
    observe({
      if (input$displayTable %% 2 == 1) {
        output$dataTable <- renderDT({
          table
        })
      } else {
        output$dataTable <- renderDT({
          NULL
        })
      }
    })
  }) # Close the observe function here
  
  # Download data as CSV
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(finalresult, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
