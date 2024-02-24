library(readxl)
library(dplyr)
library(ggplot2)
library(DT)
library(ggtext)
library(tidyverse)
library(plotly)
library(smoother)
library(shiny)
library(shinyFiles)
library(fs)

ui <- fluidPage(
  titlePanel("CNVviso: A tool for visualizing CNV data from low-coverage WGS"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      tags$h4("Sample Zip File"),
      fileInput("zipfile", "Upload Zip File"),
      verbatimTextOutput("directorypath"),
      tags$hr(),
      checkboxInput("displayTable", "Display Data", value = FALSE),
      sliderInput("xSlider", "Chromosome size", min = 0, max = 66, value = c(0, 66))
    ),
    mainPanel(
      width = 10,
      plotlyOutput("myPlot"),
      column(6, DTOutput("table", width = "50%"))
    )
  )
)

server <- function(input, output, session) {
  # Extract files from the uploaded ZIP file
  observeEvent(input$zipfile, {
    req(input$zipfile)
    # Extract the name of the uploaded ZIP file
    zip_file_name <- gsub(".zip", "", input$zipfile$name)
    # Ensure that the extraction directory exists based on the name of the ZIP file
    extraction_directory <- paste0("unzipped_files_", zip_file_name)
    if (!dir.exists(extraction_directory)) {
      dir.create(extraction_directory)
    }
    unzip(input$zipfile$datapath, exdir = extraction_directory)
    output$directorypath <- renderPrint({
      "Process complete"
    })
  })
  
  # Define reactive dataframes for CNV files
  extracted_files <- reactive({
    req(input$zipfile)
    # Extract the name of the uploaded ZIP file
    zip_file_name <- gsub(".zip", "", input$zipfile$name)
    list.files(paste0("unzipped_files_", zip_file_name), full.names = TRUE)
  })
  
  # Define reactive expressions to read CNR, CNS, and CALL files
  dfcnr <- reactive({
    req(extracted_files())
    file_path_cnr <- extracted_files() %>% keep(., grepl(".dedup.cnr", .)) %>% first()
    read.table(file_path_cnr, header = TRUE)
  })
  
  dfcns <- reactive({
    req(extracted_files())
    file_path_cns <- extracted_files() %>% keep(., grepl(".dedup.cns", .)) %>% first()
    read.table(file_path_cns, header = TRUE)
  })
  
  dfcall <- reactive({
    req(extracted_files())
    file_path_call <- extracted_files() %>% keep(., grepl(".dedup.call.cns", .)) %>% first()
    read.table(file_path_call, header = TRUE)
  })
  # Display the name of the uploaded ZIP file
  output$directorypath <- renderPrint({
    if (is.null(input$zipfile)) {
      cat("Upload a zip file")
    } else {
      input$zipfile$name
    }
  })
  
  output$myPlot <- renderPlotly({
    req(dfcnr(),dfcns(),dfcall())
    dfcnr<-dfcnr()
    dfcns<-dfcns()
    dfcall<-dfcall()
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
      select(chromosome, start, end, lrr, cncall, point_x)
    
    
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
    
    p<-data %>%
      ggplot(aes(x = point_x, y = lrr))+
      geom_vline(size = 0.3, xintercept = c(unique(sizepros$pos)), linetype = 2, color = "grey90") +
      geom_hline(size=0.3,yintercept=4,linetype="solid")+
      geom_vline(size=0.3,xintercept=66,linetype="solid")+
      geom_line(linewidth = 0.3, color = "green",show.legend = FALSE) +
      geom_point(aes(text=paste("Chromosome:",chromosome,"<br>Start:",start,"<br>End:",end,"<br>CN:",lrr)),size = 0.9, shape = 21, color = "black", fill = "chartreuse3") +
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
      theme_classic()+
      coord_cartesian(xlim = c(input$xSlider[1], input$xSlider[2]))
    p<-ggplotly(p,tooltip = "text") %>% layout(xaxis = list(tickangle = -45))
    p
  })
  output$table <- renderDT({
    if(input$displayTable){
      req(dfcnr(),dfcns(),dfcall())
      dfcnr<-dfcnr()
      dfcns<-dfcns()
      dfcall<-dfcall()
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
        mutate(lrr = smth.gaussian(lrr)) %>% na.omit() %>% ungroup()
      
      table<-data%>%select(chromosome,lrr)
      names(table)[names(table) == "lrr"] <- "cn"
      table
    }
  })
}

shinyApp(ui, server)

