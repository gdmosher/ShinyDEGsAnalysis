# module-DEGsAnalysis.R

# Shiny module definition for module.systemPipeR contains:
#   ui, previously named ui-tab-systemPipeR
#   helper functions, prev. fun-systemPipeR
#   server, previously   server-systemPipeR
# source("this file") in global.R
# connect this ui function in the main app as you would a panel
# don't forget to call the server function in the main app server.R
# also need ->   ns <- session$ns  # ns for renderUI() calls in server


#   ui, previously named ui-tab-systemPipeR
# ==============================================================================================================================
# ==============================================================================================================================
module.DEGsAnalysisInput6 <- function(id, p1='') { # p1 is optional userdef param
  # PARAMFILENAME.COMMENTS <- paste0(PARAMFILENAME, ".comments")
  # mandatory to postpend the module name with "Input"
  # create a namespace function using the provided id
  ns <- NS(id)
  # now wrap the ui in tagList()
  # and wrap ALL input and output objects in ns()
  # cw = 2
  tagList(
    
    # tabPanel(PARAMFILENAME, # namespace .NAMESPACENAME matches tabPanel title (in case and punctuation, plus dot in front)
    #          PARAMFILENAME,
    fluidPage(
      # Application title
      titlePanel(paste0(p1,"Analysis of DEGs")),
      
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          fluidRow(
            # "Text on mainPanel\n\n"
          )
          , fluidRow(
            # , textInput(ns('command.444'), label='Enter Command')
            # actionButton(ns('execute_actionButton.6660'), label='Upload counts table')
            fileInput(ns('countsTable'),'Upload counts table', placeholder = 'countDFeByg.xls')
            
            # , actionButton(ns('execute_actionButton.6661'), label='Run edgeR and biomaRt (usu. < 2 min.)')
            # , actionButton(ns('execute_actionButton.6662'), label='Load edgeR output, make filterDEGs() plot')
            , sliderInput(ns("render_sliderInput.6661"), "Fold Slider", -2, 2, 2)
            , sliderInput(ns("render_sliderInput.6662"), "FDR Slider", 1, 20, 10, step = 1)
            # , textInput(ns('envir.444'), label='Enter name=value pairs to set environment variables (quoting is weird)')
            # , verbatimTextOutput(ns("text.666"))
            # , tags$iframe(style="height:400px; width:60%; scrolling=yes",
            #               src="DEGcounts.pdf") # must be in www folder?? if so, how to copy it there now?
            # , box(plotOutput(ns("plot6661"), height = 250, width = 750),width = 12)
          )
          # , fluidRow (workDir)
          , fluidRow ("Bar plot")
          , fluidRow (plotOutput(ns("plot6661"), height = 250, width = 750))
          , fluidRow ("Venn diagram")
          , fluidRow (imageOutput(ns("plot6662"), height = 500, width = 12))
          # , fluidRow (box(imageOutput(ns("plot6663"), height = 500, width = 12))
          # )
        )
      )
    ) # end fluid page
    
  ) # tagList
} # module.DEGsAnalysisInput6

# ==============================================================================================================================
# ==============================================================================================================================
module.DEGsAnalysisInput8 <- function(id, p1='') { # p1 is optional userdef param
  # PARAMFILENAME.COMMENTS <- paste0(PARAMFILENAME, ".comments")
  # mandatory to postpend the module name with "Input"
  # create a namespace function using the provided id
  ns <- NS(id)
  # now wrap the ui in tagList()
  # and wrap ALL input and output objects in ns()
  # cw = 2
  tagList(

    # tabPanel(PARAMFILENAME, # namespace .NAMESPACENAME matches tabPanel title (in case and punctuation, plus dot in front)
    #          PARAMFILENAME,
    fluidPage(
      # Application title
      titlePanel(paste0(p1,"Clustering and heat maps")),

      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          fluidRow(
            # "Text on mainPanel\n\n"
          )
          , fluidRow(
            # , textInput(ns('command.444'), label='Enter Command')
            # actionButton(ns('execute_actionButton.6660_2'), label='Calculate dds <- DESeqDataSetFromMatrix()')
            sliderInput(ns("render_sliderInput.6661_2"), "Fold Slider", -2, 2, 2)
            , sliderInput(ns("render_sliderInput.6662_2"), "FDR Slider", 1, 20, 10, step = 1)
            # , textInput(ns('envir.444'), label='Enter name=value pairs to set environment variables (quoting is weird)')
            # , verbatimTextOutput(ns("text.666"))
            # , tags$iframe(style="height:400px; width:60%; scrolling=yes",
            #               src="DEGcounts.pdf") # must be in www folder?? if so, how to copy it there now?
            # , box(plotOutput(ns("plot6661"), height = 250, width = 750))
          )
          # , fluidRow (box(imageOutput(ns("plot6662"), height = 500, width = 12))
          # )
          # , fluidRow (box(imageOutput(ns("plot6664"), height = 500, width = 12)) # for phylo - doesn't belong here
          # )
          , fluidRow ("Heat map")
          , fluidRow (imageOutput(ns("plot6663"), height = 500, width = 12)
          )
        )
      )
    ) # end fluid page

  ) # tagList
} # module.DEGsAnalysisInput8


#   helper functions, prev. fun-systemPipeR

#   server, previously   server-systemPipeR
# Notice: wrapper for server functions

# ==============================================================================================================================
# ==============================================================================================================================
# Module server function
module.DEGsAnalysis <- function(input, output, session, p1='') {
  print("***** module-DEGsAnalysis.R - begin module.DEGsAnalysis")
  ####### code here does not run before the UI
  #######      code and variables needed to init UI belong in global.R
  library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
  
  # ##########################################################################
  # Tab 6. Analysis of DEGs
  # ##########################################################################
  runChunk_edgeR_and_biomaRt <- function() {
    start_time <- Sys.time()
    print(start_time)
    
    withProgress(message = 'Processing edgeR_and_biomaRt < 2 min. remaining', value = 0, {
      
      # print(getwd())
      setwd(workDir)
      # try(print(dir())) # works
      
      ## WE ARE RUNNING SOME SYSTEMPIPER CODE LOCAL HERE
      library(edgeR)
      print(countsTablePath)
      countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
      targets <- read.delim("targets.txt", comment="#")
      cmp <- readComp(file="targets.txt", format="matrix", delim="-")
      edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=FALSE, mdsplot="")
      print('run_edgeR() is done.')
      incProgress(.30)
      
      
      # library("biomaRt")
      ## biomaRt vs dplyr issue https://github.com/tidyverse/dplyr/issues/756
      m <- biomaRt::useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
      desc <- biomaRt::getBM(attributes=c("tair_locus", "description"), mart=m)
      desc <- desc[!duplicated(desc[,1]),]
      descv <- as.character(desc[,2]); names(descv) <- as.character(desc[,1])
      edgeDF <- data.frame(edgeDF, Desc=descv[rownames(edgeDF)], check.names=FALSE)
      write.table(edgeDF, "results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
      # write.table(edgeDF, "results/edgeRglm_allcompGORDON.xls", quote=FALSE, sep="\t", col.names = NA)
      print('biomaRt call is done.')
      # unloadNamespace("biomaRt")
      
      incProgress(.90)
      Sys.sleep(0.9)
      
      setwd('../')
    })
    end_time <- Sys.time()
    print(end_time - start_time)
    
  }

  runChunk_DESeqDataSetFromMatrix <- function() {
    start_time <- Sys.time()
    print(start_time)
    
    # print(getwd())
    setwd(workDir)
    # try(print(dir())) # works

    ## =================================================================================
    
    ## ----targetsSE, eval=TRUE------------------------------------------------
    library(systemPipeR)
    targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
    read.delim(targetspath, comment.char = "#")
    
    ## ----targetsPE, eval=TRUE------------------------------------------------
    targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
    read.delim(targetspath, comment.char = "#")[1:2,1:6]
    
    ## ----comment_lines, eval=TRUE--------------------------------------------
    readLines(targetspath)[1:4]
    
    ## ----targetscomp, eval=TRUE----------------------------------------------
    readComp(file=targetspath, format="vector", delim="-")
    
    ## ----param_structure, eval=TRUE------------------------------------------
    parampath <- system.file("extdata", "tophat.param", package="systemPipeR")
    print(parampath)
    read.delim(parampath, comment.char = "#")
    
    ## ----param_import, eval=TRUE---------------------------------------------
    args <- suppressWarnings(systemArgs(sysma=parampath, mytargets=targetspath))
    args
    
    ## ----sysarg_access, eval=TRUE--------------------------------------------
    names(args)
    modules(args)
    cores(args)
    outpaths(args)[1]
    sysargs(args)[1]
    
    ## ----sysarg_json, eval=TRUE----------------------------------------------
    systemArgs(sysma=parampath, mytargets=targetspath, type="json")
    ## =================================================================================
    withProgress(message = 'Processing DESeqDataSetFromMatrix < 2 min. remaining', value = 0, {
      ## WE ARE RUNNING SOME SYSTEMPIPER CODE LOCAL HERE
      # sample_tree code chunk copied from systemPipeRNAseq.Rmd
      library(DESeq2, quietly=TRUE); library(ape,  warn.conflicts=FALSE)
      countDF <- as.matrix(read.table("results/countDFeByg.xls"))
      colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
      dds <<- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
      save(dds, file = "results/dds.Rdata")
      print("Calculate dds <- DESeqDataSetFromMatrix() is done.")
      
      # d <- cor(assay(rlog(dds)), method="spearman")
      # hc <- hclust(dist(1-d))
      # pdf("results/sample_tree.pdf")
      # plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
      # dev.off()
      
      incProgress(.90)
      Sys.sleep(0.9)
      
      setwd('../')
    })
    end_time <- Sys.time()
    print(end_time - start_time)
    
  }
  
  countsTablePath <<- paste0(workDir, "/results/countDFeByg.xls") # force load this file on startup
  # could speed this up by loading edgeRglm_allcompTOY.xls and ddsTOY.Rdata with countDFeBygTOY.xls
  # runChunk_DESeqDataSetFromMatrix() # must run when loading new counts
  # runChunk_edgeR_and_biomaRt()      # must run when loading new counts

  countsTablePath <<- paste0(workDir, "/results/countDFeBygTOY.xls") # force load this file on startup
  load(paste0(workDir, "/results/ddsTOY.Rdata"), verbose = TRUE) # load dds from file
  edgeDF <<- read.delim(paste0(workDir, "/results/edgeRglm_allcompTOY.xls"), row.names=1, check.names=FALSE)
  
  
  observeEvent(input$countsTable, {
    print(input$countsTable$datapath) # new filename
    countsTablePath <<- input$countsTable$datapath   # 
    # updateTextAreaInput(session, "file_textAreaInput.99", value = includeText(input$filename_fileInput.99$datapath))
    countDF <- read.delim(countsTablePath, row.names=1, check.names=FALSE)
    write.table(countDF, paste0(workDir, "/results/countDFeByg.xls"), sep = "\t", quote = FALSE)
    countsTablePath <<- paste0(workDir, "/results/countDFeByg.xls")
    runChunk_DESeqDataSetFromMatrix() # must run when loading new counts
    runChunk_edgeR_and_biomaRt()      # must run when loading new counts
  })
  # click Submit button to 'Run edgeR and biomaRt (usu. < 2 min.)'
  observeEvent(input$execute_actionButton.6661, {
    runChunk_edgeR_and_biomaRt()
    
  }) # , priority = 1)
  
  # output$text.666 <- eventReactive(
  observeEvent(input$execute_actionButton.6662, {
    start_time <- Sys.time()
    print(start_time)
    
    print('Running observeEvent')
    
    p1 <- input$render_sliderInput.6661
    print(p1)
    # updateSliderInput(session, "render_sliderInput.6661", label='hello', value = p1-1 )
    p2 <- input$render_sliderInput.6662
    print(p2)
    
    print('...')
    # p1 <- 0
    # p2 <- .15
    my_filter=c(Fold=p1, FDR=p2)
    # my_filter=c(Fold=2, FDR=30)
    print(my_filter)
    
    make_DEG_list()
    # edgeDF <- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
    # DEG_list <- filterDEGs(degDF=edgeDF, filter=my_filter)
    pdf(paste0(workDir, "/results/DEGcounts.pdf"))
    make_DEG_list()
    # DEG_list <- filterDEGs(degDF=edgeDF, filter=my_filter)
    dev.off()
    write.table(DEG_list$Summary, paste0(workDir, "/results/DEGcounts.xls"), quote=FALSE, sep="\t", row.names=FALSE)
    # write.table(DEG_list$Summary, paste0(workDir, "/results/DEGcountsGORDON.xls"), quote=FALSE, sep="\t", row.names=FALSE)
    print('filterDEGs() call is done.')
    
    end_time <- Sys.time()
    print(end_time - start_time)
  })
  
  # globals for visualizations
  my_filter <- ''
  edgeDF <- ''
  DEG_list <- ''
  dds <- ''
  make_DEG_list <- function() {
    p1 <- input$render_sliderInput.6661
    print(p1)
    p2 <- input$render_sliderInput.6662
    print(p2)
    
    print('...')
    # p1 <- 0
    # p2 <- .15
    my_filter=c(Fold=p1, FDR=p2)
    # my_filter=c(Fold=2, FDR=30)

    print(getwd())
    print(getwd())
    print(getwd())
    print(getwd())
    # edgeDF <<- read.delim(paste0('', "results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
    edgeDF <<- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
    DEG_list <<- filterDEGs(degDF=edgeDF, filter=my_filter)
    print("Returning from make_DEG_list()")
  }
  
  output$plot6661 <- renderPlot({

    # new code copied from above
    # make_DEG_list()

    withProgress(message = 'Making bar plot < 1 min. remaining', value = 0, {
      # reactive inputs
      input$countsTable
      input$execute_actionButton.6661
      input$execute_actionButton.6662
      p1 <- input$render_sliderInput.6661
      print(p1)
      p2 <- input$render_sliderInput.6662
      print(p2)
      print('...')
      # p1 <- 0
      # p2 <- .15
      my_filter=c(Fold=p1, FDR=p2)
      # my_filter=c(Fold=2, FDR=30)
      print(my_filter)
      make_DEG_list() # this line seems like a dup, but it is the one that works?!?
      # edgeDF <- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
      # DEG_list <- filterDEGs(degDF=edgeDF, filter=my_filter) # this line seems like a dup, but it is the one that works?!?
      
      pdf(paste0(workDir, "/results/DEGcounts.pdf"))
      make_DEG_list()
      # DEG_list <- filterDEGs(degDF=edgeDF, filter=my_filter) # filterDEGs() makes the bar plot pdf for Rmd
      dev.off()
      write.table(DEG_list$Summary, paste0(workDir, "/results/DEGcounts.xls"), quote=FALSE, sep="\t", row.names=FALSE)
      # write.table(DEG_list$Summary, paste0(workDir, "/results/DEGcountsGORDON.xls"), quote=FALSE, sep="\t", row.names=FALSE)
      print('filterDEGs() call is done.')
    })
  })
  
  output$plot6662 <- renderImage({
    # make_DEG_list()
    withProgress(message = 'Making venn plot < 1 min. remaining', value = 0, {
      # reactive inputs
      input$countsTable
      input$execute_actionButton.6661
      input$execute_actionButton.6662
      p1 <- input$render_sliderInput.6661
      print(p1)
      p2 <- input$render_sliderInput.6662
      print(p2)
      
      print('...')
      # p1 <- 0
      # p2 <- .15
      my_filter=c(Fold=p1, FDR=p2)
      # my_filter=c(Fold=2, FDR=30)
      print(my_filter)
      make_DEG_list()
      # edgeDF <- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
      # DEG_list <- filterDEGs(degDF=edgeDF, filter=my_filter)
      
      vennsetup <- overLapper(DEG_list$Up[6:9], type="vennsets")
      vennsetdown <- overLapper(DEG_list$Down[6:9], type="vennsets")
      png(paste0(workDir, "/results/vennplot.png"))
      vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
      dev.off()
      
    })
    # vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
    filename <- paste0(workDir, "/results/vennplot.png")
    list(src = filename, alt = "vennplot.png")
    
  }, deleteFile = FALSE)

  # ##########################################################################
  # Tab 8. Clustering and heat maps
  # ##########################################################################
  # click Submit button to 'Calculate dds <- DESeqDataSetFromMatrix()'
  observeEvent(input$execute_actionButton.6660_2, {
    runChunk_DESeqDataSetFromMatrix()
  }) # , priority = 1)
  

  output$plot6663 <- renderImage({
    start_time <- Sys.time()
    print(start_time)
    print("Ready to start plot6663")
    withProgress(message = 'Making heat map < 1 min. remaining', value = 0, {
      
      # make_DEG_list()
      # reactive inputs
      input$countsTable
      input$execute_actionButton.6660_2
      p1 <- input$render_sliderInput.6661_2
      print(p1)
      p2 <- input$render_sliderInput.6662_2
      print(p2)
      
      print('...')
      # p1 <- 0
      # p2 <- .15
      my_filter=c(Fold=p1, FDR=p2)
      # my_filter=c(Fold=2, FDR=30)
      print(my_filter)
      # run filterDEGs() whenever filter changes
      # DONT CALL make_DEG_list() FROM HERE *** SLIDERS ARE DIFFERENT
      edgeDF <<- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
      DEG_list <<- filterDEGs(degDF=edgeDF, filter=my_filter)
      print(getwd())
      load(paste0(workDir, "/results/dds.Rdata"), verbose = TRUE) # load dds from file
      incProgress(.10)
      end_time <- Sys.time()
      print(end_time - start_time)
      # heatmap code chunk copied from systemPipeRNAseq.Rmd
      library(pheatmap)
      print("Ready to get geneids")
      geneids <- unique(as.character(unlist(DEG_list[[1]])))
      print("Ready to perform assay (usu. < 1 min.)")
      incProgress(.25)
      y <- assay(rlog(dds))[geneids, ]
      end_time <- Sys.time()
      print(end_time - start_time)
      incProgress(.90)
      Sys.sleep(0.9)
      
      print("Ready to make heatmap")
      print(getwd())
      # print(y)
      png(paste0(workDir, "/results/heatmap1.png"))
      pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
      dev.off()
      png(paste0(workDir, "/results/heatmap1.png"))
      pheatmap(y, scale="row", clustering_distance_rows="correlation", clustering_distance_cols="correlation")
      dev.off()
      incProgress(1)
      
    })
    end_time <- Sys.time()
    print(end_time - start_time)
    
    filename <- paste0(workDir, "/results/heatmap1.png")
    list(src = filename, alt = "heatmap1.png")
    
  }, deleteFile = FALSE)
  
  output$plot6664 <- renderImage({
    start_time <- Sys.time()
    print(start_time)
    print("Ready to start plot6664")
    
    
    # make_DEG_list()
    p1 <- input$render_sliderInput.6661_2
    print(p1)
    p2 <- input$render_sliderInput.6662_2
    print(p2)
    
    print('...')
    # p1 <- 0
    # p2 <- .15
    my_filter=c(Fold=p1, FDR=p2)
    # my_filter=c(Fold=2, FDR=30)
    print(my_filter)
    edgeDF <<- read.delim(paste0(workDir, "/results/edgeRglm_allcomp.xls"), row.names=1, check.names=FALSE)
    DEG_list <<- filterDEGs(degDF=edgeDF, filter=my_filter)
    # load(paste0(workDir, "/results/dds.Rdata"), verbose = TRUE) # load dds from file
    
    end_time <- Sys.time()
    print(end_time - start_time)
    # phylo code copied from
    print("Ready to perform assay (usu. < 1 min.)")
    d <- cor(assay(rlog(dds)), method="spearman")
    end_time <- Sys.time()
    print(end_time - start_time)
    print("Ready to perform hclust")
    hc <- hclust(dist(1-d))
    print("Ready to make phylo")
    png(paste0(workDir, "/results/sample_tree.png"))
    plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
    dev.off()
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    filename <- paste0(workDir, "/results/sample_tree.png")
    list(src = filename, alt = "sample_tree.png")
    
  }, deleteFile = FALSE)
}

# 