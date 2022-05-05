library(shiny)
library(mongolite)
library(gdata)
library(dendextend)
library(heatmaply)
library(seriation)
library(scales)

# userData <- reactiveValues()
#     userData$object <- NULL


#server <- function(input, output, session) {
  #output$square <- reactive(input$plot, {
  #  input$plot
  #})
  #outputOptions(output, 'square', suspendWhenHidden = FALSE)
observeEvent(input$about_link, {
  session$sendCustomMessage("openTitleInfoBox", runif(1))
})
  genes_collection <- c()
  cell_collection <- c()
  xumap3d_collection <- c()
  xumap2d_collection <- c()
  xtsne3d_collection <- c()
  xtsne2d_collection <- c()
  colls <- c()
  text <- c()
  xumap3d <- c()
  xumap2d <- c()
  xtsne3d <- c()
  xtsne2d <- c()
  gColi <- c()
  gColi_MTG <- c()
  cell_count <- c()
  #mongodb_url <- "mongodb://host.docker.internal:27017"
  mongodb_url <- "mongodb://mongo:27017"
  observeEvent(input$dataset,{
    if(input$dataset == "iPSC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "iPSC", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "iPSC", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "iPSC", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "iPSC", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "iPSC", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "iPSC", url = mongodb_url)
      cell_count <<- cell_count_iPSC
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "H. EC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "EC_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "EC_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "EC_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "EC_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "EC_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "EC_eis4", url = mongodb_url)
      cell_count <<- cell_count_EC
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Leng et al. H. EC_Exc"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "EC_Exc", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "EC_Exc", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "EC_Exc", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "EC_Exc", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "EC_Exc", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "EC_Exc", url = mongodb_url)
      cell_count <<- cell_count_EC_Exc
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "H. EC_Exc_0.7"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "EC_Exc", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell_0.7", db = "EC_Exc", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "EC_Exc", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "EC_Exc", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "EC_Exc", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "EC_Exc", url = mongodb_url)
      cell_count <<- cell_count_EC_Exc
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Leng et al. H. EC_Inh"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "EC_Inh", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "EC_Inh", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "EC_Inh", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "EC_Inh", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "EC_Inh", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "EC_Inh", url = mongodb_url)
      cell_count <<- cell_count_EC_Inh
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "H. SFG"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "SFG_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "SFG_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "SFG_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "SFG_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "SFG_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "SFG_eis4", url = mongodb_url)
      cell_count <<- cell_count_SFG
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Leng et al. H. SFG_Exc"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "SFG_Exc", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "SFG_Exc", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "SFG_Exc", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "SFG_Exc", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "SFG_Exc", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "SFG_Exc", url = mongodb_url)
      cell_count <<- cell_count_SFG_Exc
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Leng et al. H. SFG_Inh"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "SFG_Inh", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "SFG_Inh", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "SFG_Inh", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "SFG_Inh", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "SFG_Inh", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "SFG_Inh", url = mongodb_url)
      cell_count <<- cell_count_SFG_Inh
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Allen H. MTG"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_MTG_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_MTG_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_MTG_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_MTG_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_MTG_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_MTG_eis4", url = mongodb_url)
      cell_count <<- cell_count_MTG
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Allen H. VC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_V1_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_V1_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_V1_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_V1_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_V1_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_V1_eis4", url = mongodb_url)
      cell_count <<- cell_count_V1_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Allen H. ACC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_ACC_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_ACC_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_ACC_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_ACC_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_ACC_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_ACC_eis4", url = mongodb_url)
      cell_count <<- cell_count_ACC_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "SLC17A7")
    }
    if(input$dataset == "Allen M. VC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_VISp_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_VISp_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_VISp_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_VISp_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_VISp_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_VISp_eis4", url = mongodb_url)
      cell_count <<- cell_count_VISp_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "Apoe")
    }
    if(input$dataset == "Allen M. ALM"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_ALM_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_ALM_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_ALM_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_ALM_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_ALM_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_ALM_eis4", url = mongodb_url)
      cell_count <<- cell_count_ALM_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "Apoe")
    }
    if(input$dataset == "Allen M. ACA"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_ACA_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_ACA_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_ACA_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_ACA_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_ACA_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_ACA_eis4", url = mongodb_url)
      cell_count <<- cell_count_ACA_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "Apoe")
    }
    if(input$dataset == "Allen M. MOpC"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_MOpC_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_MOpC_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_MOpC_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_MOpC_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_MOpC_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_MOpC_eis4", url = mongodb_url)
      cell_count <<- cell_count_MOpC_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "Apoe")
    }
    if(input$dataset == "Allen M. MOpN"){
      genes_collection <<- mongolite::mongo(collection = "genes", db = "allen_MOpN_eis4", url = mongodb_url)
      cell_collection <<- mongolite::mongo(collection = "cell", db = "allen_MOpN_eis4", url = mongodb_url)
      xumap3d_collection <<- mongolite::mongo(collection = "xumap3d", db = "allen_MOpN_eis4", url = mongodb_url)
      xumap2d_collection <<- mongolite::mongo(collection = "xumap2d", db = "allen_MOpN_eis4", url = mongodb_url)
      xtsne3d_collection <<- mongolite::mongo(collection = "xtsne3d", db = "allen_MOpN_eis4", url = mongodb_url)
      xtsne2d_collection <<- mongolite::mongo(collection = "xtsne2d", db = "allen_MOpN_eis4", url = mongodb_url)
      cell_count <<- cell_count_MOpN_eis4
      updateSelectizeInput(session, 'gene', choices = genes_collection$find()[[1]], server = TRUE,
                           selected = "Apoe")
    }
    xumap3d <<- xumap3d_collection$find()
    xumap2d <<- xumap2d_collection$find()
    xtsne3d <<- xtsne3d_collection$find()
    xtsne2d <<- xtsne2d_collection$find()
    i <- 1
    colls <<- c()
    while(i < cell_count){
      colls <<- c(colls, i)
      i <- i+63
    }
    #hov = colnames(cell_collection$find()[0,])
    hov = colnames(cell_collection$find())
    text <<- c()
    for (i in 1:length(hov)){
      text <<- paste(text,hov[i],':',cell_collection$find('{}', fields = paste0('{"_id":0,"',hov[i],'":1}'))[[1]],'<br>', sep = "")
    }
    gColi <<- reorder.factor(as.factor(cell_collection$find('{}', fields = paste0('{"_id":0,"',"cluster",'":1}'))[[1]]))
    if(input$dataset == "Allen H. MTG"){
      gColi_MTG <- cell_collection$find('{}', fields = paste0('{"_id":0,"',"cluster",'":1}'))[[1]]
    }  
    for (item in c("differentialsCluster1Sel", "differentialsCluster2Sel")) {
      updateSelectizeInput(session,
                           item,
                           choices  = gColi,
                           selected = "",
                           server   = FALSE)
    }
  })
  

  m_collection <- c()
  g.ndx_new <- c()
  g <- c()
  sizei <- c()
  Expressioni <- c()
  observeEvent(input$gene,{
    g <<-input$gene
    g.ndx <- grep(paste("\\b",g,"\\b",sep=""), x = genes_collection$find()[[1]], ignore.case = T)
    #print(length(g.ndx))
    if(length(g.ndx) >= 1) { 
      label <- 0
      if(length(g.ndx) <= 10){
      for(ndx in g.ndx){
        print(toupper(genes_collection$find()[[ndx,1]]))
        print(toupper(g))
        if(toupper(genes_collection$find()[[ndx,1]])==toupper(g)){
          label <- 1
          g.ndx <- ndx
          g <<- toupper(genes_collection$find()[[ndx,1]])
        }
      }
      }
      #print(g)
      if(label == 0){
        g.ndx <- g.ndx[1]
        g <<- toupper(genes_collection$find()[[g.ndx,1]])
      }
    }
    j <- 1
    m_collection <<- c()
    g.ndx_new <<- 0
    for(i in colls){
      if(j+1 <= length(colls)){
        if((g.ndx >= i)&(g.ndx <= colls[j+1])){
          if(input$dataset == "iPSC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "iPSC", url = mongodb_url)
          }
          if(input$dataset == "H. EC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Exc"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "H. EC_Exc_0.7"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Inh"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Inh", url = mongodb_url)
          }
          if(input$dataset == "H. SFG"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Exc"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Inh"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_Inh", url = mongodb_url)
          }
          if(input$dataset == "Allen H. MTG"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MTG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. VC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_V1_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. ACC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ACC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. VC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_VISp_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ALM"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ALM_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ACA"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ACA_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MOpC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpN"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MOpN_eis4", url = mongodb_url)
          }
          g.ndx_new <<- (g.ndx-i)+1
        }
      }else{
        if(g.ndx >= i){
          if(input$dataset == "iPSC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "iPSC", url = mongodb_url)
          }
          if(input$dataset == "H. EC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Exc"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "H. EC_Exc_0.7"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Inh"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "EC_Inh", url = mongodb_url)
          }
          if(input$dataset == "H. SFG"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Exc"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Inh"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "SFG_Inh", url = mongodb_url)
          }
          if(input$dataset == "Allen H. MTG"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MTG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. VC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_V1_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. ACC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ACC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. VC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_VISp_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ALM"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ALM_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ACA"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_ACA_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpC"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MOpC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpN"){
            m_collection <<- mongolite::mongo(collection = paste0("m",i), db = "allen_MOpN_eis4", url = mongodb_url)
          }
          g.ndx_new <<- (g.ndx-i)+1
        }
      }
      j <- j+1
    }
    if(length(g.ndx)==0) {g <<- "APOE"}
    sizei <<- scale(m_collection$find('{}', fields = paste0('{"_id":0,"V',g.ndx_new,'":1}'))[[1]], center = T, scale = T)
    Expressioni <<- m_collection$find('{}', fields = paste0('{"_id":0,"V',g.ndx_new,'":1}'))[[1]]
  })

  # render.h <- eventReactive(input$heatmapPlotBtn, {
  #   gradient_col <- ggplot2::scale_fill_gradient(
  #     low = "white", high = "red"
  #   )
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. VC")){
  #     heatmaply(mtg_v1_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. ACC")){
  #     heatmaply(mtg_acc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     #rownames(mtg_ec_exc_odds.ratio) <- gColi_MTG
  #     heatmaply(mtg_ec_exc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     heatmaply(mtg_ec_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen H. ACC")){
  #     heatmaply(v1_acc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     heatmaply(v1_ec_exc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     heatmaply(v1_ec_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     heatmaply(acc_ec_exc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     heatmaply(acc_ec_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     heatmaply(ec_exc_ec_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #     heatmaply(ec_exc_sfg_exc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     heatmaply(ec_exc_sfg_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #     heatmaply(ec_inh_sfg_exc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     heatmaply(ec_inh_sfg_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. SFG_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     heatmaply(sfg_exc_sfg_inh_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(visp_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(visp_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(visp_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(visp_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(alm_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(alm_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(alm_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(aca_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(aca_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen M. MOpC")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(mopc_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
    
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. VC")){
  #     heatmaply(mtg1_visp_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(mtg1_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(mtg1_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(mtg1_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(mtg1_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
    
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. VC")){
  #     heatmaply(v11_visp_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(v11_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(v11_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(v11_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(v11_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
    
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. VC")){
  #     heatmaply(acc1_visp_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(acc1_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(acc1_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(acc1_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(acc1_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }

  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. VC")){
  #     heatmaply(ec_exc1_visp_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=1000)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(ec_exc1_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(ec_exc1_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(ec_exc1_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(ec_exc1_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
    
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. VC")){
  #     heatmaply(ec_inh1_visp_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ALM")){
  #     heatmaply(ec_inh1_alm_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ACA")){
  #     heatmaply(ec_inh1_aca_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpC")){
  #     heatmaply(ec_inh1_mopc_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  #   else if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpN")){
  #     heatmaply(ec_inh1_mopn_odds.ratio, scale_fill_gradient_fun = gradient_col, dendrogram  = "both", show_dendrogram = FALSE) %>% layout(height=600)
  #   }
  # })
  # output$heatmap <- renderPlotly({
  #   render.h()
  # })
  
  render.d <- eventReactive(input$go, {
    #req(input$gene)
    colrs <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951', '#F06AF0', '#E8AB4F', '#E69138', '#98EA74', '#EA53A3')
    #hov = colnames(cell)[!colnames(cell)%in%"cluster"]
    #hov = c("donor_id", "gender", "cluster", "sample_type")
    
    text=paste('<span style="font-size: 120%;>',text,'</span>')
    if((input$plot=="UMAP") & (input$dim=="3D") & (input$lay=="Clusters")) {
      df <- data.frame(xumap3d, gCol = gColi, size = sizei)
      # text = paste('<span style="font-size: 120%;">Age:', cell$age, '<br>Gender:', cell$gender, '<br>Genotype:', cell$genotype, '</span>')
      t <- list(
        family = "sans serif",
        size = 14,
        color = toRGB("grey50"))
      plot_ly(data = df, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, source = "clust",
              sizes=c(5,50), size=~size, color = ~gCol,
              marker = list(symbol = 'circle', sizemode = 'diameter'),
              colors = colrs, type = "scatter3d", 
              height = as.numeric(input$dimension[2]), 
              mode = "markers", showlegend=T, hoverinfo = 'text',
              text = text
              ) %>% 
        layout(title=g, 
               paper_bgcolor = 'rgb(243, 243, 243)', 
               plot_bgcolor = 'rgb(243, 243, 243)', 
               scene = list(xaxis = list(showspikes=FALSE), 
                            yaxis = list(showspikes=FALSE),
                            zaxis = list(showspikes=FALSE),
                            camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25))
               )
        )
    } else if((input$plot=="UMAP") & (input$dim=="2D") & (input$lay=="Clusters")) {
      df <- data.frame(xumap2d, gCol = gColi, size = sizei)
      
      # text = paste('<span style="font-size: 120%;">Age:', cell$age, '<br>Gender:', cell$gender, '<br>Genotype:', cell$genotype, '</span>')
      plot_ly(data = df, x = ~UMAP_1, y = ~UMAP_2, source = "clust",
              color = ~gCol,
              colors = colrs,
              sizes=c(5,50), size=~size,  marker = list(symbol = 'circle', sizemode = 'diameter'),
              type = "scattergl", mode="markers", showlegend=T, 
              height = as.numeric(input$dimension[2]),
              hoverinfo = 'text',
              text = text
              ) %>% 
        layout(title=g, paper_bgcolor = 'rgb(243, 243, 243)', plot_bgcolor = 'rgb(243, 243, 243)')
  }else if((input$plot=="UMAP") & (input$dim=="3D") & (input$lay=="Gene Expression")) {
    df <- data.frame(xumap3d, Expression = Expressioni)
    plot_ly(data = df, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, source = "exp",
            sizes=c(60,60), size=60, color = ~Expression,
            colors = c("firebrick4", "yellow"), type = "scatter3d", 
            height = as.numeric(input$dimension[2]), 
            mode = "markers", showlegend=FALSE, hoverinfo = 'text',
            text = text) %>%
      layout(title=g, plot_bgcolor= "black", paper_bgcolor = "black", scene = list(xaxis = list(showspikes=FALSE), 
                                                                                   yaxis = list(showspikes=FALSE), 
                                                                                   zaxis = list(showspikes=FALSE)))
  } else if((input$plot=="UMAP") & (input$dim=="2D") & (input$lay=="Gene Expression")) {
    df <- data.frame(xumap2d, Expression = Expressioni)
    plot_ly(data = df, x = ~UMAP_1, y = ~UMAP_2, color = ~Expression, source = "exp",
            colors = c("firebrick4", "yellow"),  sizes=c(4,4), #size=4, 
            type = "scattergl", mode="markers",
            height = as.numeric(input$dimension[2]),
            hoverinfo = 'text',
            text = text) %>% 
      layout(title=g, plot_bgcolor= "black", paper_bgcolor = "black")
  } else if((input$plot=="tSNE") & (input$dim=="3D") & (input$lay=="Clusters")) {
    df <- data.frame(xtsne3d, gCol = gColi, size = sizei)
    # text = paste('<span style="font-size: 120%;">Age:', cell$age, '<br>Gender:', cell$gender, '<br>Genotype:', cell$genotype, '</span>')
    t <- list(
      family = "sans serif",
      size = 14,
      color = toRGB("grey50"))
    plot_ly(data = df, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, source = "clust",
            sizes=c(5,50), size=~size, color = ~gCol,
            marker = list(symbol = 'circle', sizemode = 'diameter'),
            colors = colrs, type = "scatter3d", 
            height = as.numeric(input$dimension[2]), 
            mode = "markers", showlegend=T, hoverinfo = 'text',
            text = text
    ) %>% 
      layout(title=g, 
             paper_bgcolor = 'rgb(243, 243, 243)', 
             plot_bgcolor = 'rgb(243, 243, 243)', 
             scene = list(xaxis = list(showspikes=FALSE), 
                          yaxis = list(showspikes=FALSE),
                          zaxis = list(showspikes=FALSE),
                          camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25))
             )
      )
  } else if((input$plot=="tSNE") & (input$dim=="2D") & (input$lay=="Clusters")) {
    df <- data.frame(xtsne2d, gCol = gColi, size = sizei)
    
    # text = paste('<span style="font-size: 120%;">Age:', cell$age, '<br>Gender:', cell$gender, '<br>Genotype:', cell$genotype, '</span>')
    plot_ly(data = df, x = ~tSNE_1, y = ~tSNE_2, source = "clust",
            color = ~gCol,
            colors = colrs,
            sizes=c(5,50), size=~size,  marker = list(symbol = 'circle', sizemode = 'diameter'),
            type = "scattergl", mode="markers", showlegend=T, 
            height = as.numeric(input$dimension[2]),
            hoverinfo = 'text',
            text = text
    ) %>% 
      layout(title=g, paper_bgcolor = 'rgb(243, 243, 243)', plot_bgcolor = 'rgb(243, 243, 243)')
  }else if((input$plot=="tSNE") & (input$dim=="3D") & (input$lay=="Gene Expression")) {
    df <- data.frame(xtsne3d, Expression = Expressioni)
    plot_ly(data = df, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, source = "exp",
            sizes=c(30,30), size=30, color = ~Expression,
            colors = c("firebrick4", "yellow"), type = "scatter3d", 
            height = as.numeric(input$dimension[2]), 
            mode = "markers", showlegend=FALSE, hoverinfo = 'text',
            text = text) %>%
      layout(title=g, plot_bgcolor= "black", paper_bgcolor = "black", scene = list(xaxis = list(showspikes=FALSE), 
                                                                                   yaxis = list(showspikes=FALSE), 
                                                                                   zaxis = list(showspikes=FALSE)))
  } else if((input$plot=="tSNE") & (input$dim=="2D") & (input$lay=="Gene Expression")) {
    df <- data.frame(xtsne2d, Expression = Expressioni)
    plot_ly(data = df, x = ~tSNE_1, y = ~tSNE_2, color = ~Expression, source = "exp",
            colors = c("firebrick4", "yellow"),  sizes=c(4,4), #size=4, 
            type = "scattergl", mode="markers",
            height = as.numeric(input$dimension[2]),
            hoverinfo = 'text',
            text = text) %>% 
      layout(title=g, plot_bgcolor= "black", paper_bgcolor = "black")
  }
  
  })
  
  #observe({updateTabsetPanel(session, inputId = "tabs1", selected = "dataset1")})
  
  output$tsne.plot <- renderPlotly({
      render.d()
    })
  #as.numeric(strsplit(rownames(visp_alm_odds.ratio)[1], "_")[[1]][3])+1
  text1 <- function(ratio){
    dend<-as.dendrogram(hclust(dist(ratio)))
    return(dend2 <- seriate_dendrogram(dend, dist(ratio, method = "euclidean")))
    #output$text1 <- renderText({rownames(visp_alm_odds.ratio)[(click_event$pointNumber[[1]][1]-nrow(visp_alm_odds.ratio))+2]})
  }
  
  # observe({
  # #req(input$go)
  # click_event <- event_data("plotly_click", source = "A")
  # if(is.null(click_event)){
  #   output$text1 <-  renderText({""})
  # }else{
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. VC")){
  #     dend2 <- text1(mtg_v1_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg_v1_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_v1_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. ACC")){
  #     dend2 <- text1(mtg_acc_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_acc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     dend2 <- text1(mtg_ec_exc_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_ec_exc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     dend2 <- text1(mtg_ec_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_ec_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. VC")){
  #     dend2 <- text1(mtg1_visp_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_visp_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(mtg1_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(mtg1_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(mtg1_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(mtg1_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(mtg1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_mopn_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen H. ACC")){
  #     dend2 <- text1(v1_acc_odds.ratio)
  #     output$text1 <- renderText({rownames(v1_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_acc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     dend2 <- text1(v1_ec_exc_odds.ratio)
  #     output$text1 <- renderText({rownames(v1_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_ec_exc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     dend2 <- text1(v1_ec_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(v1_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_ec_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. VC")){
  #     dend2 <- text1(v11_visp_odds.ratio)
  #     output$text1 <- renderText({rownames(v11_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_visp_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(v11_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(v11_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(v11_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(v11_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(v11_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(v11_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(v11_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(v11_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_mopn_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #     dend2 <- text1(acc_ec_exc_odds.ratio)
  #     output$text1 <- renderText({rownames(acc_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc_ec_exc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     dend2 <- text1(acc_ec_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(acc_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc_ec_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. VC")){
  #     dend2 <- text1(acc1_visp_odds.ratio)
  #     output$text1 <- renderText({rownames(acc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_visp_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(acc1_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(acc1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(acc1_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(acc1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(acc1_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(acc1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(acc1_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(acc1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_mopn_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #     dend2 <- text1(ec_exc_ec_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_ec_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #     dend2 <- text1(ec_exc_sfg_exc_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc_sfg_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_sfg_exc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     dend2 <- text1(ec_exc_sfg_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_sfg_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. VC")){
  #     dend2 <- text1(ec_exc1_visp_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_visp_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(ec_exc1_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(ec_exc1_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(ec_exc1_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(ec_exc1_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_exc1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_mopn_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #     dend2 <- text1(ec_inh_sfg_exc_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh_sfg_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh_sfg_exc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     dend2 <- text1(ec_inh_sfg_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh_sfg_inh_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. VC")){
  #     dend2 <- text1(ec_inh1_visp_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_visp_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(ec_inh1_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(ec_inh1_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(ec_inh1_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(ec_inh1_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(ec_inh1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_mopn_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Leng et al. H. SFG_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #     dend2 <- text1(sfg_exc_sfg_inh_odds.ratio)
  #     output$text1 <- renderText({rownames(sfg_exc_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(sfg_exc_sfg_inh_odds.ratio))+1]})
  #   }
    
  #   if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ALM")){
  #     dend2 <- text1(visp_alm_odds.ratio)
  #     output$text1 <- renderText({rownames(visp_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_alm_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(visp_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(visp_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(visp_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(visp_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(visp_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(visp_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_mopn_odds.ratio))+1]})
  #   }
    
    
  #   if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. ACA")){
  #     dend2 <- text1(alm_aca_odds.ratio)
  #     output$text1 <- renderText({rownames(alm_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_aca_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(alm_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(alm_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(alm_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(alm_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_mopn_odds.ratio))+1]})
  #   }
    
    
  #   if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpC")){
  #     dend2 <- text1(aca_mopc_odds.ratio)
  #     output$text1 <- renderText({rownames(aca_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopc_odds.ratio))+1]})
  #   }
  #   if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(aca_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(aca_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopn_odds.ratio))+1]})
  #   }
    
  #   if((input$dataset == "Allen M. MOpC")&(input$dataseto == "Allen M. MOpN")){
  #     dend2 <- text1(mopc_mopn_odds.ratio)
  #     output$text1 <- renderText({rownames(mopc_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mopc_mopn_odds.ratio))+1]})
  #   }
  #   # else if((input$dataset == "Allen H. VC")){
  #   #   dend<-as.dendrogram(hclust(dist(v1_acc_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(v1_acc_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(v1_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_acc_odds.ratio))+1]})
  #   # }
  #   # else if((input$dataset == "Allen H. ACC")){
  #   #   dend<-as.dendrogram(hclust(dist(acc1_visp_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(acc1_visp_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(acc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_visp_odds.ratio))+1]})
  #   # }
  #   # else if((input$dataset == "Leng et al. H. EC_Exc")){
  #   #   dend<-as.dendrogram(hclust(dist(ec_exc1_visp_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(ec_exc1_visp_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(ec_exc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_visp_odds.ratio))+1]})
  #   # }
  #   # else if((input$dataset == "Leng et al. H. EC_Inh")){
  #   #   dend<-as.dendrogram(hclust(dist(ec_inh1_visp_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(ec_inh1_visp_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(ec_inh1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_visp_odds.ratio))+1]})
  #   # }
  #   #else if((input$dataset == "Allen M. VC")){
  #     #dend<-as.dendrogram(hclust(dist(visp_alm_odds.ratio)))
  #     #dend2 <- seriate_dendrogram(dend, dist(visp_alm_odds.ratio, method = "euclidean"))
  #     ##output$text1 <- renderText({rownames(visp_alm_odds.ratio)[(click_event$pointNumber[[1]][1]-nrow(visp_alm_odds.ratio))+2]})
  #     #output$text1 <- renderText({rownames(visp_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_alm_odds.ratio))+1]})
  #   #}
  #   # else if((input$dataset == "Allen M. ALM")){
  #   #   dend<-as.dendrogram(hclust(dist(alm_aca_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(alm_aca_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(alm_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_aca_odds.ratio))+1]})
  #   # }
  #   # else if((input$dataset == "Allen M. ACA")){
  #   #   dend<-as.dendrogram(hclust(dist(aca_mopc_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(aca_mopc_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(aca_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopc_odds.ratio))+1]})
  #   # }
  #   # else if((input$dataset == "Allen M. MOpC")){
  #   #   dend<-as.dendrogram(hclust(dist(mopc_mopn_odds.ratio)))
  #   #   dend2 <- seriate_dendrogram(dend, dist(mopc_mopn_odds.ratio, method = "euclidean"))
  #   #   output$text1 <- renderText({rownames(mopc_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mopc_mopn_odds.ratio))+1]})
  #   # }
  #   #print(click_event$y-nrow(visp_alm_odds.ratio))
  #   #print({rownames(visp_alm_odds.ratio)[(click_event$y-nrow(visp_alm_odds.ratio))+1]})
  # }
  # })

  toptext1 <- function(ratio){
    dend<-as.dendrogram(hclust(dist(ratio)))
    return(dend2 <- seriate_dendrogram(dend, dist(ratio, method = "euclidean")))
  }

  # output$heat.click <- renderTable({
  #   #req(input$go)
  #   click_event <- event_data("plotly_click", source = "A")
  #   req(click_event)
  #   if(is.null(click_event)){
  #     "Click on cells in the plot to view genes belonging to a particular cluster" 
  #   }else{
  #      #output$text2 <- renderText({paste("You have selected", input$var)
  #     #print(as.numeric(strsplit((rownames(visp_alm_odds.ratio)[(click_event$pointNumber[[1]][1]-nrow(visp_alm_odds.ratio))+2]), "_")[[1]][3]))
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. VC")){
  #       dend2 <- toptext1(mtg_v1_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg_v1_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_v1_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- toptext1(mtg_acc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_acc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext1(mtg_ec_exc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_ec_exc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext1(mtg_ec_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg_ec_inh_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext1(mtg1_visp_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_visp_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(mtg1_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_alm_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(mtg1_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_aca_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(mtg1_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(mtg1_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mtg1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mtg1_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MTG_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- toptext1(v1_acc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v1_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_acc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext1(v1_ec_exc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v1_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_ec_exc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext1(v1_ec_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v1_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_ec_inh_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext1(v11_visp_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_visp_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(v11_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_alm_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(v11_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_aca_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(v11_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(v11_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext1(acc_ec_exc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc_ec_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc_ec_exc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext1(acc_ec_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc_ec_inh_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext1(v11_visp_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(v11_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v11_visp_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(acc1_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_alm_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(acc1_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_aca_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(acc1_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(acc1_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(acc1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext1(ec_exc_ec_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc_ec_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_ec_inh_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- toptext1(ec_exc_sfg_exc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc_sfg_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_sfg_exc_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext1(ec_exc_sfg_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc_sfg_inh_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext1(ec_exc1_visp_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_visp_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(ec_exc1_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_alm_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(ec_exc1_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_aca_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(ec_exc1_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_mopc_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(ec_exc1_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_exc1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_mopn_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- toptext1(ec_inh_sfg_exc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh_sfg_exc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh_sfg_exc_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext1(ec_inh_sfg_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh_sfg_inh_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext1(ec_inh1_visp_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_visp_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(ec_inh1_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh1_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_alm_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(ec_inh1_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh1_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_aca_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(ec_inh1_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh1_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_mopc_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(ec_inh1_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(ec_inh1_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_mopn_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. SFG_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext1(sfg_exc_sfg_inh_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(sfg_exc_sfg_inh_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(sfg_exc_sfg_inh_odds.ratio))+1]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Exc")
  #     }
      
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext1(visp_alm_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(visp_alm_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_alm_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(visp_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(visp_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_aca_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(visp_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(visp_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(visp_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(visp_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(visp_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
      
      
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext1(alm_aca_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(alm_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_aca_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(alm_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(alm_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(alm_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(alm_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
      
      
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext1(aca_mopc_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(aca_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopc_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(aca_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(aca_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
      
      
  #     if((input$dataset == "Allen M. MOpC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext1(mopc_mopn_odds.ratio)
  #       d <- as.numeric(strsplit((rownames(mopc_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mopc_mopn_odds.ratio))+1]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     # else if((input$dataset == "Allen H. VC")){
  #     #   dend<-as.dendrogram(hclust(dist(v1_acc_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(v1_acc_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(v1_acc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(v1_acc_odds.ratio))+1]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     # }
  #     # else if((input$dataset == "Allen H. ACC")){
  #     #   dend<-as.dendrogram(hclust(dist(acc1_visp_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(acc1_visp_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(acc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(acc1_visp_odds.ratio))+1]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     # }
  #     # else if((input$dataset == "Leng et al. H. EC_Exc")){
  #     #   dend<-as.dendrogram(hclust(dist(ec_exc1_visp_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(ec_exc1_visp_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(ec_exc1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_exc1_visp_odds.ratio))+1]), "_")[[1]][4])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     # }
  #     # else if((input$dataset == "Leng et al. H. EC_Inh")){
  #     #   dend<-as.dendrogram(hclust(dist(ec_inh1_visp_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(ec_inh1_visp_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(ec_inh1_visp_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(ec_inh1_visp_odds.ratio))+1]), "_")[[1]][4])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     # }
  #     #else if((input$dataset == "Allen M. VC")){
  #     #else if((input$dataseto == "Leng et al. H. EC_Exc")&(input$dataset == "Allen M. VC")){
  #       #dend<-as.dendrogram(hclust(dist(visp_alm_odds.ratio)))
  #       #dend2 <- seriate_dendrogram(dend, dist(visp_alm_odds.ratio, method = "euclidean"))
  #       #dend2<- toptext1(t(ec_exc1_visp_odds.ratio))
  #       #d <- as.numeric(strsplit((rownames(t(ec_exc1_visp_odds.ratio)[order.dendrogram(dend2),])[abs(click_event$y-nrow(t(ec_exc1_visp_odds.ratio)))+1]), "_")[[1]][3])
  #       #topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     #}
  #     # else if((input$dataset == "Allen M. ALM")){
  #     #   dend<-as.dendrogram(hclust(dist(alm_aca_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(alm_aca_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(alm_aca_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(alm_aca_odds.ratio))+1]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     # }
  #     # else if((input$dataset == "Allen M. ACA")){
  #     #   dend<-as.dendrogram(hclust(dist(aca_mopc_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(aca_mopc_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(aca_mopc_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(aca_mopc_odds.ratio))+1]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     # }
  #     # else if((input$dataset == "Allen M. MOpC")){
  #     #   dend<-as.dendrogram(hclust(dist(mopc_mopn_odds.ratio)))
  #     #   dend2 <- seriate_dendrogram(dend, dist(mopc_mopn_odds.ratio, method = "euclidean"))
  #     #   d <- as.numeric(strsplit((rownames(mopc_mopn_odds.ratio[order.dendrogram(dend2),])[abs(click_event$y-nrow(mopc_mopn_odds.ratio))+1]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     # }

  #     topgenes <- as.matrix(setNames(topgenes_collection$find('{}', fields = '{"_id":0, "Gene":1, "LogFC":1, "Pval":1}')[1:20,],c("Target","LogFC","Pval")))
  #     topgenes
  #   }
  #   #print(rownames(visp_alm_odds.ratio)[(29-nrow(visp_alm_odds.ratio))+2])
  #   #topgenes_collection1 = mongolite::mongo(collection = paste0("topgenes", click_event$curveNumber+1), db = "allen_VISp_eis4")
  #   #topgenes_collection2 = mongolite::mongo(collection = paste0("topgenes", click_event$curveNumber+1), db = "allen_ALM_eis4")
  # },striped = TRUE)

  text2 <- function(ratio){
    dst <- dist(ratio)
    hc_row <- hclust(dst)
    col_dend <- as.dendrogram(hc_row)
    dend2 <- seriate_dendrogram(col_dend, dst)
    return(dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))])))
  }

  # observe({
  #   #req(input$go)
  #   click_event <- event_data("plotly_click", source = "A")
  #   if(is.null(click_event)){
  #     output$text2 <-  renderText({""})
  #   }else{
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. VC")){
  #       dend2 <- text2(t(mtg_v1_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg_v1_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- text2(t(mtg_acc_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- text2(t(mtg_ec_exc_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- text2(t(mtg_ec_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- text2(t(mtg1_visp_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(mtg1_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(mtg1_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(mtg1_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(mtg1_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(mtg1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- text2(t(v1_acc_odds.ratio))
  #       output$text2 <- renderText({colnames(v1_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- text2(t(v1_ec_exc_odds.ratio))
  #       output$text2 <- renderText({colnames(v1_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- text2(t(v1_ec_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(v1_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- text2(t(v11_visp_odds.ratio))
  #       output$text2 <- renderText({colnames(v11_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(v11_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(v11_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(v11_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(v11_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(v11_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(v11_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(v11_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(v11_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- text2(t(acc_ec_exc_odds.ratio))
  #       output$text2 <- renderText({colnames(acc_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- text2(t(acc_ec_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(acc_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- text2(t(acc1_visp_odds.ratio))
  #       output$text2 <- renderText({colnames(acc1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(acc1_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(acc1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(acc1_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(acc1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(acc1_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(acc1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(acc1_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(acc1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- text2(t(ec_exc_ec_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- text2(t(ec_exc_sfg_exc_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc_sfg_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- text2(t(ec_exc_sfg_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- text2(t(ec_exc1_visp_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(ec_exc1_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(ec_exc1_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(ec_exc1_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(ec_exc1_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_exc1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- text2(t(ec_inh_sfg_exc_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh_sfg_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- text2(t(ec_inh_sfg_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- text2(t(ec_inh1_visp_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(ec_inh1_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(ec_inh1_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(ec_inh1_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(ec_inh1_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(ec_inh1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Leng et al. H. SFG_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- text2(t(sfg_exc_sfg_inh_odds.ratio))
  #       output$text2 <- renderText({colnames(sfg_exc_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
      
      
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- text2(t(visp_alm_odds.ratio))
  #       output$text2 <- renderText({colnames(visp_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(visp_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(visp_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(visp_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(visp_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(visp_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(visp_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
      
      
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- text2(t(alm_aca_odds.ratio))
  #       output$text2 <- renderText({colnames(alm_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(alm_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(alm_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(alm_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(alm_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
      
      
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- text2(t(aca_mopc_odds.ratio))
  #       output$text2 <- renderText({colnames(aca_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(aca_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(aca_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
      
      
  #     if((input$dataset == "Allen M. MOpC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- text2(t(mopc_mopn_odds.ratio))
  #       output$text2 <- renderText({colnames(mopc_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     }
  #     # else if((input$dataseto == "Allen H. ACC")){
  #     #   dst <- dist(t(v1_acc_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(v1_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else 
  #     #if((input$dataseto == "Leng et al. H. EC_Exc")&(input$dataset == "Allen M. VC")){
  #       #dend2 <- text2((ec_exc1_visp_odds.ratio))
  #     #else if((input$dataseto == "Leng et al. H. EC_Exc")){
  #       #dst <- dist(t(mtg_ec_exc_odds.ratio))
  #       #hc_row <- hclust(dst)
  #       #col_dend <- as.dendrogram(hc_row)
  #       #dend2 <- seriate_dendrogram(col_dend, dst)
  #       #dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #       #output$text2 <- renderText({colnames(mtg_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #       #output$text2 <- renderText({colnames(t(ec_exc1_visp_odds.ratio)[,order.dendrogram(dend2)])[click_event$x]})
  #     #}
  #     # else if((input$dataseto == "Leng et al. H. EC_Inh")){
  #     #   dst <- dist(t(mtg_ec_inh_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(mtg_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else if((input$dataseto == "Allen M. VC")){
  #     #   dst <- dist(t(acc1_visp_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(acc1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else if((input$dataseto == "Allen M. ALM")){
  #     #   dst <- dist(t(visp_alm_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(visp_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else if((input$dataseto == "Allen M. ACA")){
  #     #   dst <- dist(t(alm_aca_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(alm_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else if((input$dataseto == "Allen M. MOpC")){
  #     #   dst <- dist(t(aca_mopc_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(aca_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     # else if((input$dataseto == "Allen M. MOpN")){
  #     #   dst <- dist(t(mopc_mopn_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   output$text2 <- renderText({colnames(mopc_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]})
  #     # }
  #     #print({colnames(visp_alm_odds.ratio)[click_event$x]})
  #     #renderText({colnames(visp_alm_odds.ratio)[click_event$pointNumber[[1]][2]+1]})
  #   }
    
  # })

  toptext2 <- function(ratio){
    dst <- dist(ratio)
    hc_row <- hclust(dst)
    col_dend <- as.dendrogram(hc_row)
    dend2 <- seriate_dendrogram(col_dend, dst)
    return(dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))])))
  }

  # output$heat2.click <- renderTable({
  #   #req(input$go)
  #   click_event <- event_data("plotly_click", source = "A")
  #   req(click_event)
  #   if(is.null(click_event)){
  #     "Click on cells in the plot to view genes belonging to a particular cluster" 
  #   }else{
  #     #output$text2 <- renderText({paste("You have selected", input$var)
  #     #print(click_event)
  #     #print(click_event$pointNumber[[1]][1])
  #     #print(click_event$pointNumber[[1]][2])
  #     #print(colnames(visp_alm_odds.ratio)[click_event$pointNumber[[1]][2]+1])
  #     #print(as.numeric(strsplit((colnames(visp_alm_odds.ratio)[click_event$pointNumber[[1]][2]+1]), "_")[[1]][3]))
  #     #d <- as.numeric(strsplit((colnames(visp_alm_odds.ratio)[click_event$pointNumber[[1]][2]+1]), "_")[[1]][3])
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. VC")){
  #       dend2 <- toptext2(t(mtg_v1_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg_v1_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_V1_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- toptext2(t(mtg_acc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext2(t(mtg_ec_exc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext2(t(mtg_ec_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext2(t(mtg1_visp_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(mtg1_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(mtg1_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(mtg1_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen H. MTG")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(mtg1_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mtg1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen H. ACC")){
  #       dend2 <- toptext2(t(v1_acc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v1_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext2(t(v1_ec_exc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v1_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext2(t(v1_ec_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v1_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext2(t(v11_visp_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v11_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(v11_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v11_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(v11_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v11_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(v11_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v11_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen H. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(v11_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(v11_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Exc")){
  #       dend2 <- toptext2(t(acc_ec_exc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext2(t(acc_ec_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext2(t(acc1_visp_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(acc1_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(acc1_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(acc1_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen H. ACC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(acc1_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(acc1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. EC_Inh")){
  #       dend2 <- toptext2(t(ec_exc_ec_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- toptext2(t(ec_exc_sfg_exc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc_sfg_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext2(t(ec_exc_sfg_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext2(t(ec_exc1_visp_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(ec_exc1_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(ec_exc1_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(ec_exc1_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Exc")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(ec_exc1_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_exc1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Exc")){
  #       dend2 <- toptext2(t(ec_inh_sfg_exc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh_sfg_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Exc")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext2(t(ec_inh_sfg_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Inh")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. VC")){
  #       dend2 <- toptext2(t(ec_inh1_visp_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh1_visp_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_VISp_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(ec_inh1_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh1_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(ec_inh1_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh1_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(ec_inh1_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh1_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. EC_Inh")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(ec_inh1_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(ec_inh1_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     if((input$dataset == "Leng et al. H. SFG_Exc")&(input$dataseto == "Leng et al. H. SFG_Inh")){
  #       dend2 <- toptext2(t(sfg_exc_sfg_inh_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(sfg_exc_sfg_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "SFG_Inh")
  #     }
      
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ALM")){
  #       dend2 <- toptext2(t(visp_alm_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(visp_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(visp_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(visp_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(visp_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(visp_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen M. VC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(visp_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(visp_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
      
      
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. ACA")){
  #       dend2 <- toptext2(t(alm_aca_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(alm_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(alm_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(alm_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen M. ALM")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(alm_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(alm_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
      
      
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpC")){
  #       dend2 <- toptext2(t(aca_mopc_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(aca_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     }
  #     if((input$dataset == "Allen M. ACA")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(aca_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(aca_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
      
      
  #     if((input$dataset == "Allen M. MOpC")&(input$dataseto == "Allen M. MOpN")){
  #       dend2 <- toptext2(t(mopc_mopn_odds.ratio))
  #       d <- as.numeric(strsplit((colnames(mopc_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #       topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     }
  #     # if(input$dataseto == "Allen H. ACC"){
  #     #   dst <- dist(t(v1_acc_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(v1_acc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACC_eis4")
  #     # }
  #     #if((input$dataseto == "Leng et al. H. EC_Exc")&(input$dataset == "Allen M. VC")){
  #       #dst <- dist(t(acc1_visp_odds.ratio))
  #       #hc_row <- hclust(dst)
  #       #col_dend <- as.dendrogram(hc_row)
  #       #dend2 <- seriate_dendrogram(col_dend, dst)
  #       #dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #       #dend2 <- toptext2((ec_exc1_visp_odds.ratio))
  #       #d <- as.numeric(strsplit((colnames(t(ec_exc1_visp_odds.ratio)[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #       #topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     #}
  #     # if((input$dataseto == "Leng et al. H. EC_Exc")){
  #     #   dst <- dist(t(acc_ec_exc_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(acc_ec_exc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Exc")
  #     # }
  #     # if((input$dataseto == "Leng et al. H. EC_Inh")){
  #     #   dst <- dist(t(acc_ec_inh_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(acc_ec_inh_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][4])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "EC_Inh")
  #     # }
  #     # if(input$dataseto == "Allen M. ALM"){
  #     #   dst <- dist(t(visp_alm_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(visp_alm_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ALM_eis4")
  #     # }
  #     # if(input$dataseto == "Allen M. ACA"){
  #     #   dst <- dist(t(alm_aca_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(alm_aca_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_ACA_eis4")
  #     # }
  #     # if(input$dataseto == "Allen M. MOpC"){
  #     #   dst <- dist(t(aca_mopc_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(aca_mopc_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpC_eis4")
  #     # }
  #     # if(input$dataseto == "Allen M. MOpN"){
  #     #   dst <- dist(t(mopc_mopn_odds.ratio))
  #     #   hc_row <- hclust(dst)
  #     #   col_dend <- as.dendrogram(hc_row)
  #     #   dend2 <- seriate_dendrogram(col_dend, dst)
  #     #   dend2 <- rotate(dend2, order = rev(labels(dst)[get_order(as.hclust(dend2))]))
  #     #   d <- as.numeric(strsplit((colnames(mopc_mopn_odds.ratio[,order.dendrogram(dend2)])[click_event$x]), "_")[[1]][3])
  #     #   topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d+1), db = "allen_MOpN_eis4")
  #     # }
  #     topgenes <- as.matrix(setNames(topgenes_collection$find('{}', fields = '{"_id":0, "Gene":1, "LogFC":1, "Pval":1}')[1:20,],c("Target","LogFC","Pval")))
  #     topgenes
  #   }
  #   #print(rownames(visp_alm_odds.ratio)[(29-nrow(visp_alm_odds.ratio))+2])
  #   #topgenes_collection1 = mongolite::mongo(collection = paste0("topgenes", click_event$curveNumber+1), db = "allen_VISp_eis4")
  #   #topgenes_collection2 = mongolite::mongo(collection = paste0("topgenes", click_event$curveNumber+1), db = "allen_ALM_eis4")
  # },striped = TRUE)

  #topgenes_collection = mongolite::mongo(collection = paste0("topgenes", 1), db = "allen")
  #dim(topgenes_collection$find('{}', fields = '{"_id":0, "Gene":1, "LogFC":1, "Pval":1}')[1:5,c(1,2)])
  
  output$tsne.click <- renderTable({
    d <- event_data("plotly_click", source = "clust")
    req(d)
    if(is.null(d)){
      "Click on cells in the plot to view genes belonging to a particular cluster" 
    }
    else {
      #if(d$curveNumber){
          if(input$dataset == "iPSC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "iPSC", url = mongodb_url)
          }
          if(input$dataset == "H. EC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "EC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Exc"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "H. EC_Exc_0.7"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes_0.7", d$curveNumber+1), db = "EC_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. EC_Inh"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "EC_Inh", url = mongodb_url)
          }
          if(input$dataset == "H. SFG"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "SFG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Exc"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "SFG_Exc", url = mongodb_url)
          }
          if(input$dataset == "Leng et al. H. SFG_Inh"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "SFG_Inh", url = mongodb_url)
          }
		      if(input$dataset == "Allen H. MTG"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_MTG_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. VC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_V1_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen H. ACC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_ACC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. VC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_VISp_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ALM"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_ALM_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. ACA"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_ACA_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpC"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_MOpC_eis4", url = mongodb_url)
          }
          if(input$dataset == "Allen M. MOpN"){
            topgenes_collection = mongolite::mongo(collection = paste0("topgenes", d$curveNumber+1), db = "allen_MOpN_eis4", url = mongodb_url)
          }
        topgenes <- as.matrix(setNames(topgenes_collection$find('{}', fields = '{"_id":0, "Gene":1, "LogFC":1, "Pval":1}')[1:20,],c("Target","LogFC","Pval")))
        topgenes
      #}
    }
      
      
    })
  
  if (!file.exists("counter.Rdata")) {
    counter <- 0
  }else
  {
    load(file="counter.Rdata")
  }
  counter  <- counter + 1
  save(counter, file="counter.Rdata")
  output$notificationMenu <- renderMenu({
    dropdownMenu(type = "notifications", notificationItem(
      text = paste(counter,"users so far"),
      icon("users")
    ))
  })
  output$counter <- 
    renderText({
      paste("You are visitor no: ", counter)
    })
  
  dataset <- c("Human Middle Temporal Gyrus (MTG)", "Human Primary Visual Cortex (V1)", "Human Anterior Cingulate 
            Cortex (ACC)", "Mouse Primary Visual Cortex (VISp)", "Mouse Anterior Lateral Motor Area 
(ALM)", "Mouse Anterior Cingulate Area (ACA)", "Mouse Primary Motor Area (MOp)", "Mouse Primary 
Motor Area (MOp)")
  nuclei <- c("15,928 nuclei", "8,998 nuclei", "7,283 nuclei", "15,413 cells", "10,068 cells", "5,028 
cells", "4,916 cells", "6,847 nuclei")
  data_summary <- data.frame("Dataset"=dataset, "Cells|Nuclei"=nuclei)
  output$data_summary <- renderTable({
    data_summary
  })
#}

