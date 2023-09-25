library(viridis)
library(RSQLite)
library(manipulateWidget)
library(highcharter)
library(tm)
library(reshape2)
library(plyr)
library(htmltools)
library(RColorBrewer)
library(shiny)
library(shinyBS)

source(file = 'functions_helpers.R')


shinyInput <- function(FUN, len, id, ...)
{
  inputs <- character(len)
  for (i in seq_len(len))
  {
    inputs[i] <- as.character(FUN(paste0(id, i), style='padding:4px;font-size:80%', ...))
  }
  inputs
}



shinyServer(

  function(input, output, session)
  {


# DATA IMPORT -------------------------------------------------------------

    #Database for the pvalue = f(age) for each gene, tissue and variable
    # DBconnection3 <- dbConnect(RSQLite::SQLite(),
    #                            dbname = "DB_pvalue_ShARP-LM.db")


    #Significance pvalue = f(age) for a single gene across tissue
    #Need to take gene-specific information from DB_pvalue_ShARP-LM.db
    #p_Alterations gathers signficiance for all tissue for a single genes and vriable
    p_Alterations <- reactive({
       
      gene <- input$gene
      variable <- input$variable_Gene_Alteration
      validate(need(!is.null(gene) & gene != "", "")) 
      tissue <- geneList$tissue[[match(gene, geneList$gene)]]
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)

      p <- data.frame()
      for (i in 1:length(modifiedTissuename))
      {
        if (tissue[i] %in% c("Ovary", "Uterus", "Vagina", "Prostate", "Fallopian Tube", "Testis") & variable %in% c("Gender", "Int_Age_Gender"))
        {
          tmp <- NULL
        } else
        {
          # #SQL database
          # tmp <- dbGetQuery(DBconnection3, paste("SELECT * FROM", paste0(modifiedTissuename[i], "_", variable)))
          # #RDS files local
          # tmp <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_pvalue_ShARP-LM/",
          #                      modifiedTissuename[i], "_", variable, ".RDS", sep = ""))
          #RDS files online
          tmp <- readRDS(paste("data/DB_pvalue_ShARP-LM/",
                               modifiedTissuename[i], "_",
                               variable, ".RDS", sep = ""))
          tmp <- tmp[tmp$gene == gene,]
          tmp <- tmp[, -1]
          tmp <- -log10(tmp)
          
    
          
          fit <- smooth.spline(as.numeric(colnames(tmp)), tmp)
          fit <- predict(fit, seq(26, 64, by = 0.5))$y
          tmp <- data.frame(t(fit))
          colnames(tmp) <- seq(26, 64, by = 0.5)
          rownames(tmp) <- tissue[i]
        }
        p <- rbind(p, tmp)
      }
      p
    })


    #alterations significance for a single gene, tissue and variable
    p_Alteration_gene <- reactive({
      gene <- input$gene
      variable <- input$variable_Gene_Alteration
      tissue <- input$tissue
      validate(need(!is.null(tissue) & tissue != "" & tissue != "All tissues", ""))
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)

      if (tissue %in% c("Ovary", "Uterus", "Vagina", "Prostate", "Fallopian Tube", "Testis") & variable %in% c("Gender", "Int_Age_Gender"))
      {
        pvalueData <- NULL
      } else
      {
        # #SQL database
        # pvalueData <- dbGetQuery(DBconnection3, paste("SELECT * FROM", paste0(modifiedTissuename, "_", variable)))
        # #RDS files local
        # pvalueData <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_pvalue_ShARP-LM/",
        #                      modifiedTissuename, "_", variable, ".RDS", sep = ""))
        #RDS files online
        pvalueData <- readRDS(paste("data/DB_pvalue_ShARP-LM/",
                                    modifiedTissuename, "_", variable, ".RDS",
                                    sep = ""))
        pvalueData <- pvalueData[pvalueData$gene == gene,]
        pvalueData <- pvalueData[, -1]
      }
      pvalueData
    })

    #p gather significance for all genes given a tissue and a variable
    p <- reactive({
      validate(need(input$tissue_2 != "All tissues" & input$tissue_2 != "", ""))
      if (input$tissue_2 %in% c("Ovary", "Uterus", "Vagina", "Prostate", "Fallopian Tube", "Testis") & input$variable_model_2 %in% c("Gender", "Int_Age_Gender"))
      {
        p <- NULL
      } else
      {
        modifiedTissuename <- gsub("-", "_", input$tissue_2)
        modifiedTissuename <- gsub(" ", "", modifiedTissuename)
        modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
        modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
        # #SQL database
        # p <- dbGetQuery(DBconnection3, paste("SELECT * FROM", paste0(modifiedTissuename, "_", input$variable_model_2)))
        # #RDS files local
        # p <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_pvalue_ShARP-LM/",
        #                    modifiedTissuename, "_", input$variable_model_2, ".RDS", sep = ""))
        #RDS files online
        p <- readRDS(paste("data/DB_pvalue_ShARP-LM/",
                           modifiedTissuename, "_", input$variable_model_2,
                           ".RDS", sep = ""))
 
        rownames(p) <- p$gene
        p <- p[,-1]
      }
      p
    })

    peak <- readRDS("data/peaks_plot_2.RDS") #to knwo the FDR for each peak for FDRvsAge
    AgeWavesList <- readRDS("data/FDR_Interpolation_Tissue.RDS")
    geneList <- readRDS("data/DT_NbTissuePerGene.RDS")
    donorCondition <- readRDS("data/donorCondition.RDS")
    conditionID <- readRDS("data/conditionID.RDS")
    technicalCondition <- readRDS("data/technicalCondition.RDS")
    Cluster_Reactome <- readRDS("data/Cluster_Reactome.RDS")
    Cluster_Reactome_affiliation <- readRDS("data/Reactome_cluster_Affiliation.RDS")
    moduleInfo <- readRDS("data/moduleInfo_Corrected_Atlas.RDS")
    # color2 <- c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    #             "#7FC97F", "#BEAED4", "#FDC086", "#E6AB02", "#386CB0", "#A6761D", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02",
    #             "#7570B3", "#E7298A", "#66A61E", "#FFFF99", "#BEBADA")
    #

    #Import of the gene information
    # DBconnection <- dbConnect(RSQLite::SQLite(),
    #                            dbname = "DB_GeneExpressionV2.db")
    geneData <- reactive({
      gene <- as.character(input$gene)
      gene <- gsub("-", "_", gene) # '-' creates issues when calling database so '-' were replaced by "_" in the database's names
      gene <- gsub("\\.", "_", gene) # '.' creates issues when calling database so '-' were replaced by "_" in the database's names

      # #SQL databases
      # a <- dbGetQuery(DBconnection, paste("SELECT * FROM", gene, sep = " "))
      #RDS files local
      # a <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_GeneExpressionV2/",
      #                    gene, ".RDS", sep = ""))
      #RDS files online
      a <- readRDS(paste("data/DB_GeneExpressionV2/",
                         gene, ".RDS", sep = ""))
      a
    })

    #Database for the peak enrichment data
    # DBconnection2 <- dbConnect(RSQLite::SQLite(),
    #                           dbname = "DB_PeakEnrichment.db")

# GENE SECTION ------------------------------------------------------------

##SIDEBAR
    #Update gene List and tissue List
    observe({
      updatePickerInput(session, "gene",
                        choices =  as.character(geneList$gene)[order(as.character(geneList$gene))],
                        selected = "CDKN2A")
    })
    observe({
      gene <- as.character(req(input$gene))

      #Keep the selected tissue if the new gene is expressed in it
      #Otherwise, put the heatmpa wil all tissues
 
      if (input$tissue %in% as.character(geneList$tissue[[match(gene, geneList$gene)]]))
      {
        selected <- input$tissue
      } else
      {
        selected <- "All tissues"
      }

      updatePickerInput(session, "tissue",
                        choices =  c("All tissues", as.character(geneList$tissue[[match(gene, geneList$gene)]])[order(as.character(geneList$tissue[[match(gene, geneList$gene)]]))]),
                        selected = selected)
      updateAwesomeRadio(session, "coloredBy", choices = list("All" = "none", "Sex" = "Sex"), selected = "none",
                         inline = T, checkbox = TRUE)
      updateAwesomeRadio(session, "shapedBy", choices = list("All" = "none", "Condition*" = "Condition"), selected = "none",
                         inline = T, checkbox = TRUE)

    })
    #Update condition list
    output$shapedByCondition <- renderUI({
      condition <- colnames(donorCondition)[-1][order(colnames(donorCondition)[-1])]
      subtext <- conditionID$ID[match(condition, conditionID$description)]
      pickerInput("condition", "Select:",
                  choices = condition,
                  choicesOpt = list(subtext = subtext),
                  options = list(`live-search` = TRUE))
    })

    #table Kruskal-Wallis test when condition is checked
    output$KruskalWallisCondition <- DT::renderDataTable({
      tissue <- input$tissue
      #when a new gene is selected, geneData is updated sooner than input$tissue generating an error message
      validate(need(nrow(geneData()[geneData()$tissue == input$tissue, ]) != 0, "")) #allow to hide the error message
      tmp <- geneData()[geneData()$tissue == tissue,]
      temp <- data.frame(SAMPID = tmp$samp_id, expression = tmp$expression)
      temp <- cbind(temp, donorCondition[match(temp$SAMPID, donorCondition$SAMPID), -1])
      tmp <- data.frame()
      for (i in 3:ncol(temp))
      {
        if (length(temp$expression[temp[,i] %in% c(0,1)])!=0 & length(unique(temp[,i][temp[,i] %in% c(0,1)])) == 2)
        {

          tmp <- rbind(tmp, data.frame(condition = colnames(temp)[i],
                                       pvalue = kruskal.test(temp$expression[temp[,i] %in% c(0,1)], temp[,i][temp[,i] %in% c(0,1)])$p.value))
        }

      }

      tmp$padj <- p.adjust(tmp$pvalue, method = "BH")
      tmp$conditionID <- conditionID$ID[match(tmp$condition, conditionID$description)]
      tmp$pvalue <- round(tmp$pvalue, 4)
      tmp$padj <- round(tmp$padj, 4)

      tmp <- tmp[order(tmp$padj),]

      condition <- ifelse(identical(as.character(input$shapedBy) == "Condition" & !is.na(input$condition), logical(0)) | as.character(input$shapedBy) == "none", "none", as.character(input$condition))

      validate(need(condition %in% colnames(donorCondition), ""))
      condition <- conditionID$ID[match(condition, conditionID$description)]
      DT::datatable(tmp[, c("conditionID", "pvalue", "padj", "condition")],
                   escape = F,
                   selection = list(mode = "none"),
                   colnames = c("Condition ID", "p-value", "adjusted p-value", "condition"),
                   rownames = F,
                   options = list(deferRender = TRUE, scrollY = 400, scroller = TRUE, scollX = T,
                                  rowCallback = JS(
                     "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
                     "var full_text = aData[3]",
                     "$('td:eq(0)', nRow).attr('title', full_text);",
                     "}"),
                     columnDefs = list(list(visible=FALSE, targets=c(3)),
                                       list(className = 'dt-center', targets = 1:3))
                     ),
                   extensions = "Scroller"
      )
      # #Not working with rowCallback
      # %>% DT::formatStyle("conditionID", 'pvalue',
      #                       color = DT::styleInterval(cuts = c(0.05), values = c('tomato', 'black'))
      # ) %>% DT::formatStyle("conditionID", 'padj',
      #                       fontWeight = DT::styleInterval(cuts = c(0.05), values = c('bold', 'normal'))
      # ) %>% DT::formatStyle(
      #   "conditionID",  target = 'row', backgroundColor = DT::styleEqual(levels = condition, values = "yellow")
      # )
    })

##GLOBAL HEATMAP
#Profile
    output$Heatmap_LoessGEvsAge <- renderHighchart({
      gene <- input$gene
 
      # Images not working
      # title <- paste(gene, " ",
      #                "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'> <img src='NCBI.png' title='NCBI' height='30px' style='padding-bottom:5px;'/></a>",
      #                "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'> <img src='geneCards.png' title='GeneCards' height='30px' style='padding-bottom:5px;'/></a>",
      #                sep = "")
      
      title <- paste(gene, " ",
                     "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>NCBI</a>",
                     " | ",
                     "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'>GeneCards</a>",
                     sep = "")
      
      
      
      #The gene list is loaded faster that the tissue list. In the meantime, an error message appears since no tissue is selcted
      validate(need(input$tissue != "", "")) #avoid the error message
      g <- Heatmap_LoessGEvsAge(geneData())
      g %>%
        hc_title(text = title,
                     useHTML = T) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Profile_", gene, "_acrossTissues"),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

#Alterations
    output$Heatmap_signficanceAlterationsvsAge <- renderHighchart({
       
      gene <- input$gene
      variable <- input$variable_Gene_Alteration
      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex"))
      # title <- paste(gene, " ",
      #                "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'> <img src='NCBI.png' title='NCBI' height='30px' style='padding-bottom:5px;'/></a>",
      #                "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'> <img src='geneCards.png' title='GeneCards' height='30px' style='padding-bottom:5px;'/></a>",
      #                sep = "")
      title <- paste(gene, " ",
                     "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>NCBI</a>",
                     " | ",
                     "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'>GeneCards</a>",
                     sep = "")
      validate(need(input$tissue != "", "")) #avoid the error message
      g <- Heatmap_signficanceAlterationsvsAge(p_Alterations())
      g %>%
        hc_title(text = title,
                     useHTML = T) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Alteration_", gene, "_across_", variable, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

##GENE-SPECIFIC
#PROFILE
#Scatter plot
    output$scatterGEvsAge <- renderHighchart({
      gene <- as.character(input$gene)
      tissue <- as.character(input$tissue)
      colored <- input$coloredBy
      geneInfo <- geneList$info[match(gene, geneList$gene)]
      condition <- ifelse(identical(as.character(input$shapedBy) == "Condition" & !is.na(input$condition), logical(0)) | as.character(input$shapedBy) == "none", "none", as.character(input$condition))
      #when a new gene is selected, geneData is updated sooner than input$tissue generating an error message
      validate(need(nrow(geneData()[geneData()$tissue == input$tissue, ]) != 0, "")) #allow to hide the error message
      #Title woth logo linking to genecards and NCBI websites
      # title <- paste(gene, " ",
      #                "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'> <img src='NCBI.png' title='NCBI' height='30px' style='padding-bottom:5px;'/></a>",
      #                "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'> <img src='geneCards.png' title='GeneCards' height='30px' style='padding-bottom:5px;'/></a>",
      #                sep = "")
      title <- paste(gene, " ",
                     "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>NCBI</a>",
                     " | ",
                     "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'>GeneCards</a>",
                     sep = "")
      g <- scatter_GEvsAge(geneData(), tissue, colored = colored, shaped = condition,
                           donorCondition = donorCondition, gene = gene, geneInfo = geneInfo)
      g %>%
        hc_title(text = title,
                 useHTML = T) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Profile_", gene, "_", tissue, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))

    })


#Data table gene expression and donor information
    output$geneExpression <- DT::renderDataTable({
      gene <- as.character(input$gene)
      gene <- gsub("-", "_", gene) # '-' creates issues when calling database so '-' were replaced by "_" in the database's names
      gene <- gsub("\\.", "_", gene) # '.' creates issues when calling database so '-' were replaced by "_" in the database's names
      tissue <- as.character(input$tissue)
      ymin <- ifelse(input$Selectedymin == "reset", input$Selectedymin, as.numeric(input$Selectedymin))
      ymax <- ifelse(input$Selectedymax == "reset", input$Selectedymax, as.numeric(input$Selectedymax))
      xmin <- ifelse(input$Selectedxmin == "reset", input$Selectedxmin, as.numeric(input$Selectedxmin))
      xmax <- ifelse(input$Selectedxmax == "reset", input$Selectedxmax, as.numeric(input$Selectedxmax))
      selectedPoint <- ifelse(is.null(input$ClickedStatus), "select", input$ClickedStatus)
      #selectedPoint <- "select" #remove possiblity of choosing a single point
      #To add this option
      #1. uncomment the previous selectedPoint
      #2. assign the value TRUE to allowPointSelect in the scatter_GEvsAge function
      #3. uncomment the events = list(click = ClickFunction) in the scatter_GEvsAge function
      indexPoint <- input$ClickedIndex + 1 #JS starts counting at 0

      condition <- ifelse(identical(as.character(input$shapedBy) == "Condition" & !is.na(input$condition), logical(0)) | as.character(input$shapedBy) == "none", "none", as.character(input$condition))

      #when a new gene is selected, geneData is updated sooner than input$tissue generating an error message
      validate(need(nrow(geneData()[geneData()$tissue == input$tissue, ]) != 0, "")) #allow to hide the error message

      if (condition == "none")
      {
        datatable_donorInfo(geneData(), tissue, selectedPoint, xmin, xmax, ymin, ymax, indexPoint, donorCondition = donorCondition, conditionID)
      } else
      {
        datatable_donorInfo(geneData(), tissue, selectedPoint, xmin, xmax, ymin, ymax, indexPoint, donorCondition = donorCondition, conditionID) %>%
          DT::formatStyle(condition, backgroundColor = "lightsalmon")
      }



    })

#Table downloading
    #Preparation of the table to download
    output$geneExpressionDownload <- DT::renderDataTable({
      tissue <- as.character(input$tissue)
      gene <- input$gene
      tmp <- geneData()
      tmp <- tmp[tmp$tissue == tissue,]
      tmp$sex[tmp$sex == 1] <- "Male"
      tmp$sex[tmp$sex == 2] <- "Female"
      #Change age format to give an interval and not the exact age
      tmp$age <- as.character(cut(tmp$age, breaks = seq(20,70, by = 5), include.lowest = T))

      #Add and rename condition
      temp <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), -1]
      temp[is.na(temp)] <- "Unknown"
      temp[temp == 96] <- "Unknown"
      temp[temp == 97] <- "Unknown"
      temp[temp == 99] <- "Unknown"
      temp[temp == 1] <- "Positive"
      temp[temp == 0] <- "Negative"

      tmp <- tmp[, c("samp_id", "age", "expression", "sex")]
      colnames(tmp) <- c("GTEx sample ID", "Age (years)", "Gene Expression (log(CPM))", "Sex of the donor")
      tmp <- cbind(tmp, temp)
      rownames(tmp) <- 1:nrow(tmp)
      filename <- paste("voyAGEr_GEdata_", gene, "_", tissue, sep = '')
      DT::datatable(tmp,
                    escape = F,
                    options = list(dom = "B",
                                   pageLength = nrow(tmp),
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename)),
                                   columnDefs = list(list(visible=FALSE, targets = 0:ncol(tmp))),
                                   headerCallback = JS(
                                     "function( thead, data, start, end, display ) {
            $(thead).closest('thead').find('th').each(function(){
              $(this).css('color', 'red').css('border','none');
            });
            }"
                                   ),
                                   initComplete = JS(
                                     "function(settings) {
            var table = settings.oInstance.api();
            $(table.table().node()).removeClass('no-footer');
            }")),
                    extensions = "Buttons"
      )


    })


#ALTERATIONS
    output$Line_signficanceAlterationsvsAge <- renderHighchart({
      #when a new gene is selected, geneData is updated sooner than input$tissue generating an error message
      validate(need(nrow(geneData()[geneData()$tissue == input$tissue, ]) != 0, "")) #allow to hide the error message

      gene <- input$gene
      variable <- input$variable_Gene_Alteration
      tissue <- input$tissue
      validate(need(!is.null(tissue) & tissue != "" & tissue != "All tissues", ""))
      pvalueData <- p_Alteration_gene()
      # title <- paste(gene, " ",
      #                "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'> <img src='NCBI.png' title='NCBI' height='30px' style='padding-bottom:5px;'/></a>",
      #                "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'> <img src='geneCards.png' title='GeneCards' height='30px' style='padding-bottom:5px;'/></a>",
      #                sep = "")
      title <- paste(gene, " ",
                     "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>NCBI</a>",
                     " | ",
                     "<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'>GeneCards</a>",
                     sep = "")

      validate(need(!is.null(pvalueData), "No possible analysis for this alteration in this tissue"))

      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename
      g <- Line_signficanceAlterationsvsAge(gene, pvalueData, geneData = geneData(), tissue, coloredBy =  variable) %>%
        hc_title(text = title, useHTML = T) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Alteration_", gene, "_across_", variable, "_", tissue, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

# TISSUE SECTION ----------------------------------------------------------

##GLOBAL HEATMAP
    output$AgeWavesOrganism <- renderHighchart({
      variable <- input$variable_model_2
      updatePickerInput(session, inputId = "tissue_2",
                        choices = c("All tissues", rownames(AgeWavesList[[variable]])[order(rownames(AgeWavesList[[variable]]))]),
                        selected = "All tissues")
      g <- Heatmap_log10FDRvsAge(variable = variable, AgeWaves = AgeWavesList)

      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename
      g %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Alteration_across_", variable, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

##TISSUE_SPECIFIC PEAKS
    #Histogram distribution samples for chosen tissue
    output$Histogram_DistributionSample <- renderHighchart({
      tissue <- input$tissue_2
      # #SQL database
      # gene <- dbGetQuery(DBconnection, paste("SELECT * FROM", "PGD", sep = " ")) #gene PGD chosen since expressed in all tisssue
      # #RDS files
      # gene <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_GeneExpressionV2/",
      #                       "PGD", ".RDS", sep = ""))
      #RDS files online
      gene <- readRDS(paste("data/DB_GeneExpressionV2/",
                            "PGD", ".RDS", sep = ""))

      validate(need(tissue != "" & tissue != "All tissues", ""))
      Histogram_DistributionSample(tissue, gene)
    })


    #% significance gene altered  = f(age)
    output$SigvsAge <- renderHighchart({
      variable <- input$variable_model_2
      tissue <- input$tissue_2
      validate(need(tissue != "All tissues" & tissue != "", ""))
      validate(need(nrow(peak[peak$tissue == tissue & peak$variable == variable,]) != 0, "No analysis for this tissue and variable"))
      g <- Line_sigGenesvsAge(tissue, variable, peak)

      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename
      g %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_%genesAltered_across_", variable, "_", tissue, sep = ''),
                     buttons = list(contextButton = list(align = "right",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

    #Data table of the gene differential expression signgificance
    output$genePeakTable <- DT::renderDataTable({
      
      validate(need(!is.null(input$peakClicked), "")) # A peak must be chosen
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      peakClicked <- input$peakClicked
      validate(need(!is.na(match(peakClicked, colnames(p()))), "")) # Matching column between age and peak must be found (avoid error when change tissue)
      validate(need(input$peakClickedVariable == variable, "")) #When a peak is clicked and the variable is changed, enable to remove the table while no peak is selected from the new polot
      tmp <- data.frame(gene = as.character(rownames(p())), 
                        pvalue = round(p()[, match(peakClicked, colnames(p()))], 6))
      tmp <- tmp[order(tmp$pvalue),]
      rownames(tmp) <- 1:nrow(tmp)

      tmp$Info <- paste(paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", tmp$gene,"' target='_blank'>", "GeneCards", "</a>"),
                             paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", tmp$gene,"' target='_blank'>", "NCBI", "</a>"), sep = ", ")
      filename <- paste("DEG_across_", variable, "_", tissue, "_", peakClicked, "yo", sep = "")

      DT::datatable(tmp,
                    escape = F,
                    caption = "Differential expression significance",
                    selection = "single", #only one row selected at a time
                    #filter = "top", #can filter the results per column
                    options = list(deferRender = TRUE, scrollY = 600, scroller = TRUE, scollX = T,
                                   dom = "Bfti",
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename))),
                    colnames = c("Gene", "p-value", "Info"),
                    rownames = T,
                    extensions = c("Scroller", "Buttons")
      )
    })

    #Name of the selected gene based on the chosen peak
    selectedGene <- reactive({
      peakClicked <- input$peakClicked
      indexRowClicked <- input$genePeakTable_rows_selected
      validate(need(!is.na(match(peakClicked, colnames(p()))), ""))
      tmp <- data.frame(gene = as.character(rownames(p())),
                        pvalue = round(p()[, match(peakClicked, colnames(p()))], 6))
      tmp <- tmp[order(tmp$pvalue),]
      gene <- as.character(tmp$gene[indexRowClicked])
      gene <- ifelse(!is.character(gene), NULL, gene)
    })

    #Import of the gene information for the chosen gene in Tissue tab
    geneData_Tissue <- reactive({
      gene <- selectedGene()
      gene <- gsub("-", "_", gene) # '-' creates issues when calling database so '-' were replaced by "_" in the database's names
      gene <- gsub("\\.", "_", gene) # '.' creates issues when calling database so '.' were replaced by "_" in the database's names
      # #SQL database
      # a <- dbGetQuery(DBconnection, paste("SELECT * FROM", gene, sep = " "))
      #RDS files local
      # a <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_GeneExpressionV2/",
      #                    gene, ".RDS", sep = ""))
      #RDS file onlibe
      a <- readRDS(paste("data/DB_GeneExpressionV2/",
                         gene, ".RDS", sep = ""))
      a
    })

    #Plot of the Lm pvalues over age
    output$Line_pvaluevsAge <- renderHighchart({
      validate(need(nrow(peak[peak$tissue == input$tissue_2 & peak$variable == input$variable_model_2,]) != 0, ""))
      validate(need((input$peakClickedState == 'hover' || input$peakClickedState == '') & input$variable_model_2 == input$peakClickedVariable, "Select a peak"))
      validate(need(!is.na(selectedGene()), "Select a gene"))
      validate(need(input$peakClickedVariable == input$variable_model_2, "")) #When a peak is clicked and the variable is changed, enable to remove the plot while no peak is selected from the new polot
      tissue <- input$tissue_2
      gene <- selectedGene() 
      pvalueData <- p()
      coloredBy <- input$variable_model_2
      geneInfo <- geneList$info[match(gene, geneList$gene)]

      variable <- input$variable_model_2
      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename

      g <- Line_pvaluevsAge(gene = gene, pvalueData = pvalueData, geneData = geneData_Tissue(), tissue = tissue, coloredBy = coloredBy) %>%
        hc_title(text = gene) %>%
        hc_subtitle(text = geneInfo)
      g %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Alteration_across_", variable, "_", gene, "_", tissue, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

##TISSUE_SPECIFIC ENRICHMENT


    output$pathwaySpecific <- DT::renderDataTable({
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQL database
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      #RDS file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files online
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))
      pathway <- unique(enrich$pathway[enrich$variable == variable])
      pathway <- pathway[order(as.character(pathway))]

      DT::datatable(as.data.frame(gsub("REACTOME_", "", pathway)),
                    escape = F,
                    colnames = c(""),
                    options = list(deferRender = TRUE, scrollY = 400, scroller = TRUE, scollX = T,
                                   rowCallback = JS(
                                     "function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {",
                                     "var full_text = aData[0]",
                                     "$('td:eq(0)', nRow).attr('title', full_text);",
                                     "}")),
                    selection = list(mode = "multiple",
                                     selected = 1),
                    rownames = F,
                    extensions = "Scroller"
      )

    })

    pathwaySelected <- reactive({
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQl database
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      # #RDS file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files onlibe
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))
      pathway <- unique(enrich$pathway[enrich$variable == variable])
      pathway <- pathway[order(as.character(pathway))]
      pathway[input$pathwaySpecific_rows_selected]
    })

    #Pathway-specific NES = f(age)
    output$EnrichmentAllOrSpecific <- renderUI({
      tagList(
        awesomeRadio(inputId = "EnrichmentAllOrSpecific", label = "Pathway", inline = T,
                     choices = c("All pathways", "Select:"), checkbox = T),
        conditionalPanel(condition = "input.EnrichmentAllOrSpecific == 'Select:'",
                         div(DT::dataTableOutput("pathwaySpecific"), style = "font-size: 70%; width: 100%",server = F)
        )
      )
    })

    #Heatmap fo the NES for a given tissue adn variable
    output$Heatmap_NESAgevsPathway <- renderCombineWidgets({
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      pathway <- input$EnrichmentAllOrSpecific
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQL databasle
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      # #RDs file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files online
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))
      SigPeak <- peak[peak$tissue == tissue & peak$variable == variable, ]

      validate(need(!is.null(pathway), ""))
      #Which patwhay tp show
      if (pathway == "Select:")
      {
        validate(need(length(pathwaySelected())!=0, "Choose pathway(s) of interest"))
      }
      if (pathway == "All pathways")
      {
        pathway <- TRUE
      } else
      {
        enrich <- enrich[enrich$pathway %in% pathwaySelected(), ]
        pathway <- FALSE
      }

      validate(need(length(which(SigPeak$sig == "sig")) != 0, "No significant peak or analysis done for this tissue and variable"))
      if (pathway == TRUE)
      {
        validate(need(length(enrich$pathway[which(enrich$padj <= 0.05 & enrich$age %in% SigPeak$age[SigPeak$sig == "sig"])]) != 0 , "No significant enrichment"))
      }

      g <- Heatmap_NESAgevsPathway(peak = peak, tissue = tissue, variable = variable, enrichment = enrich, Cluster_Reactome = Cluster_Reactome, ALL = pathway)

      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename
      g
    })

    #family presentation
    output$Heatmap_family <- renderHighchart({
      #To be linked with the enrichment heatmap, the family presentation only occurs when the heatmap is shown
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQL database
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      # #RDS file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files online
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))

      SigPeak <- peak[peak$tissue == tissue & peak$variable == variable, ]
      validate(need(length(which(SigPeak$sig == "sig")) != 0, ""))
      validate(need(length(enrich$pathway[which(enrich$padj <= 0.05 & enrich$age %in% SigPeak$age[SigPeak$sig == "sig"])]) != 0, ""))
      Heatmap_family(Cluster_Reactome)
    })
    output$wordcloud_family <- renderHighchart({
      #To be linked with the enrichment heatmap, the family presentation only occurs when the heatmap is shown
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQL database
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      #RDS file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files online
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))

      SigPeak <- peak[peak$tissue == tissue & peak$variable == variable, ]
      validate(need(length(which(SigPeak$sig == "sig")) != 0, ""))
      validate(need(length(enrich$pathway[which(enrich$padj <= 0.05 & enrich$age %in% SigPeak$age[SigPeak$sig == "sig"])]) != 0, ""))
      validate(need(!is.null(input$familyClicked), "Select a family"))
      wordcloud_family(Cluster_Reactome_affiliation, family = input$familyClicked) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_WordCloud_family-", input$familyClicked, sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))

    })

    output$listPathway_family <- DT::renderDataTable({
      #To be linked with the enrichment heatmap, the family presentation only occurs when the heatmap is shown
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      #Call to the databse
      modifiedTissuename <- gsub("-", "_", tissue)
      modifiedTissuename <- gsub(" ", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\(", "", modifiedTissuename)
      modifiedTissuename <- gsub("\\)", "", modifiedTissuename)
      # #SQL database
      # enrich <- dbGetQuery(DBconnection2, paste("SELECT * FROM", modifiedTissuename, sep = " "))
      #RDS file local
      # enrich <- readRDS(paste("/genedata/home/arthur/Documents/GTEx/shiftingLM_genes/shiny_app/RDS files/DB_PeakEnrichment.db/",
      #                         modifiedTissuename, ".RDS", sep = ""))
      #RDS files online
      enrich <- readRDS(paste("data/DB_PeakEnrichment.db/",
                              modifiedTissuename, ".RDS", sep = ""))

      SigPeak <- peak[peak$tissue == tissue & peak$variable == variable, ]
      validate(need(length(which(SigPeak$sig == "sig")) != 0, ""))
      validate(need(length(enrich$pathway[which(enrich$padj <= 0.05 & enrich$age %in% SigPeak$age[SigPeak$sig == "sig"])]) != 0, ""))
      validate(need(!is.null(input$familyClicked), "Select a family"))


      tmp <- rbind(data.frame(pathway = Cluster_Reactome_affiliation[[input$familyClicked]][!is.na(Cluster_Reactome_affiliation[[input$familyClicked]][,1]),1]),
                   data.frame(pathway = Cluster_Reactome_affiliation[[input$familyClicked]][!is.na(Cluster_Reactome_affiliation[[input$familyClicked]][,2]),2]))

      filename <- paste("voyAGEr_WordCloud_family-", input$familyClicked, sep = '')
      DT::datatable(tmp,
                    escape = F,
                    selection = list(mode = "none"),
                    colnames = paste("Family", input$familyClicked),
                    rownames = F,
                    options = list(deferRender = TRUE, scrollY = 270, scroller = TRUE, scollX = T,
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename))), #to remove the searching and show tab
                    extensions = c("Scroller", "Buttons")
      )
    })




    #Manual enrichment based on user-provided geneset
    #List of user-provided gene
    manualEnrichment <- eventReactive(input$manualEnrichment,{
      #Conversion into character vector
      toupper(as.character(unlist(strsplit(input$geneListEnrichment,"\n"))))
    })

    manualEnrichmentThreshold <- eventReactive(input$manualEnrichment, {
      as.numeric(input$manualEnrichmentThreshold)
    })

    output$Geneset <- renderText({
      intersect(manualEnrichment(), rownames(p()))
      })

    output$geneList <- renderUI({
      column(width = 12,
             HTML("<b>", length(intersect(manualEnrichment(), rownames(p()))), "genes considered in the analysis</b>"),
             verbatimTextOutput("Geneset", placeholder  =TRUE)
      )
    })

    #Contingency Table
    contingencyManualEnrichment <- reactive({
      threshold <- manualEnrichmentThreshold()
      ChosenAge <- input$ageContingencyTableClicked
      pColumn <- match(as.numeric(ChosenAge), as.numeric(colnames(p())))

      tmp <- as.data.frame(
        matrix(c(length(which(p()[match(manualEnrichment(), rownames(p())), pColumn] <= threshold)), #gene from geneset DEG
                 length(which(p()[match(manualEnrichment(), rownames(p())), pColumn] > threshold)), #genes from geneset not DEG
                 length(which(p()[,pColumn] <= threshold)), #genes from all DEG
                 length(which(p()[,pColumn] > threshold))), #genes from all not DEG
               nrow = 2))
      tmp
    })

    observeEvent(input$ageContingencyTableClicked, {
      showModal(modalDialog(
        title = "Contingency Table",
        column(width = 12,
               HTML("The table summarises the numbers of studied genes differentially expressed and not that are part of the user-defined gene set and that are not.
                    <br> It is the basis for the Fisher’s exact test for the significance of the user-specified geneset’s enrichment in differentially expressed genes."),
               DT::renderDataTable(contingencyManualEnrichment(),
                                   selection = list(mode = "none"),
                                   options = list(dom = "t"),
                                   rownames = c("Differentially Expressed", "Not differentially expressed"),
                                   colnames = c("User's genes", "Other genes"))
        ),
        easyClose = TRUE,
        fade = TRUE,
        footer = modalButton(label = "", icon = icon("times")),
        style = 'padding:4px'
      ))
    })



    output$Line_FishertestEnrichmentManualvsAge <- renderHighchart({
      gene <- manualEnrichment()
      threshold <- manualEnrichmentThreshold()
      tissue <- input$tissue_2
      variable <- input$variable_model_2
      variable <- ifelse(variable == "Age", "Age", ifelse(variable == "Gender", "Sex", "Age&Sex")) #change name for download filename
      validate(need(!is.null(gene), "Fill geneset"))
      Line_FishertestEnrichmentManualvsAge(p(), gene, threshold = threshold) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_genesetEnrichment_", tissue, "_across_", variable , sep = ''),
                     buttons = list(contextButton = list(align = "left",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })



# MODULE ------------------------------------------------------------------

##SIDEBAR
    observe({
      tissue <- input$tissue_module
      moduleList <- names(moduleInfo$moduleGeneList[[tissue]])

      updatePickerInput(session, "selectedModule",
                        choices =  c("All modules", as.character(moduleList)[order(as.character(moduleList))]),
                        selected = "All modules"
                        #choicesOpt = list(style = c("black", paste("color:", gsub("ME", "", moduleList)))) #module name colored accoridng to their name
                        )

      updateAwesomeRadio(session, "moduleColoredBy", choices = list("All" = "none", "Sex" = "Sex", "Technicality" = "Technicality"), selected = "none",
                         inline = T, checkbox = TRUE)
      updateAwesomeRadio(session, "moduleShapedBy", choices = list("All" = "none", "Condition*" = "Condition"), selected = "none",
                         inline = T, checkbox = TRUE)

    })

    output$moduleColoredByTechnicality <- renderUI({
      pickerInput("moduleTechnicality", "Select:",
                  choices = colnames(technicalCondition)[-1][order(colnames(technicalCondition)[-1])],
                  options = list(`live-search` = TRUE))
    })

    output$moduleShapedByCondition <- renderUI({
      pickerInput("moduleCondition", "Select:",
                  choices = colnames(donorCondition)[-1][order(colnames(donorCondition)[-1])],
                  options = list(`live-search` = TRUE))
    })


    #Outside bar
    output$modulesOutsideBar <- renderUI({
      if (input$selectedModule != "All modules" & input$selectedModule != "" & length(moduleInfo$moduleGeneList[[input$tissue_module]]) != 0 & !is.null(length(moduleInfo$moduleGeneList[[input$tissue_module]]) == 0))
      {
        DT::dataTableOutput(outputId = "moduleGeneList")
      } else if (input$modules == "Cell types")
      {
        column(width = 12,
               DT::dataTableOutput(outputId = "referenceCellType"),
               HTML("<small>Id refers to the first author(s) of the paper or the name given to the dataset from where the cell markers are retrieved.</small>")
        )
      } else if (input$modules == "Pathways")
      {
        column(width = 12,
               pickerInput("family_modules", label = "Family", choices = 1:length(Cluster_Reactome_affiliation),
                           options = list(`live-search` = TRUE,
                                          title = "Select a family")),
               highchartOutput("wordcloud_family_enrichment")
        )
      } else if (input$modules == "Diseases" & input$diseaseMethod == "DOSE")
      {
        disease <- as.character(unique(moduleInfo$diseaseEnrichment[[input$tissue_module]]$dataframe$Description))
        column(width = 12,
               multiInput(inputId = "selectedDiseases", label = "Diseases",  choices = disease,
                          selected = c("Age related macular degeneration", "Cardiomyopathies", "Connective tissue diseases",
                                       "Diffuse scleroderma", "Musculoskeletal diseases", "Pneumonia")),
               HTML("<small>Be aware that the more diseases are displayed the less readable the heatmap may become.</small>")
        )
      } else if (input$modules == "Diseases" & input$diseaseMethod == "Manual")
      {
        disease <- as.character(unique(moduleInfo$diseaseEnrichmentManual[[input$tissue_module]]$dataframe$disease))
        column(width = 12,
               multiInput(inputId = "selectedDiseasesManual", label = "Diseases",  choices = disease,
                          selected = disease[1:3]),
               HTML("<small>Be aware that the more diseases are displayed the less readable the heatmap will be.</small>")
        )
      }
    })


##ALL MODULES

    #Loess expression  = f(age)
    output$Heatmap_LoessMEvsAge <- manipulateWidget::renderCombineWidgets({
      tissue <- input$tissue_module
      moduleExpression <- moduleInfo$moduleExpression[moduleInfo$moduleExpression$tissue == tissue, ]
      Heatmap_LoessMEvsAge(moduleExpression)
    })

    #Cell type composition
    output$referenceCellType <- DT::renderDataTable({
      tmp <- moduleInfo$sources[moduleInfo$sources$tissue == input$tissue_module,]

      DT::datatable(tmp[, c("reference", "article")],
                    caption = "Cell types markers sources",
                    escape = F,
                    options = list(dom = "t", ordering = F),
                    selection = list(mode = "single",
                                     selected = 1), #only one row selected at a time
                    #filter = "top", #can filter the results per column
                    colnames = c("Id", "Ref."),
                    rownames = F
      )
    })

    selectedReference <- reactive({
      indexRowSelected <- input$referenceCellType_rows_selected
      tmp <- moduleInfo$sources[moduleInfo$sources$tissue == input$tissue_module,]
      as.character(tmp$reference[indexRowSelected])
    })

    output$Heatmap_FisherTest_cellComposition_All <- renderHighchart({
      tissue <- input$tissue_module
      moduleCellType <- moduleInfo$fisherTest_cellComposition[moduleInfo$fisherTest_cellComposition$tissue == tissue, ]
      moduleCellType <- moduleCellType[moduleCellType$reference == selectedReference(),]
      validate(need(length(selectedReference()) == 1, "Select a reference for the cell type markers"))
      Heatmap_FisherTest_cellComposition(moduleCellType) %>%
        hc_title(text = paste0("Cell types markers from ", selectedReference())) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_module_", tissue, "_cellType", sep = ''),
                     buttons = list(contextButton = list(align = "right",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })


    #Functional Enrichment
    output$Heatmap_moduleEnrichment_All <- renderHighchart({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$moduleEnrichment[[tissue]]$matrix
      Heatmap_moduleEnrichment(enrich = moduleEnrichment, Cluster_Reactome, tissue = tissue)

    })

    #Family presentation
    output$wordcloud_family_enrichment <- renderHighchart({

      validate(need(input$family_modules %in% 1:length(Cluster_Reactome_affiliation), "Select a family"))
      wordcloud_family(Cluster_Reactome_affiliation, family = as.numeric(input$family_modules))

    })


    #Disease Enrichment - DOSE
    output$Heatmap_moduleDiseaseEnrichment_All <- renderHighchart({
      tissue <- input$tissue_module
      moduleDiseaseEnrichment <- moduleInfo$diseaseEnrichment[[tissue]]$matrix
      validate(need(length(input$selectedDiseases) >= 1, "Select disease(s) of interest"))
      Heatmap_diseaseEnrichment(diseaseEnrichment = moduleDiseaseEnrichment, disease = input$selectedDiseases) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_module_", tissue, "_diseaseDOSE", sep = ''),
                     buttons = list(contextButton = list(align = "right",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))

    })

    #Disease Enrichment - Manual
    output$Heatmap_moduleDiseaseEnrichmentManual_All <- renderHighchart({
      tissue <- input$tissue_module
      moduleDiseaseEnrichment <- moduleInfo$diseaseEnrichmentManual[[tissue]]$matrix
      validate(need(length(input$selectedDiseasesManual) >= 1, "Select disease(s) of interest"))
      validate(need(!is.na(sum(match(input$selectedDiseasesManual, moduleDiseaseEnrichment$Disease))), ""))
      Heatmap_diseaseEnrichment(diseaseEnrichment = moduleDiseaseEnrichment, disease = input$selectedDiseasesManual) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_module_", tissue, "_diseaseManual", sep = ''),
                     buttons = list(contextButton = list(align = "right",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))

    })

##MODULE_SPECIFIC

    #Eigengene expression over time
    output$Scatter_MEvsAge <- renderHighchart({
      tissue <- input$tissue_module
      module <- input$selectedModule
      coloredBy <- input$moduleColoredBy
      if (coloredBy == "Technicality")
      {
        coloredBy <- as.character(input$moduleTechnicality)
        validate(need(!is.null(input$moduleTechnicality), ""))
      }
      condition <- ifelse(identical(as.character(input$moduleShapedBy) == "Condition" & !is.na(input$moduleCondition), logical(0)) | as.character(input$moduleShapedBy) == "none", "none", as.character(input$moduleCondition))
      moduleExpression <- moduleInfo$moduleExpression[moduleInfo$moduleExpression$tissue == tissue, ]

      #allow to hide the error message
      validate(need(nrow(moduleExpression[moduleExpression$module == module, ]) != 0, ""))
      Scatter_MEvsAge(MEexpressionData = moduleExpression, module = module, colored = coloredBy, shaped = condition,
                      donorCondition = donorCondition, technicalCondition = technicalCondition) %>%
        hc_exporting(enabled = T,
                     filename = paste("voyAGEr_Profile_", module, "_", tissue, sep = ''),
                     buttons = list(contextButton = list(align = "right",
                                                         symbolStrokeWidth = 4,
                                                         symbolSize = 20,
                                                         symbolStroke = '#4D4D4D',
                                                         symbol = "menu")))
    })

    #Cell type
    output$Heatmap_FisherTest_cellComposition_One <- manipulateWidget::renderCombineWidgets({
      tissue <- input$tissue_module
      moduleCellType <- moduleInfo$fisherTest_cellComposition[moduleInfo$fisherTest_cellComposition$tissue == tissue, ]
      module <- input$selectedModule
validate(need(nrow(moduleCellType[moduleCellType$module == module,]) !=0 , ""))
      Treemap_moduleFisherTest(pvalueFisherTest = moduleCellType, module = module)

    })

    output$referenceMarker <- renderUI({
      HTML(paste("Rectangles’ areas and color brightness are respectively proportional to the significance (-log10(p)) and the enrichment coefficient (odds ratio) of the module’s enrichment in cell type markers from the literature computed with Fisher’s exact tests.
                 <br><i><small> Sources: ",
                 paste(moduleInfo$sources$reference[moduleInfo$sources$tissue == input$tissue_module], " (",
                       moduleInfo$sources$article[moduleInfo$sources$tissue == input$tissue_module], "), ", collapse = "", sep = ""), "</small></i>"))
    })


    #Functional Enrichment
    output$DT_moduleEnrichment_One <- DT::renderDataTable({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$moduleEnrichment[[tissue]]$dataframe
      module <- input$selectedModule
      moduleEnrichment <- moduleEnrichment[moduleEnrichment$module == module,]
      moduleEnrichment <- moduleEnrichment[order(moduleEnrichment$pvalue),]
      moduleEnrichment$padj <- round(moduleEnrichment$padj, 5)
      moduleEnrichment$pvalue <- round(moduleEnrichment$pvalue, 5)

      filename <- paste("voyAGEr_Enrichment_", module, "_", tissue, sep = '')
      DT::datatable(as.data.frame(cbind(moduleEnrichment[, c("pathway", "pvalue", "padj")], table = shinyInput(actionButton, nrow(moduleEnrichment),'button_', label = "View", onclick = 'Shiny.onInputChange(\"select_Reactomebutton\",  this.id)' ))),
                    escape = F,
                    selection = list(mode = "none"),
                    colnames = c("Reactome pathway", "p-value", "adjusted p-value", "Contingency table"),
                    rownames = F,
                    options = list(deferRender = TRUE, scrollY = 300,
                                   scroller = TRUE, scrollCollapse = TRUE, scrollX = TRUE,
                                   dom = "Bfti",
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename)),
                                   columnDefs = list(list(className = 'dt-center', targets = 1:3))), #to remove the searching and show tab
                    extensions = c("Scroller", "Buttons")
      ) %>% DT::formatStyle("pathway", 'pvalue',
                            color = DT::styleInterval(cuts = c(0.05), values = c('tomato', 'black'))
      ) %>% DT::formatStyle("pathway", 'padj',
                            fontWeight = DT::styleInterval(cuts = c(0.05), values = c('bold', 'normal')))

    })
    contingencyPathway <- reactive({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$moduleEnrichment[[tissue]]$dataframe
      module <- input$selectedModule
      moduleEnrichment <- moduleEnrichment[moduleEnrichment$module == module,]
      moduleEnrichment <- moduleEnrichment[order(moduleEnrichment$pvalue),]

      SelectedRow <- as.numeric(strsplit(as.character(input$select_Reactomebutton), "_")[[1]][2])
      #Disease selected for the table
      pathway <- moduleEnrichment$pathway[SelectedRow]


      contingencyPathway <- moduleInfo$contigencyTable$moduleEnrichment
      contingencyPathway <- contingencyPathway[contingencyPathway$tissue == tissue & contingencyPathway$module == module,]
      contingencyPathway <- contingencyPathway[match(pathway, contingencyPathway$pathway), ]
      contingencyPathway <- as.data.frame(matrix(c(contingencyPathway$INmoduleINpathway,
                                                   contingencyPathway$INmoduleOUTpathway,
                                                   contingencyPathway$OUTmoduleINpathway,
                                                   contingencyPathway$OUTmoduleOUTpathway), nrow = 2))
      contingencyPathway


    })

    #TO get a pop up with the contigency table
      observeEvent(input$select_Reactomebutton, ignoreInit = TRUE,{
        showModal(modalDialog(
          title = "Contingency Table",
          column(width = 12,
                 HTML("The table summarises the numbers of studied genes belonging to the pathway and not that are part of the module and that are not.
                      It is the basis for the Fisher’s exact test for the significance of the module’s enrichment in genes from the pathway."),
                 DT::renderDataTable(contingencyPathway(),
                                     selection = list(mode = "none"),
                                     options = list(dom = "t", ordering = F),
                                     colnames = c("IN module", "OUT module"),
                                     rownames = c("IN pathway", "OUT pathway"),
                                     class = "nowrap")
          ),
          easyClose = TRUE,
          fade = TRUE,
          footer = modalButton(label = "", icon = icon("times")),
          style = 'padding:4px'
        ))
      })


    #Disease enrichment - DOSE
    output$DT_diseaseEnrichment <- DT::renderDataTable({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$diseaseEnrichment[[tissue]]$dataframe
      module <- input$selectedModule
      moduleEnrichment <- moduleEnrichment[moduleEnrichment$module == module, ]
      moduleEnrichment <- moduleEnrichment[order(moduleEnrichment$pvalue),]
      moduleEnrichment$Description <- toupper(moduleEnrichment$Description)
      moduleEnrichment$pvalue <- round(moduleEnrichment$pvalue, 5)
      moduleEnrichment$p.adjust <- round(moduleEnrichment$p.adjust, 5)

      filename <- paste("voyAGEr_diseaseDOSE_", module, "_", tissue, sep = '')
      DT::datatable(moduleEnrichment[, c("Description", "pvalue", "p.adjust")],
                    escape = F,
                    selection = list(mode = "none"),
                    colnames = c("Disease", "p-value", "adjusted p-value"),
                    rownames = F,
                    options = list(deferRender = TRUE, scrollY = 300, scroller = TRUE, scollX = T,
                                   dom = "Bfti",
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename))), #to remove the searching and show tab
                    extensions = c("Scroller", "Buttons")
      ) %>% DT::formatStyle("Description", 'pvalue',
                            color = DT::styleInterval(cuts = c(0.05), values = c('tomato', 'black'))
      ) %>% DT::formatStyle("Description", 'p.adjust',
                            fontWeight = DT::styleInterval(cuts = c(0.05), values = c('bold', 'normal')))

    })



    #Disease enrichment - Manual
    output$DT_diseaseEnrichmentManual <- DT::renderDataTable({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$diseaseEnrichmentManual[[tissue]]$dataframe
      module <- input$selectedModule
      moduleEnrichment <- moduleEnrichment[moduleEnrichment$module == module, ]
      moduleEnrichment <- moduleEnrichment[order(moduleEnrichment$pvalue),]
      moduleEnrichment$disease <- toupper(moduleEnrichment$disease)
      moduleEnrichment$pvalue <- round(moduleEnrichment$pvalue, 5)
      moduleEnrichment$padj <- round(moduleEnrichment$padj, 5)

      filename <- paste("voyAGEr_diseaseManual_", module, "_", tissue, sep = '')
      DT::datatable(as.data.frame(cbind(moduleEnrichment[, c("disease", "diseaseClass", "pvalue", "padj")],
                                        table = shinyInput(actionButton, nrow(moduleEnrichment),'button_', label = "View", onclick = 'Shiny.onInputChange(\"select_Diseasebutton\",  this.id)' ))),
                    escape = F,
                    selection = list(mode = "none"),
                    colnames = c("Disease", "Disease class", "p-value", "adjusted p-value", "Contingency table"),
                    rownames = F,
                    options = list(deferRender = TRUE, scrollY = 300, scroller = TRUE, scollX = T,
                                   dom = "Bfti",
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename)),
                                   columnDefs = list(list(className = 'dt-center', targets = 2:4))), #to remove the searching and show tab
                    extensions = c("Scroller", "Buttons")
      ) %>% DT::formatStyle("disease", 'pvalue',
                            color = DT::styleInterval(cuts = c(0.05), values = c('tomato', 'black'))
      ) %>% DT::formatStyle("disease", 'padj',
                            fontWeight = DT::styleInterval(cuts = c(0.05), values = c('bold', 'normal')))

    })
    #TO get a pop up with the contigency table
    contingencyDisease <- reactive({
      tissue <- input$tissue_module
      moduleEnrichment <- moduleInfo$diseaseEnrichmentManual[[tissue]]$dataframe
      module <- input$selectedModule
      moduleEnrichment <- moduleEnrichment[moduleEnrichment$module == module, ]
      moduleEnrichment <- moduleEnrichment[order(moduleEnrichment$pvalue),]

      SelectedRow <- as.numeric(strsplit(as.character(input$select_Diseasebutton), "_")[[1]][2])
      #Disease selected for the table
      disease <- moduleEnrichment$disease[SelectedRow]


      contingencyDisease <- moduleInfo$contigencyTable$diseaseEnrichmentManual
      contingencyDisease <- contingencyDisease[contingencyDisease$tissue == tissue & contingencyDisease$module == module,]
      contingencyDisease <- contingencyDisease[match(disease, contingencyDisease$disease), ]
      contingencyDisease <- as.data.frame(matrix(c(contingencyDisease$INmoduleINpathway,
                                                   contingencyDisease$INmoduleOUTpathway,
                                                   contingencyDisease$OUTmoduleINpathway,
                                                   contingencyDisease$OUTmoduleOUTpathway), nrow = 2))
      contingencyDisease


      })

    observeEvent(input$select_Diseasebutton, ignoreInit = TRUE,{
      showModal(modalDialog(
        title = "Contingency Table",
        column(width = 12,
               HTML("The table summarises the numbers of studied genes associated with the disease and not that are part of the module and that are not.
                     It is the basis for the Fisher’s exact test for the significance of the module’s enrichment in genes associated with the disease."),
               DT::renderDataTable(contingencyDisease(),
                                   selection = list(mode = "none"),
                                   options = list(dom = "t"),
                                   rownames = c("IN disease", "OUT disease"),
                                   colnames = c("IN module", "OUT module"))
        ),
        easyClose = TRUE,
        fade = TRUE,
        footer = modalButton(label = "", icon = icon("times")),
        style = 'padding:4px'
      ))
    })


    #Data table module's genes
    output$moduleGeneList <- DT::renderDataTable({
      tmp <- moduleInfo$moduleGeneList[[input$tissue_module]]
      tmp <- data.frame(gene = as.character(tmp[[input$selectedModule]]))
      if (nrow(tmp)!=0) #to avoid error message when change of tissue
      {
        tmp$Info <- paste(paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=", tmp$gene,"' target='_blank'>", "GeneCards", "</a>"),
                          paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", tmp$gene,"' target='_blank'>", "NCBI", "</a>"), sep = ", ")
      }
      validate(need(nrow(tmp) != 0, ""))

      filename <- paste("voyAGEr_", input$selectedModule, "_", input$tissue_module, sep = "")
      DT::datatable(tmp,
                    escape = F,
                    caption = paste("Module of", nrow(tmp), "genes"),
                    selection = "none", #only one row selected at a time
                    #filter = "top", #can filter the results per column
                    options = list(deferRender = TRUE, scrollY = 300, scroller = TRUE, scollX = T,
                                   #dom = 'Btf',
                                   buttons = list(list(extend='copy',
                                                       filename = filename),
                                                  list(extend='csv',
                                                       filename = filename),
                                                  list(extend='excel',
                                                       filename= filename))), #to remove the searching and show tab
                    colnames = c("Gene", "Info"),
                    rownames = F,
                    extensions = c("Scroller", "Buttons")
      )
    })

    geneSearchedInModules <- reactive({
    	geneSearched <- as.character(toupper(input$searchGeneInModule))
    	if (geneSearched != "")
    	{
    		tmp <- moduleInfo$moduleGeneList[[input$tissue_module]]
    		res <- lapply(tmp, function(ch) which(ch == geneSearched)) #search the list of the genes from the modules for the gene of interest
    		if (length(as.numeric((which(sapply(res, function(x) length(x) > 0) != 0)))))
    		{
    			result <- paste(input$searchGeneInModule, "found in ", names(which(sapply(res, function(x) length(x) > 0) == TRUE)))

    		} else
    		{
    			result <- paste(input$searchGeneInModule, "found in 0 module")

    		}
    	} else
    	{
    		result <- ""
    	}
    	result
    })

    output$geneInModules <- renderText({
    	geneSearchedInModules()
    })

    
    
    
    
    
    
  }
)

