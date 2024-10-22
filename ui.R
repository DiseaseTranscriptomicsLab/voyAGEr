library(shiny)
library(shinythemes)
library(shinyWidgets)
library(highcharter)
library(shinycssloaders)
library(manipulateWidget)
library(fontawesome)

sidebarPanel2 <- function (..., out = NULL, width = 4) 
{
  div(class = paste0("col-sm-", width), 
      tags$form(class = "well", ...),
      out
  )
}


shinyUI(fluidPage(
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(5) {
                           float: right;
                           }
                           "))),
  navbarPage(title = "voyAGEr", theme = shinytheme("flatly"),
             tabPanel("", #title of the tab 
                      tags$style(HTML(".big_icon_test {margin-top:-5px; font-size: 25px; }")), #change size icon
                      icon = icon(name = "home", class = "big_icon_test"),
                      mainPanel(width = 12,
                                tabsetPanel(id = "home", type = "pills",
                                            tabPanel("Overview", 
                                                     icon = icon("comments"),
                                                     column(width = 12, style = "margin-bottom:50px",
                                                            HTML('<center><img src="MainFigure_2.png" height="550"></center>')
                                                     ),
                                                     tags$footer(tags$a(href = 'http://imm.medicina.ulisboa.pt/group/distrans/',
                                                                        img(src = 'imm_logo.png', title = "Lab Webpage", height = "45px"), target = "blank"),
                                                                 HTML(""),
                                                                 tags$a(href = 'https://www.genomept.pt/',
                                                                        img(src = 'genomept.png', title = "Consortium Webpage", height = "40px",
                                                                            style = "padding-right:5px;"),
                                                                        target = "blank"), 
                                                                 align = "center", 
                                                                 style = "position: fixed;
                                                          left: 0;
                                                          bottom: 0;
                                                          width: 100%;
                                                          background-color: #ecf0f1;
                                                          color: black;
                                                          text-align: center;"),
                                                     tags$head(
                                                       tags$style(HTML('.bottom-right-text {
                                                                          position: absolute;
                                                                          bottom: 10px;
                                                                          right: 10px;
                                                                          font-size: 16px;
                                                                        }'))
                                                     ),
                                                     HTML('<div class="bottom-right-text"><i class="fa fa-github"></i><a href="https://github.com/DiseaseTranscriptomicsLab/voyAGEr/tree/main" target="_blank"> Version 2.0.2 (October 2024) </a></div>')
                                            ),
                                            tabPanel("Methodology",
                                                     icon = icon("cogs"),
                                                     column(width = 12, style = "margin-bottom:50px",
                                                            HTML('<center><img src="methodology.png" height="400"></center>'),
                                                            HTML("<font size='4'><br>Our shifting age range pipeline for linear modelling (ShARP-LM) consists in carrying out tissue-specific differential gene expression analyses across samples in age windows spanning 16 years centred in consecutive years of age. 
                                                       <br>For each window, the gene expression is linearly modelled considering the age, sex and age-sex interaction effects.
                                                       <br>The derived statistics enable to identify the age periods when major gene expression changes occur due to which of the 3 effects and to evaluate their functional enrichment.</font>")
                                                     ),
                                                     tags$footer(tags$a(href = 'http://imm.medicina.ulisboa.pt/group/distrans/',
                                                                        img(src = 'imm_logo.png', title = "Lab Webpage", height = "45px"), target = "blank"), 
                                                                 HTML(""),
                                                                 tags$a(href = 'https://www.genomept.pt/',
                                                                        img(src = 'genomept.png', title = "Consortium Webpage", height = "40px",
                                                                            style = "padding-right:5px;"),
                                                                        target = "blank"), 
                                                                 align = "center", 
                                                                 style = "position: fixed;
                                                          left: 0;
                                                          bottom: 0;
                                                          width: 100%;
                                                          background-color: #ecf0f1;
                                                          color: black;
                                                          text-align: center;")
                                            ))
                      )
             ),
             tabPanel("Gene", 
                      icon = icon("dna"),
                      sidebarLayout(
                        sidebarPanel2(width = 3,
                                      pickerInput("gene", label = "Gene", 
                                                  choices = c("Loading ..." = ""), 
                                                  options = list(`live-search` = TRUE), selected = "CDKN2A"),
                                      pickerInput("tissue", label = "Tissue", 
                                                  choices = c("Loading ..." = ""), 
                                                  options = list(`live-search` = TRUE),
                                                  choicesOpt = list(content = "<div style='color: red; </div>")),
                                      conditionalPanel(condition = "input.Gene_profile == 'Alteration'",
                                                       awesomeRadio(inputId = "variable_Gene_Alteration", label = "Alterations associated with", inline = T,
                                                                    choices = list("Age" = "Age", "Sex" = "Gender", "Age&Sex" = "Int_Age_Gender"),
                                                                    checkbox = TRUE)),
                                      conditionalPanel(condition = "input.tissue != 'All tissues' && input.tissue != '' && input.Gene_profile == 'Profile'",
                                                       awesomeRadio(inputId = "coloredBy", label = "Coloured By", 
                                                                    choices = c("Loading ...")),
                                                       awesomeRadio(inputId = "shapedBy", label = "Shaped By", 
                                                                    choices = c("Loading ...")),
                                                       HTML("<small>*A curve is independently fitted for each of the positive and negative conditions if at least 10 samples are known for each of them. </small>")
                                      ),
                                      conditionalPanel(condition = "input.shapedBy == 'Condition' && input.tissue != 'All tissues'  && input.Gene_profile == 'Profile'", uiOutput("shapedByCondition")),
                                      out = conditionalPanel(condition = "input.shapedBy == 'Condition' && input.tissue != 'All tissues'  && input.Gene_profile == 'Profile'",
                                                             column(width = 12,
                                                                    fluidRow(
                                                                      column(width = 8, tags$h4("Condition profiles")),
                                                                      column(width = 4, dropdownButton(
                                                                        HTML("A statistical test (Kruskal-Wallis) was computed to examine overall difference in gene expression medians between positive and negative conditions. 
                                                                      The p-values are corrected with the Benjamini & Hochberg approach.
                                                                      <br> Even if the method accounts for the high imbalance between positive and negative samples, results must still be viewed with caution."),
                                                                        circle = TRUE,  size = "sm", inline = T,
                                                                        icon = icon("info"), width = "300px",
                                                                        tooltip = tooltipOptions(title = "Click for more information")
                                                                      ))
                                                                    ),
                                                                    HTML("<small>Hover over a condition ID to get its associated description.</small>"),
                                                                    div(DT::dataTableOutput("KruskalWallisCondition"), style = "font-size: 90%; width: 100%",server = F)))
                        ),
                        mainPanel(width = 9,
                                  tabsetPanel(id = "Gene_profile", type = "pills",
                                              tabPanel("Profile",
                                                       icon = icon("chart-line"),
                                                       conditionalPanel(condition = "input.tissue == 'All tissues' || input.tissue == '' ", 
                                                                        column(width = 12,
                                                                               tags$h4("Gene expression over age across tissues"),
                                                                               withSpinner(highchartOutput("Heatmap_LoessGEvsAge", height = "800px"), 
                                                                                           color = "#2C3E50", type = 5, size = 0.5))
                                                       ),
                                                       conditionalPanel(condition = "input.tissue != 'All tissues' && input.tissue != ''", 
                                                                        column(width = 12,
                                                                               tags$h4(""),
                                                                               HTML("For data privacy purposes, a light randomization was attributed to each scattered point along the age axis. The plotted regression line (black) was fitted using the original values."),
                                                                               withSpinner(highchartOutput("scatterGEvsAge", height = "650px"), 
                                                                                           color = "#2C3E50", type = 5, size = 0.5),
                                                                               #downloadButton('downloadDataGeneExpression', 'Table'),
                                                                               DT::dataTableOutput("geneExpressionDownload"),
                                                                               div(DT::dataTableOutput("geneExpression"), style = "font-size: 90%; width: 100%",server = F)) 
                                                       )
                                              ),
                                              tabPanel("Alteration",
                                                       icon = icon("crosshairs"),
                                                       conditionalPanel(condition = "input.tissue == 'All tissues' || input.tissue == '' ",
                                                                        column(width = 12,
                                                                               tags$h4("Gene expression alterations over age across tissues"),
                                                                               HTML("The heatmap color scale represents the significance of the gene expression alterations (log-scaled)."),
                                                                               withSpinner(highchartOutput("Heatmap_signficanceAlterationsvsAge", height = "800px"), 
                                                                                           color = "#2C3E50", type = 5, size = 0.5))
                                                       ),
                                                       conditionalPanel(condition = "input.tissue != 'All tissues' && input.tissue != ''", 
                                                                        column(width = 12,
                                                                               tags$h4("Gene expression alterations over age"),
                                                                               #Description text different for each variable
                                                                               conditionalPanel(condition = "input.variable_Gene_Alteration == 'Age'",
                                                                                                HTML("The gene expression progression over age (top panel) is seen in parallel with its expression alteration's significance (bottom panel).
                                                                                                     <br> For data privacy purposes, a light randomization was attributed to each scattered point along the age axis. The plotted regression line (black) was fitted using the original values.
                                                                                                     <br> The overall p-value, t-statistic and logFC/year (in orange) was computed using the ShARP-LM model over the entire age range.
                                                                                                     "),
                                                                               ),
                                                                               conditionalPanel(condition = "input.variable_Gene_Alteration == 'Gender'",
                                                                                                HTML("The sex-specific gene expression progression over age (top panel) is seen in parallel with the significance of the difference in expression between sexes (bottom panel).
                                                                                                     <br> For data privacy purposes, a light randomization was attributed to each scattered point along the age axis. The plotted regression lines were fitted using the original values.
                                                                                                     <br> The overall p-value, t-statistic and logFC (pink and blue large dots, for females and males, respectively) for the average age was computed using the ShARP-LM model over the entire age range.
                                                                                                     "),
                                                                               ),
                                                                               conditionalPanel(condition = "input.variable_Gene_Alteration == 'Int_Age_Gender'",
                                                                                                HTML("The sex-specific gene expression progression over age (top panel) is seen in parallel with the significance of the difference in age associated change of expression between sexes (bottom panel).
                                                                                                      <br> For data privacy purposes, a light randomization was attributed to each scattered point along the age axis. The plotted regression lines were fitted using the original values.
                                                                                                      <br> The overall p-value, t-statistic and logFC/year was computed using the ShARP-LM model over the entire age range.
                                                                                                    "),
                                                                               ),
                                                                               
                                                                               withSpinner(highchartOutput("Line_signficanceAlterationsvsAge", height = "650"), 
                                                                                           color = "#2C3E50", type = 5, size = 0.5))
                                                       )
                                              )
                                  )
                        ))
             ),
             tabPanel(title = "Tissue",
                      icon = icon("heartbeat"),
                      sidebarLayout(
                        sidebarPanel2(width = 3,
                                      pickerInput(inputId = "tissue_2", label = "Tissue", 
                                                  choices = c("Loading ..." = ""),
                                                  options = list(`live-search` = TRUE)),
                                      conditionalPanel(condition = "input.tissue_2 != 'All tissues' && input.tissue_2 != ''",
                                                       highchartOutput("Histogram_DistributionSample", height = "160px")),
                                      awesomeRadio(inputId = "variable_model_2", label = "Alterations associated with", inline = T,
                                                   choices = list("Age" = "Age", "Sex" = "Gender", "Age&Sex" = "Int_Age_Gender"),
                                                   checkbox = TRUE),
                                      conditionalPanel(condition = "input.tissue_tab == 'Enrichment' && input.tissue_2 != 'All tissues' && input.tissue_2 != ''",
                                                       awesomeRadio(inputId = "EnrichmentDataset", label = "Geneset", inline = T,
                                                                    choices = c("Reactome", "User-specified"), checkbox = T)),
                                      conditionalPanel(condition = "input.EnrichmentDataset == 'User-specified' && input.tissue_tab == 'Enrichment' && input.tissue_2 != 'All tissues' && input.tissue_2 != ''",
                                                       textAreaInput("geneListEnrichment", "List of genes",
                                                                     height = "300px", placeholder = "e.g.\nPCNA\nmki67\nglb1\nCDKN2A\nCDKN1A\nDEC1\ncdkn2b\nTNFRSF10D"),
                                                       numericInputIcon("manualEnrichmentThreshold", label = "Significance threshold (p-value)",
                                                                        value = 0.05, min = 0, max = 1),
                                                       actionButton("manualEnrichment", "Run")),
                                      conditionalPanel(condition = "input.EnrichmentDataset == 'Reactome' && input.tissue_tab == 'Enrichment'",
                                                       uiOutput("EnrichmentAllOrSpecific")),
                                      out = conditionalPanel(condition = "(input.peakClickedState == 'hover' || input.peakClickedState == '')  && input.tissue_tab == 'Peaks' ",
                                                             div(DT::dataTableOutput(outputId = "genePeakTable", height = "800px"), style = "font-size: 90%; width: 100%",server = F))
                        ),
                        mainPanel(width = 9,
                                  conditionalPanel(condition = "input.tissue_2 == 'All tissues' || input.tissue_2 == ''", 
                                                   column(width = 12, #Title different for each variable
                                                          conditionalPanel(condition = "input.variable_model_2 == 'Age'",
                                                                           tags$h4("Global age-related gene expression alterations over age across tissues"),
                                                          ),
                                                          conditionalPanel(condition = "input.variable_model_2 == 'Gender'",
                                                                           tags$h4("Global differences in gene expression between sexes over age across tissues"),
                                                          ),
                                                          conditionalPanel(condition = "input.variable_model_2 == 'Int_Age_Gender'",
                                                                           tags$h4("Global differences in age-related gene expression alterations between sexes over age across tissues"),
                                                          ),
                                                          withSpinner(highchartOutput("AgeWavesOrganism", height = "800px"), 
                                                                      color = "#2C3E50", type = 5, size = 0.5))),
                                  conditionalPanel(condition = "input.tissue_2 != 'All tissues' && input.tissue_2 != ''",
                                                   tabsetPanel(id = "tissue_tab", type = "pills",
                                                               tabPanel("Peaks",
                                                                        icon = icon("crosshairs"),
                                                                        column(width = 12,
                                                                               #Title different for each variable
                                                                               conditionalPanel(condition = "input.variable_model_2 == 'Age'",
                                                                                                tags$h4("Tissue-specific global age-related gene expression alterations over age"),
                                                                               ),
                                                                               conditionalPanel(condition = "input.variable_model_2 == 'Gender'",
                                                                                                tags$h4("Tissue-specific differences in gene expression between sexes over age"),
                                                                               ),
                                                                               conditionalPanel(condition = "input.variable_model_2 == 'Int_Age_Gender'",
                                                                                                tags$h4("Tissue-specific differences in age-related gene expression alterations between sexes over age"),
                                                                               ),
                                                                               HTML("Click on a dot to browse the altered genes (p.value/FDR <= 0.05) associated with the respective age-window linear model."),
                                                                               highchartOutput("SigvsAge", height = "400px")
                                                                        ),
                                                                        conditionalPanel(condition = "output.selectedGene != 'none' || isNaN(output.selectedGene)",
                                                                                         column(width = 12,
                                                                                                tags$h4("Gene expression alterations over age"),
                                                                                                HTML("For a given gene, its expression progression over age (bottom panel) is seen in parallel with its expression alteration's significance (top panel).
                                                                                                   <br> The overall p-value, t-statistic and logFC/year was computed using the entire age range."),
                                                                                                highchartOutput("Line_pvaluevsAge", height = "500px"))
                                                                        )
                                                               ),
                                                               
                                                               tabPanel("Enrichment",
                                                                        icon = icon("binoculars"),
                                                                        conditionalPanel(condition = "input.EnrichmentDataset == 'Reactome'",
                                                                                         column(width = 12,
                                                                                                #Title different for each variable
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Age'",
                                                                                                                 tags$h4("Tissue-specific age-related alterations in biological pathways over age"),
                                                                                                ),
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Gender'",
                                                                                                                 tags$h4("Tissue-specific differences in biological pathways between sexes over age"),
                                                                                                ),
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Int_Age_Gender'",
                                                                                                                 tags$h4("Tissue-specific differences in age-related alterations of biological pathways between sexes over age"),
                                                                                                ),
                                                                                                HTML("Gene Set Enrichment Analyses done on <a href='https://reactome.org/' target='_blank'>REACTOME</a> pathways (excluding those with fewer than 15 genes and more than 500).
                                                                                  <br>Pathways are gathered into families (together with those from <a href='https://www.genome.jp/kegg/' target='_blank'>KEGG</a> and level 3 <a href='http://geneontology.org/' target='_blank'>Gene Ontology</a> Biological Processes) based on the proportion of genes in common.
                                                                                                   <br>Only pathways with FDR <=0.05 for the respective age, whose enrichment's adjusted p.value <=0.05 and that belong to the top 1% of adjusted p.values are shown."),
                                                                                                withSpinner(combineWidgetsOutput("Heatmap_NESAgevsPathway", height = "600px"),
                                                                                                            color = "#2C3E50", type = 5, size = 0.5)),
                                                                                         tags$h4("Families of pathways"),
                                                                                         column(2, highchartOutput("Heatmap_family")),
                                                                                         column(10, 
                                                                                                tabsetPanel(type = "pills",
                                                                                                            tabPanel("WordCloud",
                                                                                                                     withSpinner(highchartOutput("wordcloud_family"), 
                                                                                                                                 color = "#2C3E50", type = 5, size = 0.5)),
                                                                                                            tabPanel("Pathways",
                                                                                                                     div(DT::dataTableOutput("listPathway_family"), style = "font-size: 90%; width: 100%", server = F)
                                                                                                            ))
                                                                                         )),
                                                                        conditionalPanel(condition = "input.EnrichmentDataset == 'User-specified' && input.tissue_tab == 'Enrichment'",
                                                                                         column(width = 12,
                                                                                                #Title different for each variable
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Age'",
                                                                                                                 tags$h4("Tissue-specific age-related alterations in user-specified gene set over age"),
                                                                                                ),
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Gender'",
                                                                                                                 tags$h4("Tissue-specific differences in user-specified gene set between sexes over age"),
                                                                                                ),
                                                                                                conditionalPanel(condition = "input.variable_model_2 == 'Int_Age_Gender'",
                                                                                                                 tags$h4("Tissue-specific differences in age-related alterations of user-specified gene set between sexes over age"),
                                                                                                ),
                                                                                                HTML("The significance of the user-specified geneset’s enrichment in differentially expressed genes is computed with Fisher’s exact tests.
                                                                                                   <br>Enter a geneset and a significance threshold value for differential expression and click on Run.
                                                                                                   <br>Click on a dot to access the contingency table associated with the respective age-window.<br>"),
                                                                                                withSpinner(highchartOutput("Line_FishertestEnrichmentManualvsAge"), 
                                                                                                            color = "#2C3E50", type = 5, size = 0.5),
                                                                                                uiOutput("geneList")))
                                                               )
                                                               
                                                   )                    
                                                   
                                  )
                                  
                                  
                        ))
             ),
             tabPanel("Module",
                      #icon = icon("object-ungroup"),
                      icon = icon("hubspot"),
                      mainPanel(width = 12,
                                tabsetPanel(id = "home", type = "pills",
                                            tabPanel("Overview",
                                                     icon = icon("comments"),
                                                     column(width = 12, style = "margin-bottom:50px",
                                                            HTML('<center><img src="moduleMethod.png" height="530"></center>'),
                                                            HTML("<font size='4'><br>For each tissue, a gene co-expression network was built based on the pairwise correlation in expression of all pairs of genes. 
                                 Groups of genes with similar expression across samples (called modules) are then identified by clustering analysis and examined to highlight enrichment in biological pathways, cell types markers or disease-associated genes.
                                        <br>A module's behaviour over age can then be inquired with the expression of its eigengene, representative of its gene expression profile.</font>"))
                                            ),
                                            tabPanel("Results",
                                                     icon = icon("chart-line"),
                                                     sidebarLayout(
                                                       sidebarPanel2(width = 3,
                                                                     pickerInput(inputId = "tissue_module", label = "Tissue", 
                                                                                 choices = c("Brain - Cortex", "Heart - Left Ventricle", "Muscle - Skeletal", "Whole Blood"),
                                                                                 options = list(`live-search` = TRUE)),
                                                                     pickerInput("selectedModule", "Module", choices = c("Loading ..." = ""), 
                                                                                 options = list(`live-search` = TRUE)),
                                                                     conditionalPanel(condition = "input.selectedModule == 'All modules' && input.modules == 'Expression'",
                                                                                      p("Enter a gene to check its belonging to a module"),
                                                                                      textInput("searchGeneInModule", label = "Gene", value = ""),
                                                                                      textOutput("geneInModules")),
                                                                     conditionalPanel(condition = "input.selectedModule == 'All modules' && input.modules == 'Diseases'",
                                                                                      awesomeRadio(inputId = "diseaseMethod", label = "Method", inline = TRUE, checkbox = TRUE, 
                                                                                                   choices = c("DisGeNET", "Manual"))
                                                                     ),
                                                                     conditionalPanel(condition = "input.selectedModule != 'All modules' && input.selectedModule != ''",
                                                                                      awesomeRadio(inputId = "moduleColoredBy", label = "Colored By", 
                                                                                                   choices = c("Loading ...")),
                                                                                      conditionalPanel(condition = "input.moduleColoredBy == 'Technicality' && input.selectedModule != 'All modules' && input.selectedModule != '' ", uiOutput("moduleColoredByTechnicality")),
                                                                                      awesomeRadio(inputId = "moduleShapedBy", label = "Shaped By", 
                                                                                                   choices = c("Loading ...")),
                                                                                      HTML("<small>*A curve is independently fitted for each of the positive and negative conditions if at least 10 samples are known for each of them.</small>")
                                                                     ),
                                                                     conditionalPanel(condition = "input.moduleShapedBy == 'Condition' && input.selectedModule != 'All modules' && input.selectedModule != '' ", uiOutput("moduleShapedByCondition")),
                                                                     out = uiOutput("modulesOutsideBar")
                                                       ),
                                                       mainPanel(width = 9,
                                                                 conditionalPanel(condition = "input.selectedModule == 'All modules' || input.selectedModule == ''", 
                                                                                  tabsetPanel(id = "modules", type = "pills",
                                                                                              tabPanel("Expression",
                                                                                                       icon = icon("chart-line"),
                                                                                                       column(width = 12, 
                                                                                                              tags$h4("Eigengene expression"),
                                                                                                              withSpinner(manipulateWidget::combineWidgetsOutput("Heatmap_LoessMEvsAge", height = "600px"), 
                                                                                                                          color = "#2C3E50", type = 5, size = 0.5))
                                                                                                       
                                                                                              ),
                                                                                              tabPanel("Cell types",
                                                                                                       icon = icon("vial"),
                                                                                                       column(width = 12, 
                                                                                                              tags$h4("Cell type enrichment"),
                                                                                                              HTML("The significance of each module’s enrichment in each cell type markers is computed with Fisher’s exact tests.
                                                                                               <br> The labelled values have odds ratio > 1 and p-value < 0.05. "),
                                                                                                              withSpinner(highchartOutput("Heatmap_FisherTest_cellComposition_All", height = "600px"), 
                                                                                                                          color = "#2C3E50", type = 5, size = 0.5))
                                                                                              ),
                                                                                              tabPanel(title = "Pathways",
                                                                                                       icon = icon("binoculars"),
                                                                                                       column(width = 12, 
                                                                                                              tags$h4("Functional enrichment"),
                                                                                                              HTML("The significance of each module’s enrichment in <a href='https://reactome.org/' target='_blank'>REACTOME</a> pathways is computed with Fisher’s exact tests.
                                                                                               <br>Pathways are gathered into families (together with those from <a href='https://www.genome.jp/kegg/' target='_blank'>KEGG</a> and level 3 <a href='http://geneontology.org/' target='_blank'>Gene Ontology</a> Biological Processes) based on the proportion of genes in common."),
                                                                                                              withSpinner(manipulateWidget::combineWidgetsOutput("Heatmap_moduleEnrichment_All", height = "600px"), 
                                                                                                                          color = "#2C3E50", type = 5, size = 0.5))
                                                                                                       
                                                                                              ),
                                                                                              tabPanel(title = "Diseases",
                                                                                                       icon = icon("briefcase-medical"),
                                                                                                       column(width = 12, 
                                                                                                              tags$h4("Disease enrichment"),
                                                                                                              conditionalPanel(condition = "input.selectedModule == 'All modules' && input.diseaseMethod == 'DisGeNET'",
                                                                                                                               HTML("The significance of the module’s enrichment in <a href='https://www.disgenet.org/search' 
                                                                                                           target='_blank'>DisGeNET</a> gene-disease associations (from the <code>CURATED</code> set) 
                                                                                                           is calculated with the disgenet2r package <a href='https://www.disgenet.org/static/disgenet2r/disgenet2r.html' 
                                                                                                           target='_blank'>(Piñero et al.)</a>.
                                                                                                           P-values were corrected for multiple testing with Benjamini-Hochberg’s FDR.
                                                                                                           Only diseases whose enrichment is significant (adjusted p-value < 0.05) in at least one module are displayed."),
                                                                                                                               withSpinner(highchartOutput("Heatmap_moduleDiseaseEnrichment_All", height = "600px"), 
                                                                                                                                           color = "#2C3E50", type = 5, size = 0.5)
                                                                                                              ),
                                                                                                              conditionalPanel(condition = "input.selectedModule == 'All modules' && input.diseaseMethod == 'Manual'",
                                                                                                                               HTML("The enrichment is based on the gene-disease association from <a href='https://www.disgenet.org/search' target='_blank'>DisGeNET</a> and calculated with Fisher tests.
                                                                                                                <br>P-values were corrected for multiple testing with Benjamini-Hochberg's FDR.
                                                                                                                Only diseases with at least 20 and up to 500 associated genes were considered. 
                                                                                                                Only diseases whose enrichment is significant (adjusted p-value < 0.05) in at least one module are displayed."),
                                                                                                                               withSpinner(highchartOutput("Heatmap_moduleDiseaseEnrichmentManual_All", height = "600px"), 
                                                                                                                                           color = "#2C3E50", type = 5, size = 0.5)
                                                                                                              )
                                                                                                       )
                                                                                              )
                                                                                              # ,
                                                                                              # tabPanel(title = "Diseases-Manual",
                                                                                              #          icon = icon("briefcase-medical"),
                                                                                              #          column(width = 12, 
                                                                                              #                 tags$h4("Disease enrichment"),
                                                                                              #                 HTML("The enrichment is based on the gene-disease association from <a href='https://www.disgenet.org/search' target='_blank'>DisGeNET</a> and calculated with Fisher tests."),
                                                                                              #                 withSpinner(highchartOutput("Heatmap_moduleDiseaseEnrichmentManual_All", height = "600px"), 
                                                                                              #                             color = "#2C3E50", type = 5, size = 0.5))
                                                                                              # )
                                                                                  )
                                                                 ),
                                                                 conditionalPanel(condition = "input.selectedModule != 'All modules' && input.selectedModule != ''", 
                                                                                  column(width = 12, 
                                                                                         tags$h4("Eigengene expression"),
                                                                                         HTML("For data privacy purposes, a light randomization was attributed to each scattered point along the age axis. The plotted regression line (black) was fitted using the original values."),
                                                                                         withSpinner(highchartOutput("Scatter_MEvsAge", height = "400px"), color = "#2C3E50", type = 5, size = 0.5)),
                                                                                  tabsetPanel(id = "module-specific", type = "pills",
                                                                                              tabPanel("Cell types",
                                                                                                       icon = icon("vial"),
                                                                                                       column(width = 12, 
                                                                                                              uiOutput("referenceMarker"),
                                                                                                              manipulateWidget::combineWidgetsOutput("Heatmap_FisherTest_cellComposition_One", height = "700px")
                                                                                                       )),
                                                                                              tabPanel("Pathways",
                                                                                                       icon = icon("binoculars"),
                                                                                                       column(width = 12, 
                                                                                                              HTML("<small>Pathways significantly associated with the module are highlighted (<font color='tomato'>p-value <= 0.05</font> in red and <font color='tomato'><b> adjusted p-value <= 0.05</b></font> in bold).
                                                                                               <br>The significance of the module’s enrichment in <a href='https://reactome.org/' target='_blank'>REACTOME</a> pathways is computed with Fisher’s exact tests.
                                                                                               <br>P-values were corrected for multiple testing with Benjamini-Hochberg’s FDR.</small>"),
                                                                                                              div(DT::dataTableOutput("DT_moduleEnrichment_One"), style = "font-size: 80%; width: 100%", server = F))),
                                                                                              tabPanel("Diseases-DisGeNET",
                                                                                                       icon = icon("briefcase-medical"),
                                                                                                       column(width = 12,
                                                                                                              HTML("<small>Diseases significantly associated with the module are highlighted (<font color='tomato'>p-value <= 0.05</font> in red and <font color='tomato'><b> adjusted p-value <= 0.05</b></font> in bold).
                                                                                               <br>The significance of the module’s enrichment in <a href='https://www.disgenet.org/search' target='_blank'>DisGeNET</a> gene-disease associations (from the <code>CURATED</code> set) is calculated with the disgenet2r package <a href='https://www.disgenet.org/static/disgenet2r/disgenet2r.html' target='_blank'>(Piñero et al.)</a>.
                                                                                               <br>P-values were corrected for multiple testing with Benjamini-Hochberg’s FDR. Only diseases whose enrichment is significant (adjusted p-value < 0.05) in at least one module are displayed.</small>"),
                                                                                                              div(DT::dataTableOutput("DT_diseaseEnrichment"), style = "font-size: 80%; width: 100%", server = F))),
                                                                                              tabPanel("Diseases-Manual",
                                                                                                       icon = icon("briefcase-medical"),
                                                                                                       column(width = 12,
                                                                                                              HTML("<small>Diseases significantly associated with the module are highlighted (<font color='tomato'>p-value <= 0.05</font> in red and <font color='tomato'><b> adjusted p-value <= 0.05</b></font> in bold).
                                                                                               <br>The significance of the module’s enrichment in <a href='https://www.disgenet.org/search' target='_blank'>DisGeNET</a> gene-disease associations is calculated with Fisher's exact tests. Only diseases with at least 20 associated genes and less than 500 genes were used in the analysis.
                                                                                               <br>P-values were corrected for multiple testing with Benjamini-Hochberg’s FDR. Only diseases whose enrichment is significant (adjusted p-value < 0.05) in at least one module are displayed.</small>"),
                                                                                                              div(DT::dataTableOutput("DT_diseaseEnrichmentManual"), style = "font-size: 80%; width: 100%", server = F)))
                                                                                  )
                                                                                  
                                                                                  
                                                                 )
                                                       )
                                                     ))
                                )))#,
             #tabPanel("Tutorial", 
             # icon = icon("file-text"),
             # tags$iframe(style="height:800px; width:100%", src="voyAGEr-WebAppTutorial.html"))
  )
  
  
))
