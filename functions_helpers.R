
color2 <- function()
{
  color2 <- c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
              "#7FC97F", "#BEAED4", "#FDC086", "#E6AB02", "#386CB0", "#A6761D", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", 
              "#7570B3", "#E7298A", "#66A61E", "#FFFF99", "#BEBADA")
  color2
}

Heatmap_log10FDRvsAge <- function(variable, AgeWaves)
{
   
  tmp <- AgeWaves[[variable]]
  tmp[is.na(tmp)] <- 0#3 
  tmp[tmp < 0] <- 0
  tmp[sapply(tmp, is.infinite)] <- NA
  maxtmp <- max(tmp, na.rm = T)
  tmp[is.na(tmp)] <- maxtmp + 2 # to ease visualization, use Inf values as the maximum value + 2 (significant for sure)
 
  #tmp <- as.matrix(tmp[order(row.names(tmp)), ])
  
  #tmp[tmp < 0] <- 0
  i <- (hclust(dist(tmp)))
  #browser()
  tmp <- as.matrix(tmp[i$order, ])
  tmp <- melt(tmp)
  colnames(tmp) <- c("tissue", "age", "-log10(FDR)")
  
  ClickFunction <- JS("function(event) {Shiny.onInputChange('tissueClicked', event.point.tissue)  }")
  ClickFunction2 <- JS("function(event) {console.log(event)  }")
  
  g <- hchart(tmp, type = "heatmap", boost = list(useGPUTranslations = T),
              hcaes(x = age, y = tissue, value = `-log10(FDR)`)) %>% 
    hc_colorAxis(stops = color_stops(100, rev(magma(100))),
                 min = 0,
                 reversed = F) %>%
    hc_legend(layout = "vertical", verticalAlign = "top",
              align = "right", title = list(text = "-log<sub>10</sub>(FDR) <br/>  <i><small> Proportion altered genes <i/><small/>", style = list(fontSize = '16px')), useHTML = TRUE) %>%
    hc_xAxis(title = list(text = "Age (years)", style = list(fontSize = '18px')),
             labels = list(style = list(fontSize = "14px"))) %>%
    hc_yAxis(title = list(text = "Tissue", style = list(fontSize = '18px')), reversed = TRUE) %>%
    hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>Tissue</b>: ' +  this.series.yAxis.categories[this.point.y] + ',<br><b> -log<sub>10</sub>(FDR)</b>: ' + Highcharts.numberFormat(this.point.value, 2);}"), useHTML = T) %>%
    hc_plotOptions(series = list(events = list(click = ClickFunction)))
  
  return(g)
}

Line_sigGenesvsAge <- function(tissue, variable, peak)
{
  tmp <- peak[peak$tissue == tissue & peak$variable == variable,]
  tmp$pvalue[tmp$pvalue == 0] <- 1e-6
  tmp$log10p <- -log10(tmp$pvalue)
  
  
  ClickFunction <- JS("function(event) {Shiny.onInputChange('peakClicked', event.point.age)
                      Shiny.onInputChange('peakClickedState', event.point.state)
                      Shiny.onInputChange('peakClickedVariable', event.point.series.name)}")
  ClickFunction2 <- JS("function(event) {console.log(event)  }")
  
  if (max(tmp$log10p)< -log10(0.05)){
    max_scale <- -log10(0.05) 
  } else {
    max_scale <- max(tmp$log10p)
  }
  
  tmp$mycolor <- circlize::colorRamp2(breaks = c(0, max_scale/2, max_scale), colors = c("#dfe59a", "#9dca94", "#346875"))(tmp$log10p)
  
  g <- hchart(tmp, "scatter", hcaes(x = age, y = PercSigGene, 
                                    color = mycolor), name = variable) %>%
    hc_plotOptions(series = list(allowPointSelect = TRUE,
                                 lineWidth = 6,
                                 color = "black",
                                 marker = list(radius = 8, enabled = T, 
                                               states = list(select = list(fillColor = "tomato", lineWidth = 0, radius = 14))),
                                 events = list(click = ClickFunction))) %>%
    hc_colorAxis(stops = list(list(0, "#dfe59a"), list(0.5, "#9dca94"), list(1, "#346875")),
                 min = 0, max = max(tmp$log10p)) %>% 
    hc_legend(layout = "horizontal", reversed = T,
              align = "center", title = list(text = "-log<sub>10</sub>(FDR)"), useHTML = TRUE) %>%
    hc_xAxis(title = list(text = "Age (years)"), min = 20, max = 70) %>%
    hc_yAxis(title = list(text = "% altered genes")) %>%
    hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.age + ' y.o., <br><b>%alt genes</b>: ' + Highcharts.numberFormat(this.point.PercSigGene, 2) + ',<br><b> -log<sub>10</sub>(FDR)</b>: ' + Highcharts.numberFormat(this.point.options.log10p, 2);}"), useHTML = TRUE) 
  
  return(g)
}

Line_pvaluevsAge <- function(gene, pvalueData, geneData, allAges, abline, tissue, coloredBy)
{
  #GE info
  GE <- geneData
  GE <- GE[GE$tissue == tissue,] 
  GE$jittered_age <- jitter_ages(GE$age)
  #statsAllAges_text <- paste0("Overall changes: p-value = ", round(allAges$pvalue,4), " | t-statistic = ", round(allAges$tvalue,2))
  statsAllAges_text <- paste0("Overall changes: p-value = ", round(allAges$pvalue,4), " | t-statistic = ", round(allAges$tvalue,2), " | logFC/year = ", round(abline$a,4))
   
  if (coloredBy == "Age")
  {
    
    fit <- loess(expression ~ age, data = GE)
    fit <- plyr::arrange(generics::augment(fit), age)
    
    #Pvalue info
    tmp <- melt(as.matrix(pvalueData[rownames(pvalueData) == gene,]))
    colnames(tmp) <- c("gene", "Age", "p-value")
    tmp$log10p <- -log10(tmp$`p-value`)
    
    g <-   highchart() %>% 
      hc_yAxis_multiples(
        list(lineWidth = 2, title = list(useHTML = TRUE, text = "log<sub>10</sub>(p-value)", style = list(color = "#9A839A")), lineColor = "#9A839A", labels = list(style = list(color = "#9A839A")), height = "50%",
             plotLines = list(list(value = 1.3, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.01")))),
        list(showLastLabel = FALSE, opposite = TRUE, title = list(text = "GE (logCPM)", style = list(color = "yellowgreen")), height = "50%", top = "50%",
             lineColor = "yellowgreen", lineWidth = 2,
             labels = list(style = list(color = "yellowgreen")))
      ) %>% 
      hc_add_series(data = list_parse2(GE[, c("jittered_age", "expression")]), type = "scatter", showInLegend = FALSE,
                    yAxis = 1, color = "yellowgreen", marker = list(radius = 4, symbol = "round")) %>% 
      # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), yAxis = 1,
      #               type = "arearange", showInLegend = FALSE, color = "lightgrey") %>%
      hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "black") %>% 
      hc_add_series(data = list_parse2(tmp[, c("Age", "log10p")]), type = "line", showInLegend = FALSE, name = "log10p",
                    color = "#9A839A", marker = list(radius = 6, enabled = T, color = "black"), lineWidth = 5) %>%
      hc_xAxis(title = list(text = "Age (years)"), lineColor = "grey", min = 20, max = 70) %>%
      hc_tooltip(formatter = JS("function(){
                                if(this.point.series.name == 'log10p'){ 
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);}
                                else{return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}
                                }"), useHTML = TRUE)
    
  } else
  {
    fit_male <- loess(expression ~ age, data = GE[GE$sex == 1,])
    fit_male <- plyr::arrange(generics::augment(fit_male), age)
    
    fit_female <- loess(expression ~ age, data = GE[GE$sex == 2,])
    fit_female <- plyr::arrange(generics::augment(fit_female), age)
    
    #Pvalue info
    tmp <- melt(as.matrix(pvalueData[rownames(pvalueData) == gene,]))
    colnames(tmp) <- c("gene", "Age", "p-value")
    tmp$log10p <- -log10(tmp$`p-value`)
    
    g <-   highchart() %>% 
      hc_yAxis_multiples(
        list(lineWidth = 2, title = list(useHTML = TRUE, text = "log<sub>10</sub>(p-value)", style = list(color = "#9A839A")), lineColor = "#9A839A", labels = list(style = list(color = "#9A839A")), height = "50%",
             plotLines = list(list(value = 1.3, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.01")))),
        list(showLastLabel = FALSE, opposite = TRUE, height = "50%", top = "50%",
             title = list(text = "GE (logCPM)", style = list(color = "grey")), 
             lineColor = "grey", lineWidth = 2,
             labels = list(style = list(color = "grey")))
      ) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 1, c("jittered_age", "expression")]), type = "scatter", name = "Male",
                    yAxis = 1, color = "cornflowerblue", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 2, c("jittered_age", "expression")]), type = "scatter", name = "Female",
                    yAxis = 1, color = "hotpink", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "cornflowerblue") %>% 
      hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "hotpink") %>% 
      hc_add_series(data = list_parse2(tmp[, c("Age", "log10p")]), type = "line", showInLegend = FALSE, name = "log10p",
                    color = "#9A839A", marker = list(radius = 6, enabled = T, color = "black"), lineWidth = 5) %>%
      hc_xAxis(title = list(text = "Age (years)"), lineColor = "grey", min = 20, max = 70) %>%
      hc_tooltip(formatter = JS("function(){
                                if(this.point.series.name == 'log10p'){ 
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);}
                                else{return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}
                                }"), useHTML = TRUE)
    
  }
  
  g <- g%>% 
    hc_subtitle(text = statsAllAges_text, style = list(color = "black") ) 
  
  return(g)
}

jitter_ages <- function(agevec){
  
  set.seed(123456)
  
  # Create age bins, compatible to those used in the public metadata
  bins <- cut(agevec, breaks = seq(20, 70, 5), include.lowest = TRUE)
  
  # Compute bin midpoints (mean of each bin)
  bin_midpoints <- as.numeric(tapply(agevec, bins, function(x) mean(range(x))))
  
  # Replace values with bin midpoints
  new_values <- bin_midpoints[(bins)]
  
  # Add jitter within 5-year window
  jittered_values <- jitter(new_values, amount = 2.5) # half the bin width 
  # Ensure values stay within 20 and 70
  jittered_values[jittered_values < 20] <- 20
  jittered_values[jittered_values > 70] <- 70
  
  return(jittered_values)
  
}

scatter_GEvsAge <- function(data, tissue, colored, shaped, donorCondition, gene, geneInfo)
{
  #tmp <- dbGetQuery(DBconnection, paste("SELECT * FROM", gene, sep = " "))
  tmp <- data
  tmp <- tmp[tmp$tissue == tissue,]
  tmp$jittered_age <- jitter_ages(tmp$age) 
  
  #To zoom the plot onto a selected area
  #Removed so age can't be found to respect GTEx rules
  #To remove the possibility to subset some donor GE, comment , zoomType = 'xy' , events = list(selection = SelectFunction) in the functions
  SelectFunction <- JS("function(event) {
                              var x_axis_min,x_axis_max,y_axis_min,y_axis_max;

                              if (event.xAxis) {
                                x_axis_min = Highcharts.numberFormat(event.xAxis[0].min, 2),
                                x_axis_max = Highcharts.numberFormat(event.xAxis[0].max, 2);
                              }else{
                                x_axis_min = 'reset',
                                x_axis_max = 'reset';
                              }
                              if (event.yAxis) {
                                y_axis_min = Highcharts.numberFormat(event.yAxis[0].min, 2),
                                y_axis_max = Highcharts.numberFormat(event.yAxis[0].max, 2);
                              }else{
                                y_axis_min = 'reset',
                                y_axis_max = 'reset';
                              }
                              Shiny.onInputChange('Selectedxmin', x_axis_min);
                              Shiny.onInputChange('Selectedxmax', x_axis_max);
                              Shiny.onInputChange('Selectedymin', y_axis_min);
                              Shiny.onInputChange('Selectedymax', y_axis_max);
    }")
  
  
  ClickFunction <- JS("function(event) {Shiny.onInputChange('ClickedIndex', event.point.index)
                                        Shiny.onInputChange('ClickedStatus', event.point.state)}")
  
  ClickFunction2 <- JS("function(event) {console.log(event)  }")
  
  
  if (colored == "none" & shaped == "none")
  {
    fit <- loess(expression ~ age, data = tmp)
    fit <- plyr::arrange(generics::augment(fit), age)
    
    g <- highchart() %>%
      hc_chart(type = 'scatter'
               , zoomType = 'xy'
               , events = list(selection = SelectFunction)
      ) %>%
      hc_add_series(data = list_parse2(tmp[, c("jittered_age", "expression")]), showInLegend = FALSE) %>%
      # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), 
      #               type = "arearange", showInLegend = FALSE) %>%
      hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                    type = "line", showInLegend = FALSE) %>%
      hc_plotOptions(arearange = list(color = "gainsboro", opacity = 0.005),
                     line = list(lineWidth = 7, color = "black"),
                     scatter = list(marker = list(radius = 5, symbol = "line"), color = "yellowgreen"),
                     series = list(allowPointSelect = TRUE, 
                                   events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  } else if (colored == "Sex" & shaped == "none")
  {
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      fit_male <- loess(expression ~ age, data = tmp[tmp$sex == 1,])
      fit_male <- plyr::arrange(generics::augment(fit_male), age)
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      fit_female <- loess(expression ~ age, data = tmp[tmp$sex == 2,])
      fit_female <- plyr::arrange(generics::augment(fit_female), age)
    }
    g <- highchart() %>%
      hc_chart(type = 'scatter'
               , zoomType = 'xy'
               , events = list(selection = SelectFunction)
      ) 
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1, c("jittered_age", "expression")]), 
                      showInLegend = TRUE, color = "cornflowerblue", name = "Male") 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2, c("jittered_age", "expression")]), 
                      showInLegend = TRUE, color = "hotpink", name = "Female")
    }
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "cornflowerblue") 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "hotpink") 
    }
    g <- g%>%
      hc_plotOptions(line = list(lineWidth = 7),
                     scatter = list(marker = list(radius = 5, symbol = "circle", fillOpacity = 0.4)),
                     series = list(allowPointSelect = TRUE, 
                                   events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  }  else if (colored == "none" & shaped != "none")
  {
    tmp$condition <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), match(shaped, colnames(donorCondition))]
    
    #Fitting of a regression per condition if the number of positive/negative case is higher than 10
    if (length(which(tmp$condition %in% c(0,1))) == 0)
    {
      fit <- loess(expression ~ age, data = tmp)
      fit <- plyr::arrange(generics::augment(fit), age)
      
      g <- highchart() %>%
        hc_chart(type = 'scatter') %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)),
        #               type = "arearange", showInLegend = FALSE,
        #               color = "gainsboro", opacity = 0.005) %>%
        hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "black") %>%
    	hc_legend(title = list(text = shaped,
    						   style = list(fontSize = '16px'))) %>%
        hc_plotOptions(series = list(allowPointSelect = TRUE, 
                                     events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
      
    }else if (min(table(tmp$condition)[names(table(tmp$condition)) %in% c(0,1)]) >= 10)
    {
      fit_positive <- loess(expression ~ age, data = tmp[tmp$condition == 1,])
      fit_positive <- plyr::arrange(generics::augment(fit_positive), age)
      
      fit_negative <- loess(expression ~ age, data = tmp[tmp$condition == 0,])
      fit_negative <- plyr::arrange(generics::augment(fit_negative), age)
      
      median_positive <- data.frame(20:70, rep(median(tmp$expression[which(tmp$condition == 1)]), length(20:70)))
      median_negative <- data.frame(20:70, rep(median(tmp$expression[which(tmp$condition == 0)]), length(20:70)))
      
      g <- highchart() %>%
        hc_chart(type = 'scatter'
                 , zoomType = 'xy'
                 , events = list(selection = SelectFunction)
        ) %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 1, c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "triangle"), color = "rgb(60,179,113)", name = "Positive") %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 0, c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "cross"), color = "rgb(255,99,71)", name = "Negative") %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        hc_add_series(data = list_parse2(data.frame(fit_positive$age, fit_positive$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "rgb(60,179,113)") %>%
        hc_add_series(data = list_parse2(data.frame(fit_negative$age, fit_negative$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "rgb(255,99,71)") %>%
        hc_add_series(data = list_parse2(median_positive), type = "line", name = "median - Positive",
                      marker = list(enabled = F),
                      color = "rgb(60,179,113)", dashStyle = "Dash") %>%
        hc_add_series(data = list_parse2(median_negative), type = "line", name = "median - Negative",
                      marker = list(enabled = F),
                      color = "rgb(255,99,71)", dashStyle = "Dash") %>%
    	hc_legend(title = list(text = shaped,
    						   style = list(fontSize = '16px'))) %>%
        hc_plotOptions(series = list(allowPointSelect = TRUE, 
                                     events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
      
      
    } else
    {
      
      fit <- loess(expression ~ age, data = tmp)
      fit <- plyr::arrange(generics::augment(fit), age)
      
      g <- highchart() %>%
        hc_chart(type = 'scatter'
                 , zoomType = 'xy', events = list(selection = SelectFunction)
        ) %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 1, c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "triangle"), color = "rgb(60,179,113)", name = "Positive") %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 0, c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "cross"), color = "rgb(255,99,71)", name = "Negative") %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "expression")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), 
        #               type = "arearange", showInLegend = FALSE,
        #               color = "gainsboro", opacity = 0.005) %>%
        hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "black") %>%
    	hc_legend(title = list(text = shaped,
    						   style = list(fontSize = '16px'))) %>%
        hc_plotOptions(series = list(allowPointSelect = TRUE, 
                                     events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
    }
    
  } else if (colored == "Sex" & shaped != "none")
  {
    tmp$condition <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), match(shaped, colnames(donorCondition))]
    
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      fit_male <- loess(expression ~ age, data = tmp[tmp$sex == 1,])
      fit_male <- plyr::arrange(generics::augment(fit_male), age)
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      fit_female <- loess(expression ~ age, data = tmp[tmp$sex == 2,])
      fit_female <- plyr::arrange(generics::augment(fit_female), age)
    }
    g <- highchart() %>%
      hc_chart(type = 'scatter'
      , zoomType = 'xy', events = list(selection = SelectFunction))
    #Male and conditions
    if (nrow(tmp[tmp$sex == 1,])!=0)
    { 
      g <- g %>%  
        
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & tmp$condition == 1, c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Pos",
                      marker = list(radius = 7, symbol = "triangle")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & tmp$condition == 0, c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Neg",
                      marker = list(radius = 7, symbol = "cross")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & !(tmp$condition %in% c(0, 1)), c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Unknown",
                      marker = list(radius = 5, symbol = "circle")) 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    { 
      g <- g %>%
        #Female and conditions
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & tmp$condition == 1, c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Pos",
                      marker = list(radius = 7, symbol = "triangle")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & tmp$condition == 0, c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Neg",
                      marker = list(radius = 7, symbol = "cross")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & !(tmp$condition %in% c(0, 1)), c("jittered_age", "expression")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Unknown",
                      marker = list(radius = 5, symbol = "circle")) 
    }
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g %>%
        #Loess fitting for both sex
        hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "cornflowerblue", lineWidth = 7) 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "hotpink", lineWidth = 7) 
    }
    g <- g%>%
      #Fake series to create legend
      hc_add_series(type = "scatter", name = "Positive", color=  "grey", 
                    marker = list(symbol = "triangle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Negative", color=  "grey", 
                    marker = list(symbol = "cross", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Unknown", color=  "grey", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Male", color=  "cornflowerblue", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Female", color=  "hotpink", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      hc_legend(title = list(text = shaped,
    						 style = list(fontSize = '16px'))) %>%
      #Options for point click and selection
      hc_plotOptions(series = list(allowPointSelect = TRUE, 
                                   events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  }
  g <- g %>%
    hc_xAxis(title = list(text = "Age (years)"), min = 20, max = 70) %>%
    hc_yAxis(title = list(text = "Gene Expression (logCPM)")) %>%
    hc_tooltip(formatter = JS("function(){return ' <b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}"), useHTML = TRUE) %>%
    #hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}"), useHTML = TRUE) %>%
    #hc_title(text = gene) %>% 
    hc_subtitle(text = geneInfo) 
  
  return(g)
}  

datatable_donorInfo <- function(geneData, tissue, selectedPoint, xmin, xmax, ymin, ymax, indexPoint, donorCondition, conditionID)
{
  tmp <- geneData[geneData$tissue == tissue,]
  tmp$sex[tmp$sex == 1] <- "M"
  tmp$sex[tmp$sex == 2] <- "F"
  
  if (selectedPoint == "hover" | selectedPoint == "") #one point selected
  {
    tmp <- tmp[indexPoint,]
  } else if (is.numeric(xmin) | is.numeric(xmax) | is.numeric(ymin) | is.numeric(ymax))
  {
    tmp <- tmp[tmp$expression >= ymin & tmp$expression <= ymax & tmp$age >= xmin & tmp$age <= xmax,]
  } else 
  {
    tmp <- tmp
  }
  
  tmp <- cbind(tmp, donorCondition[match(tmp$samp_id, donorCondition$SAMPID), -1])
  tmp[tmp == 1] <- as.character(tags$i(
    class = "fa fa-check", 
    style = "color: rgb(60,179,113)"))
  tmp[tmp == 0] <- as.character(tags$i(
    class = "fa fa-times", 
    style = "color: rgb(255,99,71)"))
  tmp[is.na(tmp)] <- as.character(tags$i(
    class = "fa fa-question", 
    style = "color: rgb(169,169,169)"))
  tmp[tmp == 96] <- as.character(tags$i(
    class = "fa fa-question", 
    style = "color: rgb(169,169,169)"))
  tmp[tmp == 97] <- as.character(tags$i(
    class = "fa fa-question", 
    style = "color: rgb(169,169,169)"))
  tmp[tmp == 99] <- as.character(tags$i(
    class = "fa fa-question", 
    style = "color: rgb(169,169,169)"))
  
  tmp <- tmp[, c("samp_id", "age", "expression", "sex", colnames(donorCondition)[-1])]
  tmp$expression <- round(tmp$expression, 2)
  tmp$age <- as.character(cut(tmp$age, breaks = seq(20,70, by = 5), include.lowest = T))
  
  DT::datatable(tmp[, -which(colnames(tmp) == "samp_id")],
                escape = F,
                caption = htmltools::tags$caption(style = 'caption-side: top; text-align: center;',
                								  'Positive: ',
                                                  tags$i(class = "fa fa-check", style = "color: rgb(60,179,113)"),
                                                  '- Negative: ',
                                                  tags$i(class = "fa fa-times", style = "color: rgb(255,99,71)"),
                                                  '- Unknown: ',
                                                  tags$i(class = "fa fa-question", style = "color: rgb(169,169,169)")
                ),
                rownames = tmp$samp_id,
                #extensions = 'Buttons', no integrated downlaod bttn
                extensions = "FixedColumns", #freeze the 4 first columns
                selection = "none",
                colnames = c("Age interval (years)", "Gene Expression (logCPM)", "Sex", 
                             # c("Cytomegalovirus", "EB Virus", "Hepatitis B 1.", "Hepatitis B 2.", "Hepatitis C", "Abnormal WBC", "Alzheimer's", "Arthritis", "Ascites", "Asthma", "Infections", "Past Blood Donations Denied",
                             #                         "Cancer 1.", "Cancer 2.", "Respiratory Disease 1.", "Cocaine", "Respiratory Disease 2.", "Cerebrovascular Disease", "Dialysis", "Dementia", "Depression", "Heart attack", "Ischemic Heart Disease",
                             #                         "Heart Disease", "Hypertension", "Liver disease", "Open wounds", "Transplant", "Osteomyelitis", "Pill abuse", "Pneumonia 1.", "Pneumonia 2.", "+ blood culture", "Rheumatoid Arthritis",
                             #                         "Renal Failure", "Schizophrenia", "Sepsis", "Steroid", "Seizures", "Diabetes I", "Diabetes II", "Toxics", "Loss weight")
                             as.character(conditionID$ID[match(colnames(donorCondition)[-1], conditionID$description)])
                ),
                options = list(scrollX = T
                				, fixedColumns = list(leftColumns = 4) #freeze the 4 first columns
                				, autoWidth = TRUE
                               # ,
                               # columnDefs = list(list(targets = 0, render = JS(
                               #   "function(data, type, row, meta) {",
                               #   "return type === 'display' && data.length > 5 ?",
                               #   "'<span title=\"' + data + '\">' + data.substr(0, 5) + '...</span>' : data;",
                               #   "}")
                               # ))
                ),
                callback = JS("var tips = ['GTEx sample ID',
              'Age of the donor by quinquennium',
              'Gene Expression in logCPM',
              'Sex of the donor',
              'Cytomegalovirus Total Antibody blood test',
              'Epstein-Barr Virus IgG Antibody Blood test',
              'Hepatitis B core Antibody IgM',
              'Hepatitis B core Antibody Total',
              'Hepatitis C virus Antibody testing',
              'Abnormal White Blood Count',
              'Alzheimers OR Dementia',
              'Arthritis',
              'Ascites',
              'Asthma',
              'Bacterial Infections (including septicemia (bacteria in the blood), meningococcal disease, staphylococcal infection, streptococcus, sepsis)',
              'Past Blood Donations Denied',
              'Cancer Diagnosis in the last 5y',
              'History Of Non Metastatic Cancer',
              'Chronic Lower Respiratory Disease',
              'Cocaine Use In 5y',
              'Chronic Respiratory Disease (Chronic Obstructive Pulmonary Syndrome (COPD) OR Chronic Lower Respiratory Disease (CLRD) (chronic bronchitis, emphysema, asthma))',
              'Cerebrovascular Disease (stroke, TIA, embolism, aneurysm, other circulatory disorder affecting the brain)',
              'Dialysis Treatment',
              'Dementia With Unknown Cause',
              'Major depression (unipolar depression, major depressive disorder)',
              'Heart attack, acute myocardial infarction, acute coronary syndrome',
              'Ischemic Heart Disease (coronary artery disease (CAD), coronary heart disease, ischemic cardiomyopathy)',
              'Heart Disease',
              'Hypertension',
              'Liver Disease (liver abscess, failure, fatty liver syndrome, inherited liver insufficiency, acute/chronic hepatic insufficiency, necrobacillosis, rupture)',
              'Open Wounds',
              'Received Tissue Organ Transplant',
              'Osteomyelitis',
              'Prescription Pill Abuse',
              'Pneumonia',
              'Pneumonia (acute respiratory infection affecting the lungs)',
              'Positive Blood Cultures',
              'Rheumatoid Arthritis',
              'Renal Failure',
              'Schizophrenia',
              'Documented Sepsis',
              'Long Term Steroid Use',
              'Unexplained Seizures',
              'Diabetes mellitus type 1 (IDDM, formerly juvenile diabetes)',
              'Diabetes mellitus type II (NIDDM, adult onset diabetes)',
              'Exposure To Toxics',
              'Unexplained Weight Loss (Information)'],
                            header = table.columns().header();
                            for (var i = 0; i < tips.length; i++) {$(header[i]).attr('title', tips[i]);}"), 
  ) %>%
    DT::formatStyle("sex", backgroundColor = DT::styleEqual(c("M", "F"), c('aliceblue', 'mistyrose')))
}

Heatmap_LoessGEvsAge <- function(data)
{
  #a <- dbGetQuery(DBconnection, paste("SELECT * FROM", "A1BG_AS1", sep = " "))
  a <- data
  b <- data.frame()
  for (k in unique(a$tissue))
  {
    fit <- a[a$tissue == k,]
    fit <- loess(expression ~ age, data = fit)
    b <- rbind(b, data.frame(age = 20:70, expression = predict(fit, 20:70), tissue = k))
  }
  b <- dcast(b, tissue~age, value.var = "expression")
  rownames(b) <- b$tissue
  b <- b[,-1]
  b <- t(scale(t(b), scale = T))
  i <- (hclust(dist(b)))
  
  b <- b[i$order, ]
  #b[is.na(b)] <- -100
  b <- melt(b)
  colnames(b) <- c("tissue", "age", "expression")
 
  g <- hchart(b, type = "heatmap", hcaes(x = age, y = tissue, value = expression)) %>%
    hc_colorAxis(stops = list(list(0.3, "dodgerblue"), list(0.5, "white"), list(0.7, "tomato")),
                 min = -2, max = 2,
                 reversed = F) %>%
    hc_tooltip(formatter = JS("function(){if (!this.point.options.expression) 
                            {return '<b>Age</b>: ' + this.point.age + ' y.o., <br><b>Tissue</b>: ' + this.point.tissue + ',<br><b>scaled GE </b>: No data'}
  else
  {return '<b>Age</b>: ' + this.point.age + ' y.o., <br><b>Tissue</b>: ' + this.point.tissue + ',<br><b>scaled GE</b>: ' + Highcharts.numberFormat(this.point.options.expression, 2);}
                            }")) %>%
    hc_legend(title = list(text = "<center>Gene Expression<br>(Z-scores)</center>", style = list(fontSize = '16px')),align = "right", layout = "vertical", verticalAlign = "top",reversed = TRUE) %>%
    hc_xAxis(title = list(text = "Age (years)", style = list(fontSize = '18px')),
             labels = list(style = list(fontSize = "14px"))) %>%
    hc_yAxis(title = list(text = "Tissue", style = list(fontSize = '18px'))) 
  
  
  
  return(g)
}

Heatmap_NESAgevsPathway <- function(peak, tissue, variable, enrichment, Cluster_Reactome, ALL = TRUE)
{
  
  if (ALL == TRUE)
  {
    #Pick pathwyas that are signifcant for the peak to reduce their number and make the plot more easily readable
    temp <- peak[peak$tissue == tissue & peak$variable == variable,]
    i <- enrichment$pathway[which(enrichment$padj <= 0.05 
                                  & enrichment$padj <= quantile(enrichment$padj, 0.01) 
                                  & enrichment$age %in% temp$age[temp$sig == "sig"])] 
    i <- unique(i)
    while (identical(i, character(0))) #in some cases the previous filtering of intersting pathways reduces to 0 their number, so easing of the filtring
    {
    	threshold <- 0.02	
    	i <- enrichment$pathway[which(enrichment$padj <= 0.05 
                                  & enrichment$padj <= quantile(enrichment$padj, threshold) 
                                  & enrichment$age %in% temp$age[temp$sig == "sig"])] 
    	i <- unique(i)
    	threshold <- threshold + 0.01
    }
    
    
    
    #Matrix of NES over time for selcted pathwas
    i <- which(enrichment$pathway %in% i & enrichment$variable == variable)
    tmp <- dcast(enrichment[i,], pathway ~ age, value.var = "NES")
    rownames(tmp) <- tmp$pathway
    tmp <- tmp[, -1]
    tmp <- as.matrix(tmp[hclust(dist(tmp))$order,])
  } else
  {
    tmp <- dcast(enrichment[enrichment$variable == variable,], pathway ~ age, value.var = "NES")
    rownames(tmp) <- tmp$pathway
    tmp <- tmp[, -1]
    if (nrow(tmp) > 1)
    {
      tmp <- as.matrix(tmp[hclust(dist(tmp))$order,])
    } else
    {
      tmp <- as.matrix(tmp)
    }
  }
  
  #family of the selected pathways, l is created to add the family info the heatmapo
  temp <- Cluster_Reactome$cluster[match(rownames(tmp), Cluster_Reactome$pathways)]
  l <- list()
  for (i in 1:length(temp))
  {
    l <- c(l, list(list(y = 1, x = (i-1), value = temp[i], color = color2()[temp[i]], categories = rownames(tmp)[i])))
  }
  
  #Data preparaiton for the plot
  tmp <- melt(tmp)
  colnames(tmp) <- c("pathway", "age", "NES")
  
  
  g <- hchart(tmp, type = "heatmap", yAxis = 1, 
              hcaes(y = as.factor(age), x = pathway, value = NES),
              tooltip = list(pointFormat = HTML("<b>Age</b>: {point.age} <br> <b>Pathway</b> {point.pathway} <br> <b>NES</b> {point.NES:.2f}"), 
                             headerFormat = HTML(""))) %>%
    hc_chart(zoomType = "x") %>%
    hc_colorAxis(stops = list(list(0.1, "black"), 
                              list(0.3, "dodgerblue"), 
                              list(0.5, "white"), 
                              list(0.7, "tomato"), 
                              list(0.9, "black")),
                 min = -4, max = 4) %>%
    hc_legend(title = list(text = "NES")) %>%
    hc_add_series(data = l, type = "heatmap", showInLegend = F,
                  tooltip = list(pointFormat = HTML("<b>Family</b>: {point.value} <br> <b>Pathway</b> {point.categories}"), 
                                 headerFormat = HTML(""))) %>%
    hc_yAxis_multiples(
      list(height = "10%", categories = as.factor(c("Family")), labels  = list(enabled = FALSE), title = list(text = "Family")),
      list(showLastLabel = TRUE, opposite = TRUE, height = "90%", top = "10%", lineWidth = 1,
           categories = floor(unique(tmp$age)), 
           labels = list(step = 10, enabled = F),
           title = list(text = ""))
    ) %>%
    hc_xAxis(labels = list(enabled = F), tickWidth = 0, title = list(text = "Reactome pathways"))
  
  #adding of the percentage of the sig genes on the right side
  temp <- peak[peak$tissue == tissue & peak$variable == variable,]
  temp <- temp[, c("tissue", "age", "PercSigGene")]
  names(temp) <- c("x", "y", "value")
  
  
  
  g1 <- hchart(temp, type = "heatmap", 
               hcaes(x = x, y = as.character(y), value = value)) %>%
    hc_yAxis(height = "74%", top = "8.5%", bottom = "17.5%",
             opposite = T, labels = list(step = 10), categories = floor(temp$y),
             title = list(text = "Age (years)"), offset = -120, lineWidth = 0) %>%
    hc_xAxis(title = "", labels = list(enabled = F), tickWidth = 0, lineWidth = 0,
             width = '30%', right = "70%") %>%
    hc_colorAxis(stops = color_stops(10, rev(plasma(10))),
                 reversed = F) %>%
    hc_legend(title = list(text = "% altered <br> genes", style = list(transform = "rotate(90deg)", fontWeight = "bold")), 
              align = "right", useHTML = T, floating = T,
              layout = "vertical", verticalAlign = "top", x = -10) %>%
    hc_tooltip(pointFormat = HTML("<b>%alt genes</b>: {point.value:.2f} <br><b>Age:</b> {point.yf} y.o."), 
               headerFormat = HTML("")) 
  
  return(manipulateWidget::combineWidgets(g, g1, ncol = 2, colsize = c(3,1)))
}

wordcloud_family <- function(Cluster_Reactome_affiliation, family)
{
  
  #browser()
  tmp <- c(Cluster_Reactome_affiliation[[family]]$REACTOME.pathways,Cluster_Reactome_affiliation[[family]]$KEGG.GO.pathways) 
  tmp <- tmp[!is.na(tmp)]
  tmp <- gsub("REACTOME", "", tmp) 
  tmp <- gsub("GO", "", tmp)
  tmp <- gsub("KEGG", "", tmp)
  tmp <- gsub("_", " ", tmp)
  # remove uninformative words
  tmp <- gsub("PROCESS", "", tmp)
  tmp <- gsub("REGULATION", "", tmp)
  for (i in toupper(tm::stopwords("en")))
  {
    i <- paste0(" ", i, " ")
    tmp <- gsub(i, " ", tmp)
  }
  
  tmp <- tm::VCorpus(tm::VectorSource(tmp))
  tmp <- as.matrix(tm::TermDocumentMatrix(tmp))
  tmp <- sort(rowSums(tmp), decreasing = T)
  tmp <- data.frame(word = names(tmp), freq = tmp)
  tmp$color <- rev(colorRampPalette(c("black", color2()[family]))(2*max(tmp$freq)))[tmp$freq]
  
  if (max(tmp$freq) >= 3)
  {
    g <- hchart(tmp[tmp$freq >= 3,], "wordcloud", 
                hcaes(name = word, weight = freq, color = color), name = "Occurence") %>%
      hc_title(text = paste("Family", family))
  } else
  {
    g <- hchart(tmp[tmp$freq >= 2,], "wordcloud", 
                hcaes(name = word, weight = freq, color = color), name = "Occurence") %>%
      hc_title(text = paste("Family", family))
  }
  return(g)
}

Heatmap_family <- function(Cluster_Reactome)
{
  tmp <- data.frame(Family = 1:length(unique(Cluster_Reactome$cluster)))
  tmp$color <- color2()[tmp$Family]
  rownames(tmp) <- tmp$Family
  
  
  l <- list()
  for (i in 1:nrow(tmp))
  {
    l <- c(l, list(list(x = 1, y = (i-1), value = tmp$Family[i], 
                        color = tmp$color[i], categories = tmp$Family[i])))
  }
  
  
  ClickFunction <- JS("function(event) {Shiny.onInputChange('familyClicked', event.point.categories)  }")
  ClickFunction2 <- JS("function(event) {console.log(event)  }")
  
  g <- highchart() %>%
    hc_add_series(data = l, type = "heatmap", showInLegend = F,
                  tooltip = list(pointFormat = HTML("<b>Family</b>: {point.value}"), 
                                 headerFormat = HTML(""))) %>%
    hc_yAxis(title = list(text = "Family"), tickAmount = length(unique(Cluster_Reactome$cluster)),
             labels = list(formatter = JS("function(){return this.value+1;}"))) %>%
    hc_xAxis(labels = list(enabled = FALSE), tickWidth = 0)  %>%
    hc_plotOptions(series = list(events = list(click = ClickFunction)))
  return(g)
  
}

Heatmap_FisherTest_cellComposition <- function(pvalueFisherTest)
{
   
  pvalueFisherTest$value <- round(pvalueFisherTest$value,2)
  pvalueFisherTest$oddsRatio <- round(pvalueFisherTest$oddsRatio,2)
  
  pvalueFisherTest$cellType <- gsub("\\.", " ", pvalueFisherTest$cellType)
  names(pvalueFisherTest)[3] <- "value"
  #Labls authrisation: when too many modules, difficult read yaxis laebls
  legendLabels <- ifelse(length(unique(pvalueFisherTest$module)) <= 20, TRUE, FALSE)
  
  pvalueFisherTest$tooltip_info <- with(pvalueFisherTest, 
                                        paste("<b>Module:</b> ", module,
                                              "<br><b>Cell type:</b> ", cellType,
                                              "<br><b>Odds Ratio:</b> ", oddsRatio,
                                              "<br><b>-log<sub>10</sub>(p-value):</b> ", value))
  
  
  # First, we ensure that the `value` column represents -log10(p-value)
  # If it doesn't, you'll need to correct this.
  
  # Create a new column for labeling based on the conditions you mentioned
  pvalueFisherTest$label_value <- ifelse(pvalueFisherTest$oddsRatio > 1 & pvalueFisherTest$value > 1.3, pvalueFisherTest$oddsRatio, NA)
  
  # Now, you simply check for the existence of this value in the labeling function
  g <- hchart(pvalueFisherTest, type = "heatmap",
              hcaes(x = cellType, y = module, value = oddsRatio)) %>%
    hc_plotOptions(series = list(dataLabels = list(enabled = T, 
                                                   formatter = JS("function(){if(!isNaN(this.point.label_value)){return this.point.label_value;}}")))) %>%
    hc_legend(verticalAlign = "top", align = "left", layout = "vertical",
              title = list(text = "Odds Ratio"), useHTML = TRUE) %>%
    hc_colorAxis(stops = color_stops(9, RColorBrewer::brewer.pal(9, "Greens")),
                 min = 0, max = plyr::round_any(max(pvalueFisherTest$oddsRatio), 5, f = ceiling)) %>%
    hc_tooltip(pointFormat = "{point.tooltip_info}", useHTML = T) %>%
    hc_yAxis(labels = list(enabled = legendLabels), title = list(text = "Modules"),
             lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0) %>%
    hc_xAxis(title = list(text = "Cell type"), labels = list(autoRotation = F, rotation = 60)) # Add the rotation argument here
  
  
  
  
  return(g)
  # #Adding of heatmap to depict the modules
  # tmp <- data.frame(Module = levels(as.factor(pvalueFisherTest$module)))
  # tmp$color <- gsub("ME", "", tmp$Module)
  # 
  # #Since highcharter do not recognize all the color name, convert them to HEX odes
  # colorHEX <- data.frame(color=colors(), HEX = rgb(t(col2rgb(colors())), maxColorValue=255), stringsAsFactors=FALSE)
  # 
  # tmp$color <- colorHEX$HEX[match(tmp$color, colorHEX$color)]
  # rownames(tmp) <- tmp$Module
  # 
  # l <- list()
  # for (i in 1:nrow(tmp))
  # {
  #   l <- c(l, list(list(x = 1, y = (i-1), value = tmp$Module[i], 
  #                       color = tmp$color[i], categories = gsub(".*_", "", tmp$Module[i]))))
  # }
  # 
  # 
  # g1 <- highchart() %>%
  #   hc_add_series(data = l, type = "heatmap", showInLegend = F,
  #                 tooltip = list(pointFormat = HTML("<b>Module</b>: {point.categories}"), 
  #                                headerFormat = HTML(""))) %>%
  #   hc_yAxis(height = paste0(95 - as.numeric(round(predict(fit, data.frame(x = nrow(tmp))), 2)), "%"), #fit used here
  #            top = paste0(round(predict(fit, data.frame(x = nrow(tmp))), 2), "%"), bottom = "5%",
  #            title = list(text = "Modules"), lineWidth = 0,
  #            minorGridLineWidth = 0, gridLineWidth = 0,
  #            opposite = T,
  #            labels = list(enabled = FALSE)) %>%
  #   hc_xAxis(labels = list(enabled = FALSE), tickWidth = 0, lineWidth = 0) 
  # 
  # manipulateWidget::combineWidgets(g, 
  #                                  g1, ncol = 2, colsize = c(3,0.3))
  
  
  
}

Treemap_moduleFisherTest <- function(pvalueFisherTest, module)
{
  tmp <- pvalueFisherTest[pvalueFisherTest$module == module,] 
  names(tmp)[3] <- "value"
  #when oddsRatio are infinite, reduce them to 1000
  if (length(which(tmp$oddsRatio == Inf))!=0)
  {
    tmp$oddsRatio[which(tmp$oddsRatio == Inf)] <- 1000
  }
  
  l <- list()
  mycolor <- list()
  for (i in 1:length(unique(tmp$reference)))
  {
    l <- c(l, list(list(id = unique(tmp$reference)[i], name = unique(tmp$reference)[i])))
    mycolor[[i]] <- colorRampPalette(RColorBrewer::brewer.pal(9, c("Greens", "Oranges", "Purples")[i]))(max(tmp$oddsRatio))
  }
  
  for (i in 1:nrow(tmp))
  {
    l <- c(l, list(list(name = tmp$cellType[i], 
                        parent = tmp$reference[i],
                        value = tmp$value[i], 
                        oddsRatio = tmp$oddsRatio[i],
                        color = mycolor[[match(tmp$reference[i], unique(tmp$reference))]][ceiling(tmp$oddsRatio[i])]
    )))
  }
  #ClickFunction <- JS("function(event) {console.log(event)  }")
  g <- highchart()%>%
    hc_add_series(data = l, type = "treemap",
                  layoutAlgorithm =  'stripes', alternateStartingDirection = TRUE,
                  levels = list(list(level = 1, layoutAlgorithm = 'sliceAndDice',
                                     borderWidth = 2))) %>%
    #hc_plotOptions(series = list(events = list(click = ClickFunction)))  %>%
    hc_tooltip(formatter = JS("function(){return '<b>Cell type: </b>' + this.point.name + '<br><b>-log<sub>10</sub>(p): </b>' + this.point.value + '<br><b>Enrichment coefficient: </b>' + this.point.oddsRatio + '<br><b>Reference: </b>' + this.point.parent ;}"), useHTML = T)
  
  
  g1 <- list()
  for (i in 1:length(unique(tmp$reference)))
  {
    g1[[i]] <-  highchart() %>%
      hc_add_series(data  = list(), type = "treemap") %>%
      hc_colorAxis(minColor = RColorBrewer::brewer.pal(9, c("Greens", "Oranges", "Purples")[i])[1],
                   maxColor = RColorBrewer::brewer.pal(9, c("Greens", "Oranges", "Purples")[i])[9],
                   min = 0,
                   max = max(tmp$value)) %>%
      hc_legend(title = list(text = paste("Enrichment coefficient - ", unique(tmp$reference)[i], sep = "")), 
                verticalAlign = "top", align = "center", layout = "horizontal", useHTML = T)
  }
  
  
  manipulateWidget::combineWidgets(manipulateWidget::combineWidgets(list = g1, ncol = 3), g,
                                   nrow = 2, rowsize = c(1, 3))
}

Heatmap_LoessMEvsAge <- function(age)
{
  a <- data.frame()
  for (i in unique(age$module))
  {
    fit <- loess(value~age, data = age[age$module == i,])
    a <- rbind(a, 
               #data.frame(expression = scale(predict(fit, seq(20, 70, 0.5)), scale = F),
               data.frame(expression = scale(predict(fit, seq(20, 70, 0.5)), scale = T),
                          age = seq(20, 70, 0.5),
                          module = i))
  }
  
  fit <- dcast(a, module ~ age, value.var = "expression")
  fit <- fit[hclust(dist(fit[,-1]))$order,]
  rownames(fit) <- fit$module
  fit <- fit[,-1]
  fit <- melt(as.matrix(fit))
  names(fit) <- c("module", "age", "expression")
  
  g <- hchart(fit, type = "heatmap", hcaes(x = age, y = module, value = expression))    %>%
    hc_colorAxis(stops = list(list(0.3, "dodgerblue"), list(0.5, "white"), list(0.7, "tomato")),
                 reversed = F) %>%
    hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.age + ' y.o., <br><b>Module</b>: ' + this.point.module + ',<br><b>scaled GE</b>: ' + Highcharts.numberFormat(this.point.options.expression, 2);}")) %>% 
    hc_legend(title = list(text = "<center>Eigengene Expression<br>(Z-scores)</center>", style = list(fontSize = '16px')),align = "right", layout = "vertical", verticalAlign = "top",reversed = TRUE) %>%
    
    
     hc_xAxis(title = list(text = "Age (years)", style = list(fontSize = '18px')),
             labels = list(style = list(fontSize = "14px")),
             lineWidth = 0) %>%
    hc_yAxis(title = list(text = "", style = list(fontSize = '18px')),
             labels = list(enabled = FALSE),
             lineWidth = 0,
             minorGridLineWidth = 0, gridLineWidth = 0) 
  
  #Adding of heatmpa on the left to corroborate labels
  tmp <- data.frame(Module = levels(as.factor(fit$module)))
  tmp$color <- gsub("ME", "", tmp$Module)
  
  #Since highcharter do not recognize all the color name, convert them to HEX odes
  colorHEX <- data.frame(color=colors(), HEX = rgb(t(col2rgb(colors())), maxColorValue=255), stringsAsFactors=FALSE)
  
  tmp$color <- colorHEX$HEX[match(tmp$color, colorHEX$color)]
  rownames(tmp) <- tmp$Module
  
  l2 <- list()
  for (i in 1:nrow(tmp))
  {
    l2 <- c(l2, list(list(x = 1, y = (i-1), value = tmp$Module[i], 
                          color = tmp$color[i], categories = gsub(".*_", "", tmp$Module[i]))))
  }
  
  #The blank space at the top of the second heatmap is dependent of the tissue. A fixed one gives a bad alignment across heatmanp
  
  g1 <- highchart() %>%
    hc_add_series(data = l2, type = "heatmap", showInLegend = F,
                  tooltip = list(pointFormat = HTML("<b>Module</b>: {point.categories}"), 
                                 headerFormat = HTML(""))) %>%
    hc_tooltip(enabled = FALSE) %>%
    hc_yAxis(height = "92.5%", #fit used here
             top = "0%", bottom = "7.5%",
             title = list(text = "Modules"), lineWidth = 0,
             minorGridLineWidth = 0, gridLineWidth = 0,
             opposite = F,
             labels = list(enabled = TRUE),
             categories = rownames(tmp)) %>%
    hc_xAxis(labels = list(enabled = FALSE), tickWidth = 0, lineWidth = 0) 
  
  return(manipulateWidget::combineWidgets(g1, g, ncol = 2, colsize = c(0.9, 3)))
}

Scatter_MEvsAge <- function(MEexpressionData, module, colored, shaped, donorCondition, technicalCondition)
{ 
  tmp <- MEexpressionData[MEexpressionData$module == module,]
  tmp$jittered_age <- jitter_ages(tmp$age)
  
  ClickFunction2 <- JS("function(event) {console.log(event)  }")
  
  
  if (colored == "none" & shaped == "none")
  {
    fit <- loess(value ~ age, data = tmp)
    fit <- plyr::arrange(generics::augment(fit), age)
    
    #Not all color names are recognized by highcarther, conversion in HEX code
    colorHEX <- data.frame(color=colors(), HEX = rgb(t(col2rgb(colors())), maxColorValue=255), stringsAsFactors=FALSE)
    colorHEX <- colorHEX$HEX[match(gsub("ME", "", module), colorHEX$color)]
    
    g <- highchart() %>%
      hc_chart(type = 'scatter') %>%
      hc_add_series(data = list_parse2(tmp[, c("jittered_age", "value")]), showInLegend = FALSE) %>%
      # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), 
      #               type = "arearange", showInLegend = FALSE) %>%
      hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                    type = "line", showInLegend = FALSE) %>%
      hc_plotOptions(arearange = list(color = "gainsboro", opacity = 0.005),
                     line = list(lineWidth = 7, color = colorHEX),
                     scatter = list(marker = list(radius = 5, symbol = "line"), color = "grey"),
                     series = list(allowPointSelect = FALSE, 
                                   #events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  } else if (colored == "Sex" & shaped == "none")
  {
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      fit_male <- loess(value ~ age, data = tmp[tmp$sex == 1,])
      fit_male <- plyr::arrange(generics::augment(fit_male), age)
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      fit_female <- loess(value ~ age, data = tmp[tmp$sex == 2,])
      fit_female <- plyr::arrange(generics::augment(fit_female), age)
    }
    g <- highchart() %>%
      hc_chart(type = 'scatter') 
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1, c("jittered_age", "value")]), 
                      showInLegend = TRUE, color = "cornflowerblue", name = "Male") 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    { 
      g <- g %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2, c("jittered_age", "value")]), 
                      showInLegend = TRUE, color = "hotpink", name = "Female")
    }
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g %>%
        hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "cornflowerblue") 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      g <- g %>%
        hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "hotpink") 
    }
    g <- g %>%
      hc_plotOptions(line = list(lineWidth = 7),
                     scatter = list(marker = list(radius = 5, symbol = "circle", fillOpacity = 0.4)),
                     series = list(allowPointSelect = FALSE, 
                                   #events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  } else if (colored != "none" & colored != "Sex" & shaped == "none")
  {
    tmp$coloredBy <- technicalCondition[match(tmp$samp_id, technicalCondition$SAMPID), match(colored, colnames(technicalCondition))]
    if (class(tmp$coloredBy) == "character")
    {
      g <- highchart() %>%
        hc_chart(type = 'scatter') 
      
      mysymbol <- c('circle', 'square','diamond', 'triangle', 'triangle-down')
      mycolor <- RColorBrewer::brewer.pal(8, "Dark2")
      for (i in 1:length(unique(tmp$coloredBy)[!is.na(unique(tmp$coloredBy))]))
      {
        g <- g %>%
          hc_add_series(data = list_parse2(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i], c("jittered_age", "value")]),
                        showInLegend = TRUE,
                        marker = list(radius = 5, symbol = mysymbol[i]), 
                        color = mycolor[i], name = unique(tmp$coloredBy)[i])
        # if (nrow(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i],]) >= 10)
        # {
        #   fit <- loess(value ~ age, data = tmp[tmp$coloredBy == unique(tmp$coloredBy)[i],])
        #   fit <- plyr::arrange(generics::augment(fit), age)
        #   g <- g %>% hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
        #                            type = "line", showInLegend = FALSE,
        #                            lineWidth = 7, color = mycolor[i])
        # }
        
      }
      g <- g 
    } else if (class(tmp$coloredBy) == "numeric")
    {
      mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(length(seq(floor(round(min(tmp$coloredBy), 1)), ceiling(round(max(tmp$coloredBy), 1)), by = 0.1)))
      tmp$colorCondition <- mycolor[match(tmp$coloredBy[order(tmp$coloredBy)], round(seq(min(tmp$coloredBy), max(tmp$coloredBy), by = 0.1), 1))]
      
      fit <- loess(value ~ age, data = tmp)
      fit <- plyr::arrange(generics::augment(fit), age)
      tmp$age <- tmp$jittered_age
      
      g <- hchart(tmp, hcaes(x = age, y = value, color = colorCondition), type = "scatter") %>%
        hc_colorAxis(stops = color_stops(9, colors = RColorBrewer::brewer.pal(9, "Reds")),
                     min = min(tmp$coloredBy), max = max(tmp$coloredBy)) %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)),
        #               type = "arearange", showInLegend = FALSE,
        #               color = "gainsboro", opacity = 0.005) %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
        #               type = "line", showInLegend = FALSE,
        #               lineWidth = 7, color = "black") %>%
        hc_legend(title = list(text = colored))
      
      
    }
  } else if (colored == "none" & shaped != "none")
  {
    tmp$condition <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), match(shaped, colnames(donorCondition))]
    
    
    #Fitting of a regression per condition if the number of positive/negative case is higher than 10
    if (length(which(tmp$condition %in% c(0,1))) == 0)
    {
      fit <- loess(value ~ age, data = tmp)
      fit <- plyr::arrange(generics::augment(fit), age)
      
      g <- highchart() %>%
        hc_chart(type = 'scatter') %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)),
        #               type = "arearange", showInLegend = FALSE,
        #               color = "gainsboro", opacity = 0.005) %>%
        hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "black") %>%
        hc_plotOptions(series = list(allowPointSelect = FALSE, 
                                     #events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
      
    } else if (min(table(tmp$condition)[names(table(tmp$condition)) %in% c(0,1)]) >= 10)
    {
      fit_positive <- loess(value ~ age, data = tmp[tmp$condition == 1,])
      fit_positive <- plyr::arrange(generics::augment(fit_positive), age)
      
      fit_negative <- loess(value ~ age, data = tmp[tmp$condition == 0,])
      fit_negative <- plyr::arrange(generics::augment(fit_negative), age)
      
      g <- highchart() %>%
        hc_chart(type = 'scatter') %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 1, c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "triangle"), color = "rgb(60,179,113)", name = "Positive") %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 0, c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "cross"), color = "rgb(255,99,71)", name = "Negative") %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        hc_add_series(data = list_parse2(data.frame(fit_positive$age, fit_positive$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "rgb(60,179,113)") %>%
        hc_add_series(data = list_parse2(data.frame(fit_negative$age, fit_negative$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "rgb(255,99,71)") %>%
        hc_plotOptions(series = list(allowPointSelect = FALSE, 
                                     #events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
      
      
    } else
    {
      
      fit <- loess(value ~ age, data = tmp)
      fit <- plyr::arrange(generics::augment(fit), age)
      
      g <- highchart() %>%
        hc_chart(type = 'scatter') %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 1, c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "triangle"), color = "rgb(60,179,113)", name = "Positive") %>%
        hc_add_series(data = list_parse2(tmp[tmp$condition == 0, c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 7, symbol = "cross"), color = "rgb(255,99,71)", name = "Negative") %>%
        hc_add_series(data = list_parse2(tmp[!(tmp$condition %in% c(0, 1)), c("jittered_age", "value")]), showInLegend = TRUE,
                      marker = list(radius = 5, symbol = "circle"), color = "grey", name = "Unknown") %>%
        # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), 
        #               type = "arearange", showInLegend = FALSE,
        #               color = "gainsboro", opacity = 0.005) %>%
        hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
                      type = "line", showInLegend = FALSE,
                      lineWidth = 7, color = "black") %>%
        hc_plotOptions(series = list(allowPointSelect = FALSE, 
                                     #events = list(click = ClickFunction), 
                                     marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
    }
    
  } else if (colored == "Sex" & shaped != "none")
  {
    
    tmp$condition <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), match(shaped, colnames(donorCondition))]
    
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      fit_male <- loess(value ~ age, data = tmp[tmp$sex == 1,])
      fit_male <- plyr::arrange(generics::augment(fit_male), age)
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      fit_female <- loess(value ~ age, data = tmp[tmp$sex == 2,])
      fit_female <- plyr::arrange(generics::augment(fit_female), age)
    }
    g <- highchart() %>%
      hc_chart(type = 'scatter')
    #Male and conditions
    if (nrow(tmp[tmp$sex == 1,])!=0)
    { 
      g <- g %>%  
        
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & tmp$condition == 1, c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Pos",
                      marker = list(radius = 7, symbol = "triangle")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & tmp$condition == 0, c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Neg",
                      marker = list(radius = 7, symbol = "cross")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 1 & !(tmp$condition %in% c(0, 1)), c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "cornflowerblue", name = "Male - Unknown",
                      marker = list(radius = 5, symbol = "circle")) 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    { 
      g <- g %>%
        #Female and conditions
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & tmp$condition == 1, c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Pos",
                      marker = list(radius = 7, symbol = "triangle")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & tmp$condition == 0, c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Neg",
                      marker = list(radius = 7, symbol = "cross")) %>%
        hc_add_series(data = list_parse2(tmp[tmp$sex == 2 & !(tmp$condition %in% c(0, 1)), c("jittered_age", "value")]), 
                      showInLegend = FALSE, color = "hotpink", name = "Female - Unknown",
                      marker = list(radius = 5, symbol = "circle")) 
    }
    if (nrow(tmp[tmp$sex == 1,])!=0)
    {
      g <- g %>%
        #Loess fitting for both sex
        hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "cornflowerblue", lineWidth = 7) 
    }
    if (nrow(tmp[tmp$sex == 2,])!=0)
    {
      g <- g%>%
        hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), 
                      type = "line", showInLegend = FALSE, color = "hotpink", lineWidth = 7) 
    }
    g <- g%>%
      #Fake series to create legend
      hc_add_series(type = "scatter", name = "Positive", color=  "grey", 
                    marker = list(symbol = "triangle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Negative", color=  "grey", 
                    marker = list(symbol = "cross", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Unknown", color=  "grey", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Male", color=  "cornflowerblue", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      hc_add_series(type = "scatter", name = "Female", color=  "hotpink", 
                    marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) %>%
      #Options for point click and selection
      hc_plotOptions(series = list(allowPointSelect = FALSE, 
                                   #events = list(click = ClickFunction), 
                                   marker = list(states = list(select = list(fillColor = "lightgrey", lineWidth = 1, radius = 6, lineColor = "black"))))) 
  } else if (colored != "none" & colored != "Sex" & shaped != "none")
  {
    tmp$coloredBy <- technicalCondition[match(tmp$samp_id, technicalCondition$SAMPID), match(colored, colnames(technicalCondition))]
    
    tmp$condition <- donorCondition[match(tmp$samp_id, donorCondition$SAMPID), match(shaped, colnames(donorCondition))]
    if (class(tmp$coloredBy) == "character")
    {
      g <- highchart() %>%
        hc_chart(type = 'scatter') 
      
      mysymbol <- c('circle', 'square','diamond', 'triangle', 'triangle-down')
      mycolor <- RColorBrewer::brewer.pal(8, "Dark2")
      for (i in 1:length(unique(tmp$coloredBy)[!is.na(unique(tmp$coloredBy))]))
      {
        g <- g %>%
          hc_add_series(data = list_parse2(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i] & tmp$condition == 1, 
                                               c("jittered_age", "value")]),
                        marker = list(radius = 5, symbol = "triangle"), 
                        color = mycolor[i], showInLegend = F) %>%
          hc_add_series(data = list_parse2(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i] & tmp$condition == 0, 
                                               c("jittered_age", "value")]),
                        marker = list(radius = 5, symbol = "cross"), 
                        color = mycolor[i], showInLegend = F) %>%
          hc_add_series(data = list_parse2(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i] & !(tmp$condition %in% c(0, 1)), 
                                               c("jittered_age", "value")]),
                        marker = list(radius = 5, symbol = "circle"), 
                        color = mycolor[i], showInLegend = F)
        # if (nrow(tmp[tmp$coloredBy == unique(tmp$coloredBy)[i],]) >= 10)
        # {
        #   fit <- loess(value ~ age, data = tmp[tmp$coloredBy == unique(tmp$coloredBy)[i],])
        #   fit <- plyr::arrange(generics::augment(fit), age)
        #   g <- g %>% hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), 
        #                            type = "line", showInLegend = FALSE,
        #                            lineWidth = 7, color = mycolor[i])
        # }
        
      }
      #Fake series for the legend
      g <- g %>%
        hc_add_series(type = "scatter", name = "Positive", color=  "grey", 
                      marker = list(symbol = "triangle", radius = 8), showInLegend = TRUE) %>%
        hc_add_series(type = "scatter", name = "Negative", color=  "grey", 
                      marker = list(symbol = "cross", radius = 8), showInLegend = TRUE) %>%
        hc_add_series(type = "scatter", name = "Unknown", color=  "grey", 
                      marker = list(symbol = "circle", radius = 8), showInLegend = TRUE) 
      for (i in 1:length(unique(tmp$coloredBy)))
      {
        g <- g %>%
          hc_add_series(type = "scatter", name = unique(tmp$coloredBy)[i], color=  mycolor[i], 
                        marker = list(symbol = "circle", radius = 8), showInLegend = TRUE)
      }
      
    } else if (class(tmp$coloredBy) == "numeric")
    {
      mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(length(seq(floor(round(min(tmp$coloredBy), 1)), ceiling(round(max(tmp$coloredBy), 1)), by = 0.1)))
      tmp$colorCondition <- mycolor[match(tmp$coloredBy[order(tmp$coloredBy)], round(seq(min(tmp$coloredBy), max(tmp$coloredBy), by = 0.1), 1))]
      
      tmp$age <- tmp$jittered_age
      
      g <- highchart() %>%
        hc_chart(type = "scatter") %>%
        hc_colorAxis(stops = color_stops(9, colors = RColorBrewer::brewer.pal(9, "Reds")),
                     min = min(tmp$coloredBy), max = max(tmp$coloredBy)) %>%
        hc_legend(title = list(text = colored))
      
      tmp$condition[which(!(tmp$condition %in% c(0,1)))] <- 99
      
      for (j in unique(tmp$condition))
      {
        l <- list()
        for (i in 1:nrow(tmp[tmp$condition == j,]))
        {
          l<- c(l, list(list(y = tmp[tmp$condition == j,]$value[i], 
                             x = tmp[tmp$condition == j,]$age[i], 
                             color = tmp[tmp$condition == j,]$colorCondition[i])))
        }
        
        g <- g %>%
          hc_add_series(data = l, color = "grey",
                        marker = list(symbol = ifelse(j == 1, "triangle", ifelse(j == 0, "cross", "circle")),
                                      radius = 7),
                        name = ifelse(j == 1, "Positive", ifelse(j == 0, "Negative", "Unknown")))
      }
    }
  }
  
  g <- g %>%
    hc_xAxis(title = list(text = "Age (years)"), min = 20, max = 70) %>%
    hc_yAxis(title = list(text = "Eigengene Expression (logCPM)")) %>%
    hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>Expression</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}"), useHTML = TRUE) %>%
    hc_title(text = module) 
  
  return(g)
}  

Heatmap_moduleEnrichment <- function(enrich, Cluster_Reactome, module = F, tissue)
{
  
  enrichv0 <- enrich
  
  i <- hclust(dist(enrich))$order
  j <- hclust(dist(t(enrich)))$order
  enrich <- enrich[i,j]
  
  #family of the selected pathways, l is created to add the family info the heatmapo
  temp <- Cluster_Reactome$cluster[match(colnames(enrich), Cluster_Reactome$pathways)]
  l <- list()
  for (i in 1:length(temp))
  {
    l <- c(l, list(list(y = 1, x = (i-1), value = as.character(temp[i]), color = color2()[temp[i]], categories = colnames(enrich)[i])))
  }
  
  #Data preparaiton for the plot
  enrich <- melt(enrich)
  
  colnames(enrich) <- c("module", "pathway", "-log10(p)")
  enrich$module <- gsub(".*_", "" ,enrich$module)
  
  if (module != F)
  {
    enrich <- enrich[enrich$module == module, ]
  }
  
  g <- hchart(enrich, type = "heatmap", yAxis = 1, name = "fisher",
              hcaes(y = module, x = pathway, value = `-log10(p)`),
              tooltip = list(pointFormat = HTML("<b>Module</b>: {point.module} <br> <b>Pathway</b> {point.pathway} <br> <b>-log10(p)</b> {point.-log10(p):.2f}"), 
                             headerFormat = HTML(""))) %>%
    hc_chart(zoomType = "x") %>%
    hc_colorAxis(stops = color_stops(9, RColorBrewer::brewer.pal(9, "Oranges"))) %>%
    hc_plotOptions(heatmap = list(dataLabels = list(enabled = T, 
                                                    formatter = JS("function(){if(this.point.value >= 2 && this.point.series.name == 'fisher'){return '**';} else if (this.point.value >= 1.3 && this.point.series.name == 'fisher'){return '*';}}")))) %>%
    hc_legend(title = list(text = "-log<sub>10</sub>(p)"), useHTML = TRUE,
              verticalAlign = "top", align = "left", layout = "vertical") %>%
    hc_add_series(data = l, type = "heatmap", showInLegend = F,
                  tooltip = list(pointFormat = HTML("<b>Family</b>: {point.value} <br> <b>Pathway</b> {point.categories}"), 
                                 headerFormat = HTML(""))) %>%
    hc_yAxis_multiples(
      list(height = "10%", categories = as.factor(c("Family")), 
           lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0,
           labels  = list(enabled = FALSE), title = list(text = "Family")),
      list(showLastLabel = TRUE, opposite = TRUE, height = "90%", top = "10%", 
           lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0,
           title = list(text = ""),
           labels  = list(enabled = FALSE),
           categories = levels(as.factor(enrich$module)))
    ) %>%
    hc_xAxis(labels = list(enabled = F), 
             tickWidth = 0, lineWidth = 0,
             title = list(text = "Reactome pathways"))
  
  #Adding of heatmpa on the right to replace labels
  tmp <- data.frame(Module = rownames(enrichv0))
  tmp$color <- gsub("ME", "", tmp$Module)
  
  #Since highcharter do not recognize all the color name, convert them to HEX odes
  colorHEX <- data.frame(color=colors(), HEX = rgb(t(col2rgb(colors())), maxColorValue=255), stringsAsFactors=FALSE)
  
  tmp$color <- colorHEX$HEX[match(tmp$color, colorHEX$color)]
  rownames(tmp) <- tmp$Module
  
  l2 <- list()
  for (i in 1:nrow(tmp))
  {
    l2 <- c(l2, list(list(x = 1, y = (i-1), value = tmp$Module[i], 
                          color = tmp$color[i], categories = gsub(".*_", "", tmp$Module[i]))))
  }
  
  #The blank space at the top of the second heatmap is dependent of the tissue. A fixed one gives a bad alignment across heatmanp
  fit <- data.frame(tissue = c("Brain - Cortex", "Muscle - Skeletal", "Whole Blood", "Heart - Left Ventricle"),
                    height = c(5, 3.5, 10, 10))
  g1 <- highchart() %>%
    hc_add_series(data = l2, type = "heatmap", showInLegend = F,
                  tooltip = list(pointFormat = HTML("<b>Module</b>: {point.categories}"), 
                                 headerFormat = HTML(""))) %>%
    hc_yAxis(height = paste0(96.5 - fit$height[fit$tissue == tissue], "%"), #fit used here
             top = paste0(fit$height[fit$tissue == tissue], "%"), bottom = "3.5%",
             title = list(text = "Modules"), lineWidth = 0,
             minorGridLineWidth = 0, gridLineWidth = 0,
             opposite = T,
             labels = list(enabled = FALSE)) %>%
    hc_xAxis(labels = list(enabled = FALSE), tickWidth = 0, lineWidth = 0) 
  
  return(manipulateWidget::combineWidgets(g, g1, ncol = 2, colsize = c(3,0.3)))
  
  
}

Heatmap_diseaseEnrichment <- function(diseaseEnrichment, disease)
{
  #browser()
  tmp <- diseaseEnrichment
  tmp$log10p <- -log10(tmp$padj)
  tmp <- tmp[tmp$Disease %in% disease,]
  tmp <- dcast(data = tmp, Module ~ Disease, value.var = "log10p")
  
  tmp[is.na(tmp)] <- 0
  
  if (length(disease) > 1)
  {
    rownames(tmp) <- tmp$Module
    tmp <- tmp[, -1]
    i <- hclust(dist(tmp))$order
    j <- hclust(dist(t(tmp)))$order
    tmp <- tmp[i,j]
  } else
  {
    tmp <- matrix(tmp[,2], dimnames = list(tmp$Module, colnames(tmp)[2]))
  }
  tmp <- melt(as.matrix(tmp))
  names(tmp) <- c("Modules", "Disease", "log10(p)")
  
  g <- hchart(tmp, type = "heatmap",
              hcaes(y = Modules, x = Disease, value = `log10(p)`),
              tooltip = list(pointFormat = HTML("<b>Module</b>: {point.Modules} <br> <b>Disease</b> {point.Disease} <br> <b>-log10(p)</b> {point.log10(p):.2f}"), 
                             headerFormat = HTML(""))) %>%
    hc_chart(zoomType = "x") %>%
    hc_colorAxis(stops = color_stops(9, RColorBrewer::brewer.pal(9, "Oranges"))) %>%
    hc_plotOptions(heatmap = list(dataLabels = list(enabled = T, 
                                                    formatter = JS("function(){if(this.point.value >= 2){return '**';} else if (this.point.value >= 1.3){return '*';}}")))) %>%
    hc_legend(title = list(text = "-log<sub>10</sub>(p)"), useHTML = TRUE,
              verticalAlign = "top", align = "left", layout = "vertical",
              reversed = T) %>%
    hc_yAxis_multiples(showLastLabel = TRUE, 
                       labels = list(enabled = T),
                       lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0,
                       title = list(text = "Modules"),
                       categories = levels(as.factor(tmp$Modules))) %>%
    hc_xAxis(tickWidth = 0, lineWidth = 0,
             labels = list(enabled = ifelse(length(disease) > 1, TRUE, FALSE)),
             title = list(text = "Diseases"))
  
  return(g)
}

Line_FishertestEnrichmentManualvsAge <- function(p, gene, threshold)
{
  fisher <- data.frame()
  for (i in 1:ncol(p))
  {
    fisher <- rbind(fisher, 
                    data.frame(age = as.numeric(colnames(p)[i]),
                               pvalue = fisher.test(
                                 matrix(c(length(which(p[match(gene, rownames(p)),i] <= threshold)), #gene from geneset DEG
                                          length(which(p[match(gene, rownames(p)),i] > threshold)), #genes from geneset not DEG
                                          length(which(p[,i] <= threshold)), #genes from all DEG
                                          length(which(p[,i] > threshold))), #genes from all not DEG
                                        nrow = 2), alternative = "greater")$p.value)) 
  }
  
  ClickFunction <- JS("function(event) {Shiny.onInputChange('ageContingencyTableClicked', event.point.x)  }")
  fisher$pvalue <- -log10(fisher$pvalue)
  fisher$color <- circlize::colorRamp2(breaks = c(0, 1.3, 2), colors = c("#dfe59a", "#9dca94", "#346875"))(fisher$pvalue)
  
  g <- hchart(fisher, type = "scatter",
              hcaes(x= age, y= pvalue, color = color)) %>%
    hc_xAxis(title = list(text = "Age (years)"), min = 20, max = 70) %>%
    hc_yAxis(title = list(text = "Fisher tests -log<sub>10</sub>(p-value)", 
                          useHTML = T),
             plotLines = list(list(value = 1.3, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.01")))) %>%
    hc_plotOptions(series = list(events = list(click = ClickFunction),
                                 lineWidth = 4, color = "#4D4D4D",
                                 marker = list(radius = 6, enabled = T))) %>%
    hc_tooltip(formatter = JS("function(){
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>-log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);
                                }"), useHTML = TRUE) %>% 
    hc_legend(layout = "horizontal", reversed = T,
              align = "center", title = list(text = "-log<sub>10</sub>(p-value)"), useHTML = TRUE) %>%
    hc_colorAxis(stops = list(list(0, "#dfe59a"), list(0.4, "#9dca94"), list(1.3, "#346875")),
                 min = 0, max = max(fisher$pvalue))
  return(g)
  
  
  
  
  
}

Line_signficanceAlterationsvsAge <- function(gene, pvalueData, geneData, allAges, abline, tissue, coloredBy)
{ 
 
  
  #GE info
  GE <- geneData
  GE <- GE[GE$tissue == tissue,]
  GE$jittered_age <- jitter_ages(GE$age)
  statsAllAges_text <- paste0("Overall changes: p-value = ", round(allAges$pvalue,4), " | t-statistic = ", round(allAges$tvalue,2), " | logFC/year = ", round(abline$a,4))
    
    
   
 
  if (coloredBy == "Age")
  {
    
    fit_abline <- data.frame(x=seq(20,70,0.1),
                             y=abline$a * (seq(20,70,0.1) - mean(GE$age))+ abline$b)
    
    
    fit <- loess(expression ~ age, data = GE)
    fit <- plyr::arrange(generics::augment(fit), age)
    

    
    #Pvalue info
    tmp <- melt(as.matrix(pvalueData))
    colnames(tmp) <- c("tissue", "Age", "p-value")
    tmp$log10p <- -log10(tmp$`p-value`)
    
    g <-   highchart() %>% 
      hc_yAxis_multiples(
        list(lineWidth = 2, opposite = TRUE, 
             title = list(useHTML = TRUE, text = "-log<sub>10</sub>(p-value)", 
                          style = list(color = "#9A839A")), 
             lineColor = "#9A839A", labels = list(style = list(color = "#9A839A")), height = "50%", top = "50%",
             plotLines = list(list(value = 1.3, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 3, dashStyle = "dot", label = list(text = "p = 0.01")))),
        list(showLastLabel = FALSE, opposite = FALSE, 
             title = list(text = "GE (logCPM)", style = list(color = "yellowgreen")), 
             height = "50%",
             lineColor = "yellowgreen", lineWidth = 2,
             labels = list(style = list(color = "yellowgreen")))
      )  %>% 
      # hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted - 2*fit$.se.fit, fit$.fitted + 2*fit$.se.fit)), yAxis = 1,
      #               type = "arearange", showInLegend = FALSE, color = "lightgrey") %>%
      hc_add_series(data = list_parse2(data.frame(fit$age, fit$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "black") %>% 
      hc_add_series(data = list_parse2(tmp[, c("Age", "log10p")]), type = "line", 
                    showInLegend = FALSE, name = "log10p",
                    color = "#9A839A", marker = list(radius = 6, enabled = T, color = "black"), lineWidth = 5) %>%
      hc_xAxis(title = list(text = "Age (years)"), lineColor = "grey", min = 20, max = 70) %>%
      hc_tooltip(formatter = JS("function(){
                                if(this.point.series.name == 'log10p'){ 
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);}
                                else{return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}
                                }"), useHTML = TRUE)%>% 
      hc_add_series(data = list_parse2(GE[, c("jittered_age", "expression")]), type = "scatter", showInLegend = FALSE,
                    yAxis = 1, color = "yellowgreen", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(
        data = list_parse2(data.frame(fit_abline$x, fit_abline$y)),
        type = "line",
        showInLegend = FALSE,
        yAxis = 1,
        name = "Entire Age Range",
        color = "#E9724C",  # Set the line color to blue
        lineWidth = 10,    # Set the line width
        dashStyle = "dot", # Set the line style to dotted
        marker = list(radius = 0, symbol = "round")
      )  %>% 
      hc_subtitle(text = statsAllAges_text, style = list(color = "#E9724C") ) 
    
  } else if (coloredBy == "Sex") {
     
      statsAllAges_text <- paste0("Overall changes: p-value = ", round(allAges$pvalue,4), " | t-statistic = ", round(allAges$tvalue,2), " | logFC = ", round(abline$a,4)) #Change the logFC/year to logFC
 
    
    
    fit_male <- loess(expression ~ age, data = GE[GE$sex == 1,]) 
    fit_male <- plyr::arrange(generics::augment(fit_male), age)
    
    fit_female <- loess(expression ~ age, data = GE[GE$sex == 2,]) 
    fit_female <- plyr::arrange(generics::augment(fit_female), age) 
     
    x <- mean(GE$age) 
    
     stats_male <- data.frame(x=x,
                              y=abline$a * (1 - mean(GE$sex)) + abline$b)
     
     
     stats_female <- data.frame(x=x,
                              y=abline$a * (2 - mean(GE$sex)) + abline$b)
    
    
    
    #Pvalue info
    tmp <- melt(as.matrix(pvalueData))
    colnames(tmp) <- c("tissue", "Age", "p-value")
    tmp$log10p <- -log10(tmp$`p-value`)
    
    g <-   highchart() %>% 
      hc_yAxis_multiples(
        list(lineWidth = 2, opposite = TRUE, 
             title = list(useHTML = TRUE, text = "log<sub>10</sub>(p-value)", 
                          style = list(color = "#9A839A")), lineColor = "#9A839A", 
             labels = list(style = list(color = "#9A839A")), 
             height = "50%", top = "50%",
             plotLines = list(list(value = 1.3, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.01")))),
        list(showLastLabel = FALSE, opposite = FALSE, 
             height = "50%", 
             title = list(text = "GE (logCPM)", style = list(color = "grey")), 
             lineColor = "grey", lineWidth = 2,
             labels = list(style = list(color = "grey")))
      ) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 1, c("jittered_age", "expression")]), type = "scatter", name = "Male",
                    yAxis = 1, color = "cornflowerblue", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 2, c("jittered_age", "expression")]), type = "scatter", name = "Female",
                    yAxis = 1, color = "hotpink", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "cornflowerblue") %>% 
      hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "hotpink") %>% 
      hc_add_series(data = list_parse2(tmp[, c("Age", "log10p")]), 
                    type = "line", showInLegend = FALSE, name = "log10p",
                    color = "#9A839A", marker = list(radius = 6, enabled = T, color = "black"), lineWidth = 5) %>%
      hc_xAxis(title = list(text = "Age (years)"), lineColor = "grey", min = 20, max = 70) %>%
      hc_tooltip(formatter = JS("function(){
                                if(this.point.series.name == 'log10p'){ 
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);}
                                else{return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}
                                }"), useHTML = TRUE) %>% 
      hc_add_series(
        data = list_parse2(data.frame(stats_male$x, stats_male$y)),
        type = "scatter",
        showInLegend = FALSE,
        yAxis = 1,
        name = "Entire Age Range (male)",
        color = "blue", 
        marker = list(radius = 10, symbol = "cross")
      )%>% 
      hc_add_series(
        data = list_parse2(data.frame(stats_female$x, stats_female$y)),
        type = "scatter",
        showInLegend = FALSE,
        yAxis = 1,
        name = "Entire Age Range (female)",
        color = "mediumvioletred",  # Set the line color to blue 
        marker = list(radius = 10, symbol = "cross")
      ) %>%  
      hc_subtitle(text = statsAllAges_text, style = list(color = "black") ) 
    
  } else {
    
    statsAllAges_text <- paste0("Overall changes: p-value = ", round(allAges$pvalue,4), " | t-statistic = ", round(allAges$tvalue,2), " | logFC/year = ", round(abline$a,4)) #Change the logFC/year to logFC
    
    
    
    fit_male <- loess(expression ~ age, data = GE[GE$sex == 1,]) 
    fit_male <- plyr::arrange(generics::augment(fit_male), age)
    
    fit_female <- loess(expression ~ age, data = GE[GE$sex == 2,]) 
    fit_female <- plyr::arrange(generics::augment(fit_female), age)
    
    
    fit_abline_interaction <- data.frame(x=seq(20,70,0.1),
                                         y=abline$a * (seq(20,70,0.1) - mean(GE$age)*mean(GE$sex))+ abline$b)
    
    
     
    
    #Pvalue info
    tmp <- melt(as.matrix(pvalueData))
    colnames(tmp) <- c("tissue", "Age", "p-value")
    tmp$log10p <- -log10(tmp$`p-value`)
    
    g <-   highchart() %>% 
      hc_yAxis_multiples(
        list(lineWidth = 2, opposite = TRUE, 
             title = list(useHTML = TRUE, text = "log<sub>10</sub>(p-value)", 
                          style = list(color = "#9A839A")), lineColor = "#9A839A", 
             labels = list(style = list(color = "#9A839A")), 
             height = "50%", top = "50%",
             plotLines = list(list(value = 1.3, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.05")),
                              list(value = 2, color = "black", width = 1, dashStyle = "dot", label = list(text = "p = 0.01")))),
        list(showLastLabel = FALSE, opposite = FALSE, 
             height = "50%", 
             title = list(text = "GE (logCPM)", style = list(color = "grey")), 
             lineColor = "grey", lineWidth = 2,
             labels = list(style = list(color = "grey")))
      ) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 1, c("jittered_age", "expression")]), type = "scatter", name = "Male",
                    yAxis = 1, color = "cornflowerblue", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(GE[GE$sex == 2, c("jittered_age", "expression")]), type = "scatter", name = "Female",
                    yAxis = 1, color = "hotpink", marker = list(radius = 4, symbol = "round")) %>% 
      hc_add_series(data = list_parse2(data.frame(fit_male$age, fit_male$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "cornflowerblue") %>% 
      hc_add_series(data = list_parse2(data.frame(fit_female$age, fit_female$.fitted)), yAxis = 1, lineWidth = 5,
                    type = "line", showInLegend = FALSE, color = "hotpink") %>% 
      hc_add_series(data = list_parse2(tmp[, c("Age", "log10p")]), 
                    type = "line", showInLegend = FALSE, name = "log10p",
                    color = "#9A839A", marker = list(radius = 6, enabled = T, color = "black"), lineWidth = 5) %>%
      hc_xAxis(title = list(text = "Age (years)"), lineColor = "grey", min = 20, max = 70) %>%
      hc_tooltip(formatter = JS("function(){
                                if(this.point.series.name == 'log10p'){ 
                                return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>log<sub>10</sub>(p-value)</b>: ' + Highcharts.numberFormat(this.point.y, 2);}
                                else{return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>GE</b>: ' + Highcharts.numberFormat(this.point.y, 2) + ' logCPM';}
                                }"), useHTML = TRUE)%>% 
      # hc_add_series(
      #   data = list_parse2(data.frame(fit_abline_interaction$x, fit_abline_interaction$y)),
      #   type = "line",
      #   showInLegend = FALSE,
      #   yAxis = 1,
      #   name = "Entire Age Range",
      #   color = "black",  # Set the line color to blue
      #   lineWidth = 10,    # Set the line width
      #   dashStyle = "dot", # Set the line style to dotted
      #   marker = list(radius = 0, symbol = "round")
      # ) %>%  
      hc_subtitle(text = statsAllAges_text, style = list(color = "black") ) 
    
    
  }
   
  return(g)
}

Heatmap_signficanceAlterationsvsAge <- function(p)
{
  tmp <- p
  tmp[is.na(tmp)] <- 0
  tmp[tmp < 0] <- 0
  #tmp <- -log10(tmp)
  if (nrow(tmp) > 1)
  {
    i <- (hclust(dist(tmp)))
    
    tmp <- as.matrix(tmp[i$order, ])
  } else
  {
    tmp <- as.matrix(tmp)
  }
  tmp <- melt(tmp)
  colnames(tmp) <- c("tissue", "age", "-log10(p)")
  #When the heatmap has only one row, the tissue name does not appear (just its first letter)
  yAxisTitle <- ifelse(length(unique(tmp$tissue)) == 1, as.character(tmp$tissue[1]), "Tissue")
  
  g <- hchart(tmp, type = "heatmap", boost = list(useGPUTranslations = T),
              hcaes(x = age, y = tissue, value = `-log10(p)`)) %>% 
    hc_colorAxis(stops = color_stops(20, rev(magma(20))),
                 min = 0,
                 reversed = F) %>%
    hc_legend(layout = "vertical", verticalAlign = "top",
              align = "right", title = list(text = "-log<sub>10</sub>(p)", style = list(fontSize = '16px')), useHTML = TRUE) %>%
    hc_xAxis(title = list(text = "Age (years)", style = list(fontSize = '18px')),
             labels = list(style = list(fontSize = "14px")),
             lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0) %>%
    hc_yAxis(title = list(text = yAxisTitle, style = list(fontSize = '18px')), reversed = TRUE,
             lineWidth = 0, minorGridLineWidth = 0, gridLineWidth = 0) %>%
    hc_tooltip(formatter = JS("function(){return '<b>Age</b>: ' + this.point.x + ' y.o., <br><b>Tissue</b>: ' +  this.series.yAxis.categories[this.point.y] + ',<br><b> -log<sub>10</sub>(p)</b>: ' + Highcharts.numberFormat(this.point.value, 2);}"), useHTML = T) 
  
  return(g)
}

Histogram_DistributionSample <- function(tissue, gene)
{
  a <- gene
  tmp <- a[a$tissue == tissue,]
  temp <- table(tmp$age, tmp$sex)
  
  
  g <- highchart() %>%
    hc_chart(type = "column") %>%
    hc_plotOptions(column = list(stacking = "normal"), 
                   series = list(borderWidth = 0)) %>%
    hc_xAxis(categories = seq(as.numeric(rownames(temp)[1]), as.numeric(rownames(temp)[nrow(temp)])),
             labels = list(style = list(fontSize = "5px"), enabled = T),
             tickLength = 0,
             title = list(text = "Age (years)", style = list(fontSize = "8px"))) %>%
    hc_yAxis(max = max(temp),
             labels = list(style = list(fontSize = "8px"), enabled = F),
             tickLength = 0,
             title = list(text = "# samples", style = list(fontSize = "8px"))) %>%
    hc_tooltip(pointFormat = HTML("{series.name}: {point.y}<br>Total: {point.stackTotal}"), 
               headerFormat = HTML("<b>{point.x} y.o.</b><br>")) %>%
    hc_title(text = paste("Total number of samples:", sum(temp), sep = " "),
             style = list(fontSize = "12px", fontWeight = "bold")) %>%
    hc_legend(verticalAlign = "top")
  
  for (i in 1:ncol(temp))
  {
    g <- g %>%
      hc_add_series(name = ifelse(colnames(temp)[i] == 1, "Male", "Female"), 
                    data = list_parse2(as.data.frame((temp[,i]))), 
                    color = ifelse(colnames(temp)[i] == 1, "#C6E2FF", "#FFB6C1"),
                    showInLegend = F)
  }
  g
}



