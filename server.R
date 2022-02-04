library(shiny)
library(ggplot2) 
library(readr) 
library(stringr)
library(tidyr)
library(Hmisc)
library(ggrepel)
library(DT)
library(ggrepel)
library(dplyr)
library(tm)
library(shinycssloaders)
library(shinythemes)
library(drawProteins)
library(magrittr)
library(writexl)
library(plotly)

############### LOAD DATA ###############
clean <- read_delim("clinvar.Jul14.pf", "\t", escape_double = FALSE, trim_ws = TRUE,
                            col_types = cols(pos_aa = col_integer()))

# clean <- read_csv("clinvar.Nov27.test.pf", col_types = cols(RowNumber = col_skip(), 
#                                                 X1 = col_skip()))
#genes     <- read_delim("genes.refSeq", "\t", escape_double = FALSE, trim_ws = TRUE)

genes     <- read_csv("HGNC_gene_names",col_types = cols_only(gene = col_character()))

sequences <- read_delim("gene-ccds-seq-length-uniprot.txt", "\t",
                        escape_double = FALSE,
                        col_types = cols(Length = col_number()),
                        trim_ws = TRUE)


shinyServer(function(input, output, session) {

############### PANEL DECISION ###############
  
  # Cambios de pestagna
  # Principal 
  observeEvent(input$goButton,  { updateTabsetPanel(session, "inTabset", selected = "panel2")})
  observeEvent(input$newsearch, { updateTabsetPanel(session, "inTabset", selected = "panel1")})
  # Resultados
  observeEvent(input$seevariants,    { updateTabsetPanel(session, "displaywhat", selected = "1")})
  observeEvent(input$seegenes,       { updateTabsetPanel(session, "displaywhat", selected = "2")})
  observeEvent(input$seephenotypes,  { updateTabsetPanel(session, "displaywhat", selected = "3")})
  observeEvent(input$seetable,       { updateTabsetPanel(session, "displaywhat", selected = "4")})  
  
############### QUERY DATASET: subset_one ###############
  query <- eventReactive(input$goButton,{
    if (input$genename %in% genes$gene) {
      return(str_trim(input$genename, side="both"))
    } else {
    return(str_trim(toupper(input$genename), side="both"))
      }
  })
  
    # poner busqueda en mayuscula
  output$displayquery <- renderText({
    paste0("\"", query(), "\"")
  })
  
  subset_one <- reactive({
      if (query() == "CLINVAR") {
      clean
    } else if (sum(grepl(paste0("^",query(),"$"),genes$gene)) >0) {
      subset(clean, grepl(paste0("^",query(),"$"), clean$GeneSymbol, ignore.case = TRUE))
    } else if (sum(grepl(paste0("^",query(),"$"), clean$Name)) >0 ) {  #variant name in HGVS format
      subset(clean, grepl(paste0("^",query(),"$"), clean$Name))
    } else {
      validate(need(grepl(query(), clean$PhenotypeList, ignore.case=TRUE), "Please repeat your search."))
      subset(clean, grepl(query(), clean$PhenotypeList, ignore.case=TRUE))
    }
  }) # es gen o no. primer grepeo
  
############### FILTER SUBSET: subset_two ###############
  clinicalinput    <- reactive({input$clinicalinput})
  consequenceinput <- reactive({input$consequenceinput})
  typeinput <- reactive({input$typeinput})
  reviewinput <- reactive({input$reviewinput})
  gnomadinput <- reactive({input$gnomad})
  caddinput <- reactive({input$cadd})
  
  query_two <- eventReactive(input$filter|input$goButton,{
    return(c(clinicalinput(),consequenceinput(),typeinput(),reviewinput(),gnomadinput(),caddinput()))
  })
  
  clinical_index <- eventReactive(input$goButton|input$filter,{
    clin_index <- c()
    if (length(query_two()) == 0) {
      clin_index <- seq(1,nrow(subset_one()))}
    else if (length(clinicalinput())== 0){
      clin_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in clinicalinput()){
        if (i == "Pathogenic"){
          clin_index <- c(clin_index, grep(7,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Likely pathogenic"){
          clin_index <- c(clin_index, grep(6,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Risk factor/Association"){
          clin_index <- c(clin_index, grep(5,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Uncertain/conflicting"){
          clin_index <- c(clin_index, grep(4,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Likely benign"){
          clin_index <- c(clin_index, grep(3,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Benign"){
          clin_index <- c(clin_index, grep(2,subset_one()$ClinicalSignificance_grouped))}
        if (i == "Protective"){
          clin_index <- c(clin_index, grep(1,subset_one()$ClinicalSignificance_grouped))}
        }}
    return(clin_index)
  })
  
  consequence_index <- eventReactive(input$goButton|input$filter,{
    con_index <- c()
    if (length(query_two()) == 0){
      con_index <- seq(1,nrow(subset_one()))}
    else if (length(consequenceinput())== 0){
      con_index <- seq(1,nrow(subset_one()))}
    else {
    for (i in consequenceinput()){
      if (i == "Missense"){
        con_index <- c(con_index, grep("Missense", subset_one()$consequence))}
      if (i == "Stop gain"){
        con_index <- c(con_index, grep("Stop gain", subset_one()$consequence))}
      if (i == "In frame indel"){
        con_index <- c(con_index, grep("In frame indel", subset_one()$consequence))}
      if (i == "Frameshift"){
        con_index <- c(con_index, grep("Frameshift", subset_one()$consequence))}
      if (i == "Synonymous"){
        con_index <- c(con_index, grep("Synonymous", subset_one()$consequence))}
      if (i == "Intronic"){
        con_index <- c(con_index, grep("Intronic", subset_one()$consequence))}
      if (i == "Non-coding"){
        con_index <- c(con_index, grep("Non-coding", subset_one()$consequence))}
      if (i == "Splice-D/A"){
        con_index <- c(con_index, grep("Splice-D/A", subset_one()$consequence))}
      if (i == "3-UTR"){
        con_index <- c(con_index, grep("3-UTR", subset_one()$consequence))}
      if (i == "5-UTR"){
        con_index <- c(con_index, grep("5-UTR", subset_one()$consequence))}
        }}
    return(con_index)
  })
  
  type_index <- eventReactive(input$goButton|input$filter,{
    type_index <- c()
    if (length(query_two()) == 0){
      type_index <- seq(1,nrow(subset_one()))}
    else if (length(typeinput())== 0){
      type_index <- seq(1,nrow(subset_one()))}
    else {
    for (i in typeinput()){
      if (i == "SNV"){
        type_index <- c(type_index, grep("SNV", subset_one()$Type))}
      if (i == "Deletion"){
        type_index <- c(type_index, grep("Deletion", subset_one()$Type))}
      if (i == "Duplication"){
        type_index <- c(type_index, grep("Duplication", subset_one()$Type))}
      if (i == "Indel"){
        type_index <- c(type_index, grep("Indel", subset_one()$Type))}
      if (i == "Inversion"){
        type_index <- c(type_index, grep("Inversion", subset_one()$Type))}
      if (i == "Insertion"){
        type_index <- c(type_index, grep("Insertion", subset_one()$Type))}
      if (i == "CNV loss"){
        type_index <- c(type_index, grep("CNV loss", subset_one()$Type))}
      if (i == "CNV gain"){
        type_index <- c(type_index, grep("CNV gain", subset_one()$Type))}
      if (i == "Complex"){
        type_index <- c(type_index, grep("Complex", subset_one()$Type))}
        }}
    return(type_index)
  })
  
  review_index <- eventReactive(input$goButton|input$filter,{
    review_index <- c()
    if (length(query_two()) == 0){
      review_index <- seq(1,nrow(subset_one()))}
    else if (length(reviewinput())== 0){
      review_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in reviewinput()){
        if (i == "No assertion criteria provided"){
          review_index <- c(review_index, grep("No assertion criteria provided", subset_one()$review))}
        if (i == "Criteria provided/ conflicting interpretations"){
          review_index <- c(review_index, grep("Criteria provided/ conflicting interpretations", subset_one()$review))}
        if (i == "Criteria provided/ single submitter"){
          review_index <- c(review_index, grep("Criteria provided/ single submitter", subset_one()$review))}
        if (i == "Practice guideline"){
          review_index <- c(review_index, grep("Practice guideline", subset_one()$review))}
        if (i == "Reviewed by expert panel"){
          review_index <- c(review_index, grep("Reviewed by expert panel", subset_one()$review))}
        if (i == "Criteria provided/ multiple submitters/ no conflicts"){
          review_index <- c(review_index, grep("Criteria provided/ multiple submitters/ no conflicts", subset_one()$review))}
      }}
    return(review_index)
  })
  
  gnomad_index <- eventReactive(input$goButton|input$filter,{
    gnomad_index <- c()
    if (length(query_two()) == 0) {
      gnomad_index <- seq(1,nrow(subset_one()))}
    else if (length(gnomadinput())== 0){
      gnomad_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in gnomadinput()){
        if (i == "Yes"){
          gnomad_index <- c(gnomad_index, grep(1,subset_one()$gnomAD_binary))}
        if (i == "No"){
          gnomad_index <- c(gnomad_index, grep(0,subset_one()$gnomAD_binary))}
      }}
    return(gnomad_index)
  })
  
  cadd_index <- eventReactive(input$goButton|input$filter,{
    cadd_index <- c()
    if (length(query_two()) == 0) {
      cadd_index <- seq(1,nrow(subset_one()))}
    else if (length(caddinput())== 0){
      cadd_index <- seq(1,nrow(subset_one()))}
    else {
      for (i in caddinput()){
        if (i == ">30"){
          cadd_index <- c(cadd_index, which(subset_one()$CADD_phred > 30))}
        if (i == ">25"){
          cadd_index <- c(cadd_index, which(subset_one()$CADD_phred > 25))}
        if (i == ">20"){
          cadd_index <- c(cadd_index, which(subset_one()$CADD_phred > 20))}
      }}
    return(cadd_index)
  })
  
  merged_index <- reactive({
    sub_index <- c()
      sub_index <- intersect(intersect(intersect(intersect(intersect(consequence_index(),clinical_index()),type_index()),review_index()),gnomad_index()),cadd_index())
    return(sub_index)
  })
  
  subset_two <- eventReactive(input$goButton|input$filter,{
    if (length(merged_index()) == 0) {
      h <- NULL
      x <- as.data.frame(h)
      return(x)
    } else {
      x <- subset_one()[merged_index(),]
      return(x)}
    })
  
  output$displayobservation <- renderText({
    if (nrow(subset_two()) == 0) {
      return(paste0("0 variants"))
    } else {
      return(paste0(nrow(subset_two())," variants"))
    }
  })

############### PANEL 1 STATISTICS ###############
  # histogram "review" ----------------------------------------------------------
  #JOINING ALL REVIEW of subset_two
  reviewgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$review, collapse=',')})
  #GGPLOT REVIEW
  output$review <- renderPlot({
      if (nrow(subset_two())==0){
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.34, y = 0.9, paste(""), 
             cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
      #all possible levels
      allreview <-   unique(clean$review)
      reviewgene_2 <- as.data.frame(t(str_split_fixed(reviewgene(), ',', Inf)))
      reviewgene_3 = data.frame(factor(reviewgene_2[,1], levels = allreview))
      reviewgene_4 = names(sort(table(reviewgene_3)))
      reviewgene_5 = data.frame(factor(reviewgene_3[,1], levels = reviewgene_4))
      ggplot(reviewgene_5, mapping = aes(x=reviewgene_5[,1])) +
        ylab("Counts") + 
        xlab("") +
        geom_bar(fill = "#43AC6A", show.legend = F) +
        geom_text(aes(label=..count..),stat="count", hjust=-0.05) +
        scale_y_continuous(labels = scales::comma, expand = expand_scale(mult = c(0, .25)))+
        scale_x_discrete(drop=FALSE, labels=c("No assertion criteria provided"="No assertion\ncriteria provided", 
                                              "no assertion criteria provided"="no assertion\ncriteria provided",
                                              "Criteria provided/ conflicting interpretations"= "Criteria provided/\nconflicting\ninterpretations", 
                                              "Criteria provided/ single submitter"="Criteria provided/\nsingle submitter", 
                                              "Practice guideline"="Practice\nguideline",
                                              "Reviewed by expert panel"="Reviewed by\nexpert panel",
                                              "Criteria provided/ multiple submitters/ no conflicts"="Criteria provided/\nmultiple submitters/\nno conflicts"))+
        coord_flip()+
        theme_light()+
        theme(title = element_blank(), text = element_text(size=15), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
  })
  
  # histogram type ----------------------------------------------------------
  #JOINING ALL TYPE of subset_two 
  typegene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$Type, collapse=',')})
  ####GGPLOT TYPE
  output$Type <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste("No variants 
                    remaining after filtering."), 
           cex = 1.5, col = "black", family="sans", font=1, adj=0.5)
    } else {
      alltype <-   unique(clean$Type)
      typegene_2 <- as.data.frame(t(str_split_fixed(typegene(), ',', Inf)))
      typegene_3 = data.frame(factor(typegene_2[,1], levels = alltype))
      typegene_4 = names(sort(table(typegene_3)))
      typegene_5 =data.frame(factor(typegene_3[,1], levels = typegene_4))
      ggplot(typegene_5,mapping = aes(x=typegene_5[,1])) +
        ylab("Counts") + 
        xlab("") +
        geom_bar(fill = "#43AC6A", show.legend = F) +
        geom_text(aes(label=..count..),stat="count", hjust=-0.05) +
        scale_y_continuous(labels = scales::comma, expand = expand_scale(mult = c(0, .25)))+
        scale_x_discrete(drop=FALSE, labels=c("Undetermined variant"="Undetermined\nvariant"))+
        coord_flip()+
        theme_light()+
        theme(title = element_blank(), text = element_text(size=15), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
  })
  
  # histogram consequence ---------------------------------------------------
  # JOINING ALL CONSEQUENCE ANNOTATIONS  (FROM ALL OBSERVATIONS OF ONE GENE/DISEASE)
  consgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$consequence, collapse=',')})   
  # GGPLOT CONSEQUENCE
  output$consequence <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste(""), 
           cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
      allcons <- as.character(unique(clean$consequence))
      consgene_2 <- as.data.frame(t(str_split_fixed(consgene(), ',', Inf)))
      #consgene_2[,1] <- gsub(NA, "NA", consgene_2[,1])
      consgene2a <- factor(consgene_2$V1, levels = allcons)
      consgene_3 = data.frame(consgene2a)
      consgene_4 = names(sort(table(consgene_3)))
      consgene_5 =data.frame(factor(consgene_3[,1], levels = consgene_4))
      #consgene_5 = data.frame(consgene_5a[is.na(consgene_5a)== F])
      
    ggplot(consgene_5,mapping = aes(x=consgene_5[,1], na.rm=TRUE)) +
      ylab("Counts") + 
      xlab("") +
      geom_bar(fill = "#43AC6A", show.legend = F)+
      geom_text(aes(label=..count..),stat="count", hjust=-0.05) +
      scale_x_discrete(drop=FALSE)+
      scale_y_continuous(labels = scales::comma, expand = expand_scale(mult = c(0, .25)))+
      coord_flip() +
      theme_light() +
      theme(title = element_blank(), text = element_text(size=15), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
  })
  
  # histogram clinical significance -----------------------------------------
  # JOINING ALL CLINICAL SIGNIFICANCE ANNOTATIONS FROM ONE GENE: clinicalgene_3$V1
  clinicalgene <- eventReactive(input$filter|input$goButton,{paste0(subset_two()$ClinicalSignificance_grouped, collapse=',')})
  # GGPLOT CLINICAL SIGNIFICANCE
  output$ClinicalSignificance <- renderPlot({
    if (nrow(subset_two())==0){
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.34, y = 0.9, paste(""), 
           cex = 1.5, col = "black", family="serif", font=2, adj=0.5)
    } else {
      clinicalgene_2 <- as.data.frame(t(str_split_fixed(clinicalgene(), ',', Inf)))
      clinicalgene_3 = data.frame(factor(clinicalgene_2[,1], levels = c("1","2","3","4","5","6","7","NA")))
      clinicalgene_4 = names(sort(table(clinicalgene_3)))
      clinicalgene_5a = factor(clinicalgene_3[,1], levels = clinicalgene_4)
      clinicalgene_5 = data.frame(clinicalgene_5a[is.na(clinicalgene_5a)== F])
      
    ggplot(clinicalgene_5, aes(x=clinicalgene_5[,1])) +
      ylab("Counts") + 
      xlab("") +
      geom_bar(fill = "#43AC6A", show.legend = F) +
      geom_text(aes(label=..count..),stat="count", hjust=-0.05) +
      scale_y_continuous(labels = scales::comma, expand = expand_scale(mult = c(0, .25)))+
      scale_x_discrete(drop=FALSE, labels=c("1"="Protective","2"="Benign", "3"="Likely\nbenign", "4"="Uncertain/\nconflicting", 
                                            "5"="Risk factor/\nAssociation","6"="Likely\npathogenic","7"="Pathogenic","NA"="NA"))+
      theme_light() +
      theme(title = element_blank(), text = element_text(size=15), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      coord_flip()
    }
  })
############### PANEL 2 PROTEINPLOT/ GENECOUNT PLOT  ###############
  ### decision: gene/HGVS vs. disease query
  output$displayngenes <- renderText({
    if ((sum(grepl(paste0("^",query(),"$"),genes$gene)) >0)|(sum(grepl(paste0("^",query(),"$"),clean$Name)) >0)) {
      return(paste0("Protein Mapping"))
    } else {
      geneNAME0  <- subset(subset_two(), !grepl("subset of|covers", subset_two()$GeneSymbol, ignore.case=TRUE)) 
      geneNAME <-  gsub("^-$","Intergenic", paste0(geneNAME0$GeneSymbol, collapse = ", "))
      gene0    <- as.data.frame(t(str_split(geneNAME, pattern = ", ", simplify = T)), stringsAsFactors = T)
      gene1    <- as.data.frame(table(gene0[,1]))
      return(paste0(nrow(gene1), " genes"))
      }
     })
  
    proteinsequence <- reactive({
      if (sum(grepl(paste0("^",query(),"$"),genes$gene)) >0){ #decision between gene name and HGVS nomenclature
        return(sequences[sequences$Gene_ID == query(),])
      } else{
        return(sequences[sequences$Gene_ID == subset_two()$GeneSymbol[1],])}
    })
    uniprot_code    <- reactive({
      if (sum(grepl(paste0("^",query(),"$"),genes$gene)) >0){ #decision between gene name and HGVS nomenclature
        sequences[sequences$Gene_ID == query(),7]
      } else{
        return(sequences[sequences$Gene_ID == subset_two()$GeneSymbol[1],7])}
    })
  
  output$genefield <- renderUI({
    
    if ((sum(grepl(paste0("^",query(),"$"),genes$gene)) > 0)|(sum(grepl(paste0("^",query(),"$"),clean$Name)) > 0)) { 

      #uniprot_data
        uniprot_code()  %>%
        drawProteins::get_features() %>%
        drawProteins::feature_to_dataframe() -> prot_data
        prot_data = prot_data[prot_data$type == "DOMAIN"|prot_data$type == "REGION"|
                              prot_data$type == "REPEAT"|prot_data$type == "MOTIF"| prot_data$type == "NP_BIND" , ]
        prot_data = prot_data[prot_data$type == "DOMAIN"|prot_data$type == "REGION"|
                                prot_data$type == "REPEAT"|prot_data$type == "MOTIF"| prot_data$type == "NP_BIND" , ]
        prot_data = prot_data[order(prot_data$length, decreasing=T),]
        
        #https://www.w3.org/TR/css-color-3/#svg-color
        #cols <- c("crimson", "darkseagreen",  "slateblue","mediumvioletred", "mediumturquoise","goldenrod", "sandybrown", "thistle")
        
        cols <- c("#96ceb4", "ffcc5c", 	"#92b2ff", "crimson","goldenrod", "sandybrown", "thistle",  "ffeead")
        prot_data$color <- factor(prot_data$type, labels = cols[1:length(unique(prot_data$type))])
        
        #uniprot_data_end
        length= sequences$Length[sequences$Gene_ID == query()]
        axis_template_x <- list(showgrid = F,zeroline = F , showline = FALSE,
                                 showticklabels = T, range = c(-10, length+10))
        
        axis_template_y <- list(showgrid = F,zeroline = F ,showline = FALSE,
                                 showticklabels = FALSE, range = c(-2.5, 2.5))
    
      #only consequences which can be plotted: missense,PTV,synonymous -> column: consequence_plot
      consequence_plot = as.character(subset_two()$consequence_plot)
      #Y-axis values of variants: variant dots in different height depending on consequence
      height_variant_factor <- factor(consequence_plot,levels=c("Synonymous","Missense","PTV"), labels = c("0.6", "1", "1.4"))
      #only variants as input for the plot which are PTV, missense or synonymous, no NA
      height_variant_na <- as.numeric(as.character(height_variant_factor))[!is.na(consequence_plot)]
      subset_twoconsequence_K_na <- as.factor(consequence_plot)[!is.na(consequence_plot)]
      labelposition_na <- as.numeric(subset_two()$pos_aa)[!is.na(consequence_plot)] #X-axis values of variants
      subset_plot = subset_two()[!is.na(consequence_plot),]
      subset_plot <- data.frame(subset_plot, height_variant_na, subset_twoconsequence_K_na, labelposition_na)
      subset_plot = subset_plot[!(subset_plot$labelposition_na == "NA"|is.na(subset_plot$labelposition_na)),]
      
      subset_plot$height_variant_na[subset_plot$gnomAD_binary == "1" & !is.na(subset_plot$gnomAD_binary)] <- as.numeric(as.character(subset_plot$height_variant_na[subset_plot$gnomAD_binary == "1" & !is.na(subset_plot$gnomAD_binary)]))*(-1)
        
      
      #subset_plot$height_variant_na <- ifelse((subset_plot$gnomAD_binary == "1" & !is.na(subset_plot$gnomAD_binary)), as.numeric(as.character(subset_plot$height_variant_na[subset_plot$gnomAD_binary == "1"]))*(-1), as.numeric(as.character(subset_plot$height_variant_na[subset_plot$gnomAD_binary == "0"])))
      
      
      #shapes
      lines <- list()
      
        if (nrow(subset_plot)>0){
          for (i in 1:nrow(subset_plot)){
            #if(is.na(subset_plot$gnomAD_binary[i])|(subset_plot$gnomAD_binary[i] == "0")){
              newlist = list(type = "line",line = list(color = "black"),xref = "x",yref = "y",layer="below",
                             y0 = 0, y1 = subset_plot$height_variant_na[i], x0=subset_plot$labelposition_na[i], x1=subset_plot$labelposition_na[i])
              lines = c(lines, list(newlist)) 
            #} else {
            #  newlist = list(type = "line",line = list(color = "black"),xref = "x",yref = "y",layer="below",
            #                 y0 = 0, y1 = -subset_plot$height_variant_na[i], x0=subset_plot$labelposition_na[i], x1=subset_plot$labelposition_na[i])
            #  lines = c(lines, list(newlist)) 
            #} 
          }}
       
       if (nrow(prot_data)>0){
         for (i in 1:length(prot_data$color)){
           if (prot_data$type[i]== "REPEAT"){
             width=0.12
           } else if (prot_data$type[i]== "REGION") {
             width=0.12
           }else if (prot_data$type[i]== "DOMAIN") {
             width=0.12
           }else if (prot_data$type[i]== "MOTIF") {
             width=0.12
           }else if (prot_data$type[i]== "NP_BIND") {
             width=0.12}
           newlist = list(type = "rect", fillcolor = prot_data$color[i], line = list(color ='rgba(67,67,67,1)' ), opacity = 1,
                          x0 = prot_data$begin[i], x1 = prot_data$end[i], xref = "x",
                          y0 = width*(-1), y1 = width, yref = "y", name= paste0("rect",i))
           lines = c(lines, list(newlist))
         }}
       
       output$GeneCountorProteinseq = renderPlotly({ 
       
         p <- plot_ly(subset_plot)%>% layout(xaxis = axis_template_x ,
                                             yaxis = axis_template_y ,
                                             annotations = list(x= labelposition_na, y= height_variant_na,
                                                          font = list(size = 20),xref="x",yref="y",text= subset_plot$Link_space,showarrow=FALSE),
                                             shapes = c(list(list(type = "rect", fillcolor = "silver", line = list(color ='rgba(67,67,67,1)'), opacity = 1,
                                                                  x0 = 0, x1 = length, xref = "x",
                                                                  y0 = -0.18, y1 = 0.18, yref = "y", name="proteinseq.")), lines)) %>% 
                                      config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "hoverClosestGeo", "hoverClosestGl2d","toImage",
                                             "hoverClosestCartesian","lasso2d","select2d","resetScale2d",
                                             "hoverCompareCartesian", "hoverClosestPie","toggleSpikelines"), displaylogo = FALSE) 
         
         subset_plot$gnomAD = rep(NA, length(subset_plot$gnomAD_binary))
         
         #if((subset_plot$gnomAD_binary[i] == 1) & (!is.na(subset_plot$gnomAD_freq[i]))) {
         
         for(i in 1:nrow(subset_plot)){
           if ((subset_plot$gnomAD_binary[i] == 1) & (!is.na(subset_plot$gnomAD_freq[i]))) {
               if (as.numeric(subset_plot$gnomAD_freq[i]) < 1) {
                  subset_plot$gnomAD[i] = paste0("Yes (freq= ",subset_plot$gnomAD_freq[i],")") }
                else { subset_plot$gnomAD[i] = "Yes"}
            } else if (is.na(subset_plot$gnomAD_binary[i])) {
              subset_plot$gnomAD[i] = "N/A"
            } else if (subset_plot$gnomAD_binary[i] == 0){
             subset_plot$gnomAD[i] = "No"
             }
           }

          
        #Symbols to show gnomad Y/N
        #p <- add_trace(p, x= subset_plot$labelposition_na[subset_plot$gnomAD_binary == 1],
        #                y=subset_plot$height_variant_na[subset_plot$gnomAD_binary == 1],
        #                type = "scatter", mode="markers", name="In gnomAD",hoverinfo = "none",
        #                marker = list(symbol = 1 ,size = 10, color="black",opacity = 1)) 
         
        # p <- add_trace(p, x= subset_plot$labelposition_na[subset_plot$gnomAD_binary == 0],
        #                y=subset_plot$height_variant_na[subset_plot$gnomAD_binary == 0],
        #                type = "scatter", mode="markers", name="Not in gnomAD",hoverinfo = "none",
        #                marker = list(symbol = 13,size = 14, color="black",opacity = 1)) 
         
         p <- add_trace(p, x= subset_plot$labelposition_na[subset_plot$subset_twoconsequence_K_na == "Synonymous"],
                        y=subset_plot$height_variant_na[subset_plot$subset_twoconsequence_K_na == "Synonymous"],
                        type = "scatter", mode="markers",marker = list(size = 10, color="skyblue"),name="Synonymous",size=10, text=paste0("Variant: ",subset_plot$Name[subset_plot$subset_twoconsequence_K_na == "Synonymous"], "\n",
                                                                                                                                         "Clinical Significance: ", subset_plot$ClinicalSignificance[subset_plot$subset_twoconsequence_K_na == "Synonymous"], "\n",
                                                                                                                                         "Review Status: ", subset_plot$review[subset_plot$subset_twoconsequence_K_na == "Synonymous"],"\n",
                                                                                                                                         "CADD score: ", subset_plot$CADD_phred[subset_plot$subset_twoconsequence_K_na == "Synonymous"],"\n",
                                                                                                                                         "In gnomAD: ", subset_plot$gnomAD[subset_plot$subset_twoconsequence_K_na == "Synonymous"]), showlegend = T, hoverinfo = 'text', opacity = 1)
         p <- add_trace(p, x= subset_plot$labelposition_na[subset_plot$subset_twoconsequence_K_na == "Missense"], 
                        y=subset_plot$height_variant_na[subset_plot$subset_twoconsequence_K_na == "Missense"],
                        type = "scatter", mode="markers",marker = list(size = 10, color="salmon"),name="Missense",size=10, text=paste0("Variant: ",subset_plot$Name[subset_plot$subset_twoconsequence_K_na == "Missense"], "\n",
                                                                                                                                        "Clinical Significance: ", subset_plot$ClinicalSignificance[subset_plot$subset_twoconsequence_K_na == "Missense"], "\n",
                                                                                                                                        "Review Status: ", subset_plot$review[subset_plot$subset_twoconsequence_K_na == "Missense"],"\n",
                                                                                                                                         "CADD score: ", subset_plot$CADD_phred[subset_plot$subset_twoconsequence_K_na == "Missense"],"\n",
                                                                                                                                        "In gnomAD: ", subset_plot$gnomAD[subset_plot$subset_twoconsequence_K_na == "Missense"]),showlegend = T,hoverinfo = 'text', opacity = 1) 
         p <- add_trace(p, x= subset_plot$labelposition_na[subset_plot$subset_twoconsequence_K_na == "PTV"], 
                        y=subset_plot$height_variant_na[subset_plot$subset_twoconsequence_K_na == "PTV"],
                        type = "scatter", mode="markers",marker = list(size = 10, color="dimgray"),name="PTV",size=10, text=paste0("Variant: ",subset_plot$Name[subset_plot$subset_twoconsequence_K_na == "PTV"], "\n",
                                                                                                                                  "Clinical Significance: ", subset_plot$ClinicalSignificance[subset_plot$subset_twoconsequence_K_na == "PTV"], "\n",
                                                                                                                                   "Review Status: ", subset_plot$review[subset_plot$subset_twoconsequence_K_na == "PTV"],"\n",
                                                                                                                                  "CADD score: ", subset_plot$CADD_phred[subset_plot$subset_twoconsequence_K_na == "PTV"],"\n",
                                                                                                                                  "In gnomAD: ", subset_plot$gnomAD[subset_plot$subset_twoconsequence_K_na == "PTV"]),showlegend = T,hoverinfo = 'text', opacity = 1) 

         #domain description
         if ( nrow(prot_data)>0){
           for (i in 1:length(prot_data)){
             p <- add_trace(p, x=(prot_data$begin[i]+prot_data$end[i])/2, y=0, mode="markers", type="scatter",
                            marker = list(size= 8, color = prot_data$color[i]),hoverinfo = 'text',text=prot_data$description[i], name=prot_data$description[i])
           }}
         
         return(p)})
      
     fluidRow(column(width = 12, h5(strong(paste0("Coding variant mapping and domain annotation for ",query())),br(), h6("(Above sequence: Variants not present in gnomAD/ Below: Variants also in gnomAD)", align = "center"), align = "center",br(),br(),
                    withSpinner(plotlyOutput("GeneCountorProteinseq")))))
    } else {
      
      ## if disease: gene count 
      geneNAME0  <- subset(subset_two(), !grepl("subset of|covers", subset_two()$GeneSymbol, ignore.case=TRUE)) 
      geneNAME <-  gsub("^-$","Intergenic", paste0(geneNAME0$GeneSymbol, collapse = ", "))
      gene0    <- as.data.frame(t(str_split(geneNAME, pattern = ", ", simplify = T)), stringsAsFactors = T)
      gene1    <- as.data.frame(table(gene0[,1]))
      gene2    <- head(gene1[order(-gene1$Freq),], 10)
      gene3    <- gene1[order(-gene1$Freq),]
    
      output$GeneCountor = renderPlot({ 
        ggplot(gene2, mapping = aes(x = reorder(gene2$Var1, gene2$Freq), y = gene2$Freq) ) +
        geom_bar(stat = "identity", fill = "#F04124", show.legend = F, width = 0.7) +
        ylab("Counts") +
        theme_light() +
        theme(title = element_blank(),
              axis.title.x=element_blank(), text = element_text(size=15))})
      
      output$GeneCountTotal <- DT::renderDataTable({
        colnames(gene3) <- c("Gene", "Frequency")
        DT::datatable(gene3, 
                      #extensions = c('Scroller'), 
                      rownames = F, fillContainer = T,
                      options = list(dom = "t",
                                     scrollY = 350,
                                     scroller = TRUE,
                                     pageLength = nrow(gene1)))})
        
       fluidRow(column(width = 8, h5(strong("Top 10 genes associated")), align = "center", br(), 
                      withSpinner(plotOutput("GeneCountor"))),
               column(width = 4, h5(strong("Total genes associated")), align = "center",
                      wellPanel(DT::dataTableOutput("GeneCountTotal"),
                                style = "background-color: #ffffff;
                                border-color: #ffffff;
                                box-shadow: 0px 0px 0px #ffffff;
                                margin-bottom: 5px")
                      )
                      )
    }
    })
  
  
  
  
  # histogram: https://plot.ly/r/shinyapp-linked-brush/
  
  
############### PANEL 3 PHENOTYPE ###############
  # phenotype count 
  phenoNAME <- reactive({ paste0(subset_two()$PhenotypeList, collapse = ";") })
  pheno0    <- reactive({ as.data.frame(t(str_split(phenoNAME(), pattern = ";", simplify = T)), stringsAsFactors = T) })
  pheno1    <- reactive({ as.data.frame(table(pheno0()$V1)) })
  pheno2    <- reactive({ head( pheno1()[order(- pheno1()$Freq),], 10) })
  pheno3    <- reactive({ pheno1()[order(- pheno1()$Freq),] })
  
  output$displayphenotypes <- renderText({
    return(paste0(nrow(pheno1()), " phenotypes"))
  })
  
  output$DiseaseCount <- renderPlot({
    ggplot(pheno2(), mapping = aes(x = reorder(stringr::str_wrap(pheno2()$Var1, 15), pheno2()$Freq), y = pheno2()$Freq) ) +
      geom_bar(stat = "identity", fill = "#E99002", show.legend = F, width = 0.7) +
      ylab("Counts") +
      theme_light() +
      theme(title = element_blank(),
            axis.title.x= element_blank(), text = element_text(size=15))
  })
  
  output$DiseaseCountTotal <- DT::renderDataTable({
    dat3 <- pheno3()
    colnames(dat3) <- c("Phenotype", "Frequency")
    DT::datatable(dat3, 
                  #extensions = c('Scroller'), 
                  rownames = F, 
                  fillContainer = T,
                  options = list(dom = "t",
                                 scrollY = 350,
                                 scroller = TRUE,
                                 pageLength = nrow(pheno1())
                  ))
  })

############### PANEL 4 TABLE ############### 
  #Table Final: Variantstable
  
  subset_two_gnomAD = reactive({data.frame(subset_two(),"gnomAD_binary_char" = ifelse(subset_two()$gnomAD_binary == 0, "No","Yes"))})
  subset_two_table <- reactive({select(subset_two_gnomAD(),"Link", "GeneSymbol", "Type", "consequence", "ClinicalSignificance","review","PhenotypeList", "Name", "ref_aa", "alt_aa", "pos_aa", "CADD_phred", "gnomAD_binary_char") })
  subset_two_SNV <- reactive({subset_two_table()[subset_two_table()$Type == "SNV",]})
  subset_two_noSNV <- reactive({subset_two_table()[subset_two_table()$Type != "SNV",]})
  
  
  output$clinvartable1 <- DT::renderDataTable({
    DT::datatable(subset_two_SNV(), 
                  extensions = c("Buttons"), 
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    'Table 1: ', htmltools::em('ClinVar single nucleotide variants (SNV) subset')),
                  rownames = F, colnames = c("VariationID (Clinvar)",'Gene', 'Type', "Consequence", "Clinical significance","Review status", 'Phenotype', 'HGVS', "Ref_aa", "Alt_aa", "Pos_aa", "CADD score", "gnomAD"), 
                  fillContainer = T,
                  escape=FALSE,
                  options = list(dom = "t", 
                                scrollY = 250,
                                scrollX= TRUE,
                                scroller = TRUE,
                                #pageLength = nrow(subset_two()),
                                pageLength = 5000,
                                autoWidth = TRUE
                                #,
                                #columnDefs = list(list(width = '10%', targets = c(1,5)))
                       ))
  }, escape = FALSE)
  
  output$clinvartable2 <- DT::renderDataTable({
    DT::datatable(subset_two_noSNV(), 
                  #extensions=c("Buttons","Scroller"),
                  extensions = c("Buttons"),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    'Table 2: ', htmltools::em('ClinVar Copy number variants (CNV) subset')),
                  rownames = F, colnames = c("VariationID (Clinvar)",'Gene', 'Type', "Consequence", "Clinical significance","Review status",'Phenotype', 'HGVS', "Ref_aa", "Alt_aa", "Pos_aa", "CADD score", "gnomAD"), 
                  fillContainer = T,
                  escape=FALSE,
                  options = list(dom = "t", 
                                 scrollY = 250,
                                 scrollX= TRUE,
                                 scroller = TRUE,
                                 #pageLength = nrow(subset_two()),
                                 pageLength = 5000,
                                 autoWidth = TRUE
                                 #, 
                                 #columnDefs = list(list(width = "50px" , targets = list(2,8)))
                  ))
    #https://stackoverflow.com/questions/55165477/r-shiny-dashboard-datatable-column-width
  }, escape = FALSE)
  
  #Download of tables
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste0("ShinyClinVar_", query(), "_" ,Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(subset_two_SNV(), file)
    }
  )
  
  output$downloadData1_xlsx <- downloadHandler(
    filename = function() {
      paste0("ShinyClinVar_", query(), "_" ,Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write_xlsx(subset_two_SNV(), file)
    }
  )

  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0("ShinyClinVar_", query(), "_" ,Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(subset_two_noSNV(), file)
    }
  )
  
  output$downloadData2_xlsx <- downloadHandler(
    filename = function() {
      paste0("ShinyClinVar_", query(), "_" ,Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write_xlsx(subset_two_noSNV(), file)
    }
  )
############### ABOUT ###############
  output$image1 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure1.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$image21 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure2.1.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$image22 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure2.2.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$image23 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure2.3.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$image31 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure3.1.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$image32 <- renderImage({
    imagefilename <- normalizePath(file.path('./www/Figure3.2.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  output$imageNCBI <- renderImage({
    imagefilename <- normalizePath(file.path('./www/NCBI_powered.png'))
    list(src = imagefilename)
  }, deleteFile = FALSE)
  
############### TEST ###############
   # output$test <- renderPrint({
   #  #if (query() %in% genes$gene | (sum(grepl(query(), clean$Name, fixed = TRUE))>0) ) {
   #  #uniprot_code()  %>%
   #  #drawProteins::get_features() %>%
   #  #drawProteins::feature_to_dataframe() -> prot_data}
   # 
   #   #str(subset_two()$consequence_plot)
   #   #return(levels(as.factor(prot_data$type)))
   #   #return(levels(as.factor(prot_data$description)))
   #   #print(max(prot_data$length))
   #   #return(head(prot_data[prot_data$type == "MOD_RES",]))
   #  
   #   consequence_plot = as.character(subset_two()$consequence_plot)
   #   #Y-axis values of variants: variant dots in different height depending on consequence
   #   height_variant_factor <- factor(consequence_plot,levels=c("Synonymous","Missense","PTV"), labels = c("0.6", "1", "1.4"))
   #   #only variants as input for the plot which are PTV, missense or synonymous, no NA
   #   height_variant_na <- as.numeric(as.character(height_variant_factor))[!is.na(consequence_plot)]
   #   subset_twoconsequence_K_na <- as.factor(consequence_plot)[!is.na(consequence_plot)]
   #   labelposition_na <- as.numeric(subset_two()$pos_aa)[!is.na(consequence_plot)]
   #  
   #    
   #   subset_plot = subset_two()[!is.na(consequence_plot),]
   #   subset_plot <- data.frame(subset_plot, height_variant_na, subset_twoconsequence_K_na, labelposition_na)
   #   
   #   print(height_variant_na)
     
  #})
})
