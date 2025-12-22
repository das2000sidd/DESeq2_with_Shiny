#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(ggplot2)
library(ggrepel)
library(grid)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)
library(gridExtra)

`%ni%` <- Negate(`%in%`);

## Add enrichment analysis
safe_mapIds <- function(db, keys, keytype, column) {
  
  if (is.null(keys)) return(rep(NA, length(keys)))
  
  keys <- as.character(keys)
  #keys <- keys[!is.na(keys) & keys != ""]
  
  res <- rep(NA_character_, length(keys))   # initialize with NAs
  
  mapped <- AnnotationDbi::mapIds(
    db,
    keys = keys,
    keytype = keytype,
    column = column,
    multiVals = "first"
  )
  
  # initialize result vector with NA
  #res <- rep(NA_character_, length(keys))
  
  # fill in mapped values, matching names
  #names(res) <- keys
  #res[names(mapped)] <- mapped
  res[match(names(mapped), keys)] <- mapped
  res
}

## Actual logic of implementation is here
server <- function(input, output, session) {
  
  # Robust file reader for DE table (CSV or tab-delim)
  read_de <- function(file, sep_choice) {
    req(file)
    ext <- tolower(tools::file_ext(file$name))
    
    reader <- if (ext == "csv") read.csv else read.table
    
    # prefer explicit .csv check, otherwise use sep control
    df <- tryCatch(
      reader(
        file = file$datapath,
        sep = input$sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        row.names = NULL   # 🔑 CRITICAL
      ),
      error = function(e) {
        showNotification(
          "Failed to read file. Please upload a file without row names.",
          type = "error"
        )
        return(NULL)
      }
    )
    df
  }
  
  # Load results (CSV or TXT)
  res <- reactive({
    req(input$de_file)
    df <- read_de(input$de_file, 
                  input$sep)
    
    df
})
  
  cpm <- reactive({
    req(input$cpmfile)
    
    file <- input$cpmfile
    ext <- tools::file_ext(file$name)
    if (tolower(ext) == "csv") {
      read.csv(file$datapath, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      read.table(file$datapath, sep = input$sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    }
  })
  
  
  filtered_res <- reactive({
    req(res())
    de_tab <- res();
    
    if (nrow(de_tab) == 0) {
      return(list(de = de_tab))
    }
    
    validate(
      need("Ensembl" %in% colnames(de_tab),
           "Filtered DE table must contain an 'Ensembl' column")
    )
    
    validate(
      need(any(nzchar(as.character(de_tab$Ensembl))),
           "No valid Ensembl IDs found after filtering")
    )
    
    
    validate(
      need("padj" %in% colnames(de_tab) && "log2FoldChange" %in% colnames(de_tab),
           "Uploaded DE table must contain columns 'padj' and 'log2FoldChange'.")
    )
    
    de_tab_up <- subset(de_tab, de_tab$log2FoldChange > input$log2FoldChange & de_tab$padj < input$padj);
    de_tab_dn <- subset(de_tab, de_tab$log2FoldChange < - input$log2FoldChange & de_tab$padj < input$padj);
    
    de_tab_sig <- rbind(de_tab_up, de_tab_dn);
    
    if(input$org == "Homo sapiens"){
      species <- org.Hs.eg.db
    } else {
      species <- org.Mm.eg.db
    }
    
    validate(
      need("Ensembl" %in% colnames(de_tab_sig),
           "DE table must contain an 'Ensembl' column for annotation")
    )
    
    #keys_ensembl <- as.character(de_tab_sig$Ensembl)
    #keys_ensembl <- keys_ensembl[!is.na(keys_ensembl) & keys_ensembl != ""]
    
    #de_tab_sig$Entrez <- NA
    de_tab_sig$Entrez <- tryCatch(
      safe_mapIds(
        species,
        keys = de_tab_sig$Ensembl,
        keytype = "ENSEMBL",
        column = "ENTREZID"
      ),
      error = function(e) {
        showNotification(
          paste0(
            "It seems the Ensembl IDs in your DE table do not match the selected species (",
            input$org, ").\n",
            "Please check your species selection or convert your gene IDs."
          ),
          type = "error",
          duration = 10
        )
        # return NA vector of correct length to avoid breaking the table
        rep(NA_character_, nrow(de_tab_sig))
      }
    )
    
    de_tab_sig$Symbol <- tryCatch(
      safe_mapIds(
        species,
        keys = de_tab_sig$Ensembl,
        keytype = "ENSEMBL",
        column = "SYMBOL"
      ),
      error = function(e) {
        showNotification(
          paste0(
            "It seems the Ensembl IDs in your DE table do not match the selected species (",
            input$org, ").\n",
            "Please check your species selection or convert your gene IDs."
          ),
          type = "error",
          duration = 10
        )
        # return NA vector of correct length to avoid breaking the table
        rep(NA_character_, nrow(de_tab_sig))
      }
    )
    
    de_tab_sig$Gene_Name <- tryCatch(
      safe_mapIds(
        species,
        keys = de_tab_sig$Ensembl,
        keytype = "ENSEMBL",
        column = "GENENAME"
      ),
      error = function(e) {
        showNotification(
          paste0(
            "It seems the Ensembl IDs in your DE table do not match the selected species (",
            input$org, ").\n",
            "Please check your species selection or convert your gene IDs."
          ),
          type = "error",
          duration = 10
        )
        # return NA vector of correct length to avoid breaking the table
        rep(NA_character_, nrow(de_tab_sig))
      }
    )
    
    de_tab_sig$Entrez <- unname(as.character(de_tab_sig$Entrez))
    de_tab_sig$Symbol <- unname(as.character(de_tab_sig$Symbol))
    de_tab_sig$Gene_Name <- unname(as.character(de_tab_sig$Gene_Name))
    
    de_tab_sig <- de_tab_sig %>% mutate(across(everything(), ~ as.character(.)))
    
    rownames(de_tab_sig) <- NULL
    
    list(
      de = de_tab_sig
    )
    
  })
  
  output$results_table <- renderDT({
    req(filtered_res())
    ## This function creates an HTML widget to display rectangular data
    datatable(
      filtered_res()$de,
      rownames = FALSE,
      extensions = c("Scroller", "Buttons"),
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        deferRender = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      )
    )
  })
  
  volcano_plot <- reactive({
    de_tab <- res()
    
    validate(need(nrow(de_tab) > 0, "DE table is empty"))
    validate(need("padj" %in% colnames(de_tab) && "log2FoldChange" %in% colnames(de_tab),
                  "DE table must contain 'padj' and 'log2FoldChange' for volcano plot."))
    
    # create gene and Symbol fallbacks
    if (!"Symbol" %in% colnames(de_tab) && "Ensembl" %in% colnames(de_tab)) {
      
      species <- if (input$org == "Homo sapiens") org.Hs.eg.db else org.Mm.eg.db
      
      keys <- as.character(de_tab$Ensembl)
      keys <- keys[!is.na(keys) & keys != ""]
      
      if (length(keys) > 0) {
        de_tab$Entrez <- safe_mapIds(
          species,
          keys,
          keytype = "ENSEMBL",
          column = "ENTREZID"
        )
        
        de_tab$Symbol <- safe_mapIds(
          species,
          as.character(de_tab$Entrez),
          keytype = "ENTREZID",
          column = "SYMBOL"
        )
      }
    }

    # coerce numeric columns
    de_tab$padj <- as.numeric(de_tab$padj)
    de_tab$log2FoldChange <- as.numeric(de_tab$log2FoldChange)
    
    # safe - replace NA padj with 1 so -log10 won't produce Inf
    de_tab$padj[is.na(de_tab$padj)] <- 1
    
    de_tab$Neg_log_p_val <- -log10(ifelse(is.na(de_tab$padj), 1, de_tab$padj))
    
    de_tab_up <- subset(de_tab, de_tab$log2FoldChange > input$log2FoldChange & de_tab$padj < input$padj)
    de_tab_dn <- subset(de_tab, de_tab$log2FoldChange < - input$log2FoldChange & de_tab$padj < input$padj)
    
    if(nrow(de_tab_up) > 0 & nrow(de_tab_dn) > 0){
    
    de_tab_sig <- rbind(de_tab_up, de_tab_dn)
    
    } else if(nrow(de_tab_up)== 0 & nrow(de_tab_dn) > 0){
      
      de_tab_sig <- de_tab_dn
      
    } else if(nrow(de_tab_up)!= 0 & nrow(de_tab_dn) == 0){
      
      de_tab_sig <- de_tab_up
      
    }
    
    #de_tab <- left_join(de_tab, de_tab_sig, by = "gene")
    
    #de_tab$Direction[is.na(de_tab$Direction)] <- "No sig change"
    
    #de_tab$baseMean_log <- log2(de_tab$baseMean + 1)
    
    #de_tab_sig <- subset(de_tab, abs(log2FoldChange) > input$log2FoldChange & padj < input$padj)
    #de_tab_sig$de_tab_Direction <- ifelse(de_tab_sig$log2FoldChange > 0, "Up", "Down")
    
    #de_tab_up <- subset(de_tab_sig, log2FoldChange > 0)
    #de_tab_dn <- subset(de_tab_sig, log2FoldChange < 0)
    
    de_tab_up <- de_tab_up[order(de_tab_up$padj),]
    de_tab_dn <- de_tab_dn[order(de_tab_dn$padj),]
    
    top_15_up <- head(de_tab_up, 10)
    top_15_dn <- head(de_tab_dn, 10)

    top_15_up_dn <- rbind(top_15_up, top_15_dn)
    top_15_up_dn <- subset(top_15_up_dn, !is.na(Symbol) & Symbol != "NULL")
    
    # compute finite ymax for plotting (avoid Inf/NA)
    finite_vals <- de_tab$Neg_log_p_val[is.finite(de_tab$Neg_log_p_val)]
    if (length(finite_vals) == 0) {
      ymax <- 1
    } else {
      ymax <- max(finite_vals, na.rm = TRUE)
      # cap extremely large values to avoid crazy axes (optional)
      ymax <- min(ymax, quantile(finite_vals, 0.999, na.rm = TRUE) * 1.2)
    }
    
    ## Coord for no of DE annotation
    # choose safe x,y positions for annotation using quantiles (robust)
    x_left  <- as.numeric(quantile(de_tab$log2FoldChange, probs = 0.02, na.rm = TRUE)) ;
    x_right <- as.numeric(quantile(de_tab$log2FoldChange, probs = 0.98, na.rm = TRUE)) ;
    x_annot <- x_left + 0.2 * (x_right - x_left)  # slightly inset from left edge
    
    y_annot_top <- ymax * 0.95 ;
    y_annot_bot <- ymax * 0.87 ;
    
    p <- ggplot(de_tab, aes(x = log2FoldChange, y = Neg_log_p_val)) +
      theme_classic(base_size = 16) +
      geom_point(color = "grey70", size = 2, na.rm = TRUE)
    
    if (nrow(de_tab_up) > 0) {
      p <- p + geom_point(data = de_tab_up, aes(x = log2FoldChange, y = Neg_log_p_val),
                          size = 3, color = "red", na.rm = TRUE)
    }
    if (nrow(de_tab_dn) > 0) {
      p <- p + geom_point(data = de_tab_dn, aes(x = log2FoldChange, y = Neg_log_p_val),
                          size = 3, color = "blue", na.rm = TRUE)
    }
    
    p <- p +
      ggtitle("Volcano plot") +
      theme(plot.title = element_text(hjust = 0.5, size = 12)) +
      xlab("Log2FC") + ylab("-Log10(adj. p)")
    
    # add annotation texts (counts) only if ymax > 0 so they sit inside the plot
    if (is.finite(ymax) && ymax > 0) {
      p <- p +
        annotate("text",
                 x = x_annot,
                 y = y_annot_top,
                 label = paste0(nrow(de_tab_up), " genes up"),
                 color = "red",
                 size = 4,
                 hjust = 0) +
        annotate("text",
                 x = x_annot,
                 y = y_annot_bot,
                 label = paste0(nrow(de_tab_dn), " genes down"),
                 color = "blue",
                 size = 4,
                 hjust = 0)
    }
    
    # add label repels only if we have top labels
    if (nrow(top_15_up_dn) > 0) {
      p <- p + geom_text_repel(
        data = top_15_up_dn,
        aes(x = log2FoldChange, y = Neg_log_p_val, label = Symbol),
        color = "black",
        arrow = arrow(ends = "last", type = "open"),
        max.overlaps = 50,
        na.rm = TRUE
      )
    }
    
    # final coordinate limits using finite ymax
    p + coord_cartesian(
      xlim = c(min(de_tab$log2FoldChange, na.rm = TRUE), max(de_tab$log2FoldChange, na.rm = TRUE)),
      ylim = c(0, ymax * 1.2)
    )
  })
  
  ma_plot <- reactive({
    
    de_tab <- res()
    
    validate(need(nrow(de_tab) > 0, "DE table is empty"))
    validate(need("log2FoldChange" %in% colnames(de_tab) && "baseMean" %in% colnames(de_tab),
                  "DE table must contain 'baseMean' and 'log2FoldChange' for MA plot."))
    
    if (!"gene" %in% colnames(de_tab)) de_tab$gene <- rownames(de_tab)
    if (!"Symbol" %in% colnames(de_tab)) de_tab$Symbol <- de_tab$gene
    
    de_tab$baseMean <- as.numeric(de_tab$baseMean)
    de_tab$log2FoldChange <- as.numeric(de_tab$log2FoldChange)
    de_tab$padj <- as.numeric(ifelse(is.na(de_tab$padj), 1, de_tab$padj))
    
    # log2 transform of baseMean
    de_tab$baseMean_log <- log2(de_tab$baseMean + 1)
    
    # Handle NA padj values
    de_tab$padj[is.na(de_tab$padj)] <- 1
    
    # Identify significant DE genes
    #lfc_cutoff <- as.numeric(input$log2FoldChange)
    #padj_cutoff <- as.numeric(input$padj)
    
    
    de_up <- subset(de_tab, log2FoldChange > input$log2FoldChange & padj < input$padj)
    de_dn <- subset(de_tab, log2FoldChange < - input$log2FoldChange & padj < input$padj)
    
    de_sig <- subset(de_tab, abs(log2FoldChange) > input$log2FoldChange & padj < input$padj);
    
    top_n_up_dn <- head(de_sig[order(de_sig$padj), ], 20)
    if (!"Symbol" %in% colnames(top_n_up_dn)) top_n_up_dn$Symbol <- top_n_up_dn$gene
    
    # Plot
    p <- ggplot(de_tab, aes(x = baseMean_log, y = log2FoldChange)) +
      geom_point(color = "grey70", size = 2) +
      theme_classic(base_size = 16) +
      ggtitle("MA Plot") +
      xlab("log2(BaseMean + 1)") +
      ylab("Log2 Fold Change")
    
    if (nrow(de_up) > 0) {
      p <- p + geom_point(data = de_up, aes(x = baseMean_log, y = log2FoldChange), color = "red", size = 3)
    }
    if (nrow(de_dn) > 0) {
      p <- p + geom_point(data = de_dn, aes(x = baseMean_log, y = log2FoldChange), color = "blue", size = 3)
    }
    
    if (nrow(top_n_up_dn) > 0) {
      p <- p + geom_text_repel(data = top_n_up_dn,
                               aes(x = baseMean_log, y = log2FoldChange, label = Symbol),
                               color = "black", arrow = arrow(ends = "last", type = "open"),
                               max.overlaps = 50, na.rm = TRUE)
    }
    
    # safe coords
    xlim <- range(de_tab$baseMean_log, na.rm = TRUE)
    ylim <- range(de_tab$log2FoldChange, na.rm = TRUE)
    #if (!all(is.finite(xlim))) xlim <- c(0, 1)
    #if (!all(is.finite(ylim))) ylim <- c(-5, 5)
    
    p + coord_cartesian(xlim = c(min(de_tab$baseMean_log), max(de_tab$baseMean_log)+1), ylim = c(min(de_tab$log2FoldChange)-1, max(de_tab$log2FoldChange)+1))
  })
  
  ## Volcano plot
  output$volcano <- renderPlot({
    volcano_plot()
  });
  
  ## MA plot
  output$maplot <- renderPlot({
    ma_plot()
  });
  
  ## Download volcano plot
  output$download_volcano <- downloadHandler(
    filename = function(){
      paste0("volcano_", Sys.Date(), ".png")
    },
    content = function(file){
      req(res())
      png(file, width = 12, height = 7, res = 300, units = "in")
      print(volcano_plot())
      dev.off()
    }
  )
  
  ## Download MA plot
  output$download_maplot <- downloadHandler(
    filename = function(){
      paste0("maplot_", Sys.Date(), ".png")
    },
    content = function(file){
      req(res())
      png(file, width = 12, height = 7, res = 300, units = "in")
      print(ma_plot())
      dev.off()
    }
  )
  
  # Download DE table
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("filtered_DE_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_res()$de, file, row.names = FALSE)
    }
  )
  

}
