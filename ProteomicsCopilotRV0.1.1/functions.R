#V0.1.1
library(shiny)
library(readr)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(plotly)
library(tidyverse)
library(gprofiler2)
library(corrplot)
library(pheatmap)
library(reticulate)
library(stringr)
library(xtable)
library(ggeasy)
library(zip)
library(DT)
library(rsvg)
library(png)
library(protr)
library(r3dmol)
library(UniprotR)
library(protti)
library(KSEAapp)
library(rstudioapi)
library(seqinr)
library(openxlsx)

extract_id <- function(sample_name) {
  str_extract(sample_name, "(?<=_)[0-9]+(?=\\.d)")
}

coverage_plot <- function(data, meta, id = TRUE, color_package = TRUE, header = TRUE, legend = TRUE) {
  data[data == 0] <- NA
  conditions = unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id == TRUE) {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  data_filtered[!is.na(data_filtered)] <- 1
  data_melted <- melt(data_filtered, variable.name = "Sample", value.name = "Value")
  data_annotated <- merge(data_melted, meta, by.x = "Sample", by.y = "new_sample")
  data_annotated$Sample <- factor(data_annotated$Sample, levels = meta$new_sample)
  data_annotated$condition <- factor(data_annotated$condition, levels = conditions)
  
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(conditions)) {
      plot_colors <- rep(plot_colors, length.out = length(conditions))
    }
  }
  
  plot_title <- ""
  y_label <- "Number"
  if (header) {
    if ("ProteinNames" %in% colnames(data)) {
      plot_title <- "Proteins per sample"
      y_label <- "Number of proteins"
    } else if ("PTM_Collapse_key" %in% colnames(data)) {
      plot_title <- "Phosphosites per sample"
      y_label <- "Number of phosphosites"
    }
  }
  
  p <- ggplot(data_annotated, aes(x = Sample, y = Value, fill = condition)) +
    stat_summary(fun = sum, geom = "bar", position = "dodge") +
    theme_minimal() +
    labs(title = plot_title, x = "Condition", y = y_label, fill = "Condition") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      legend.position = if (legend) "right" else "none"
    )
  
  if ("ProteinNames" %in% colnames(data)) {
    p <- p + geom_hline(yintercept = length(data$ProteinNames), linetype = "dashed", color = "red")
  } else if ("PTM_Collapse_key" %in% colnames(data)) {
    p <- p + geom_hline(yintercept = length(data$PTM_Collapse_key), linetype = "dashed", color = "red")
  }
  
  if (exists("plot_colors")) {
    p <- p + scale_fill_manual(values = plot_colors)
  }
  
  print(p)
}

missing_value_plot <- function(data, meta, bin = 0, header = TRUE, text = TRUE, text_size = 3.88) {
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  data_filtered[data_filtered == 0] <- NA
  na_count <- rowSums(is.na(data_filtered))
  
  if (bin > 0) {
    na_count[na_count > bin] <- paste0(">", bin)
  }
  
  miss_vals <- as.data.frame(table(na_count))
  miss_vals <- miss_vals[miss_vals$na_count != length(annotated_columns), ]
  
  if (bin > 0) {
    levels_vec <- c(as.character(0:bin), paste0(">", bin))
  } else {
    levels_vec <- as.character(0:max(as.numeric(as.character(miss_vals$na_count)), na.rm = TRUE))
  }
  
  miss_vals$na_count <- factor(miss_vals$na_count, levels = levels_vec)
  
  plot_title <- if ("ProteinNames" %in% colnames(data)) {
    "Missing Value Plot - Protein Level"
  } else if ("PTM_Collapse_key" %in% colnames(data)) {
    "Missing Value Plot - Phosphosite Level"
  } else {
    "Missing value plot"
  }
  
  p <- ggplot(data = miss_vals, aes(x = na_count, y = Freq)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(x = "Number of Missing Values", y = "Frequency")
  
  if (text) {
    p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = text_size)
  }
  
  if (header) {
    p <- p + ggtitle(plot_title)
  }
  
  p <- p +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  return(p)
}

histo_int <- function(data, meta, color_package = TRUE, header = TRUE, legend = TRUE) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  mean_intensities <- data.frame(condition = character(), mean_intensity = numeric(), stringsAsFactors = FALSE)
  for (condition in unique(meta$condition)) {
    columns <- meta$sample[meta$condition == condition]
    mean_intensity <- rowMeans(data_filtered[, columns, drop = FALSE], na.rm = TRUE)
    mean_intensity <- mean_intensity[is.finite(mean_intensity)]
    mean_intensities <- rbind(mean_intensities, data.frame(condition = condition, mean_intensity = mean_intensity))
  }
  
  p <- ggplot(mean_intensities, aes(x = mean_intensity, color = condition)) +
    geom_density() +
    theme_minimal() +
    labs(
      x = "log2 Intensity",
      y = "Density"
    )
  
  if (header) {
    p <- p + labs(title = "Distribution of measured intensity (log2)")
  }
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + labs(color = "Condition")
  }
  
  if (exists("plot_colors")) {
    p <- p + scale_color_manual(values = plot_colors) 
  }
  
  print(p)
}

boxplot_int <- function(data, meta, outliers = FALSE, header = TRUE, legend = TRUE, color_package = TRUE) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  mean_intensities <- data.frame(condition = character(), mean_intensity = numeric(), stringsAsFactors = FALSE)
  
  for (condition in unique(meta$condition)) {
    columns <- meta$sample[meta$condition == condition]
    mean_intensity <- rowMeans(data_filtered[, columns, drop = FALSE], na.rm = TRUE)
    mean_intensity <- mean_intensity[is.finite(mean_intensity)]
    mean_intensities <- rbind(mean_intensities, data.frame(condition = condition, mean_intensity = mean_intensity))
  }
  
  mean_intensities$condition <- factor(mean_intensities$condition, levels = unique(meta$condition))
  
  p <- ggplot(mean_intensities, aes(x = condition, y = mean_intensity, fill = condition)) +
    geom_boxplot(outlier.shape = if (outliers) 16 else NA) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    theme_minimal()
  
  if (header) {
    p <- p + labs(title = "Measured protein intensity values (log2)")
  }
  
  p <- p + labs(x = "Condition", y = "log2 Intensity", fill = "Condition")
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (exists("plot_colors")) {
    p <- p + scale_fill_manual(values = plot_colors)
  }
  
  return(p)
}


cov_plot <- function(data, meta, outliers = FALSE, color_package = TRUE, header = TRUE, legend = TRUE) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  cv_results <- data.frame(condition = character(), CV = numeric(), stringsAsFactors = FALSE)
  
  for (condition in unique(meta$condition)) {
    columns <- meta$sample[meta$condition == condition]
    condition_data <- data_filtered[, columns, drop = FALSE]
    
    means <- rowMeans(condition_data, na.rm = TRUE)
    sds <- apply(condition_data, 1, sd, na.rm = TRUE)
    cv <- (sds / means) * 100
    
    if (length(columns) >= 2) {
      cv_results <- rbind(cv_results, data.frame(condition = condition, CV = cv, stringsAsFactors = FALSE))
    }
  }
  
  cv_results <- cv_results[is.finite(cv_results$CV), ]
  cv_results$condition <- factor(cv_results$condition, levels = unique(meta$condition))
  
  p <- ggplot(cv_results, aes(x = condition, y = CV, fill = condition)) +
    geom_boxplot(outlier.shape = if (outliers) 16 else NA) +
    theme_minimal() +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    labs(
      title = if (header) "Coefficient of Variation" else NULL,
      x = "Condition",
      y = "Coefficient of Variation (%)"
    ) +
    theme(
      legend.position = if (legend) "right" else "none"
    )
  
  if (exists("plot_colors")) {
    p <- p + scale_fill_manual(values = plot_colors)
  }
  
  print(p)
}


boxplot_int_single <- function(data, meta, outliers = FALSE, id = TRUE, header = TRUE, legend = TRUE, color_package = TRUE) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  meta$id <- sapply(meta$sample, extract_id_or_number)
  if (id == TRUE) {
    meta <- meta %>%
      dplyr::group_by(condition) %>%
      dplyr::mutate(new_sample = paste0(condition, "_", dplyr::row_number(), "\n (", id, ")")) %>%
      dplyr::ungroup()
  } else {
    meta <- meta %>%
      dplyr::group_by(condition) %>%
      dplyr::mutate(new_sample = paste0(condition, "_", dplyr::row_number())) %>%
      dplyr::ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  intensities <- data.frame(sample = character(), intensity = numeric(), condition = character(), stringsAsFactors = FALSE)
  
  for (condition in unique(meta$condition)) {
    columns <- meta$new_sample[meta$condition == condition]
    for (column in columns) {
      intensity <- rowMeans(data_filtered[, column, drop = FALSE], na.rm = TRUE)
      intensity <- intensity[is.finite(intensity)]
      intensities <- rbind(intensities, data.frame(sample = column, intensity = intensity, condition = condition))
    }
  }
  
  intensities$sample <- factor(intensities$sample, levels = meta$new_sample)
  
  p <- ggplot(intensities, aes(x = sample, y = intensity, fill = condition)) +
    geom_boxplot(outlier.shape = if (outliers) 16 else NA) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (header) {
    p <- p + labs(title = "Measured protein intensity values (log2)")
  }
  
  p <- p + labs(x = "Sample", y = "log2 Intensity", fill = "Condition")
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  
  if (exists("plot_colors")) {
    p <- p + scale_fill_manual(values = plot_colors)
  }
  
  return(p)
}


pca_plot <- function(data, meta, color_package = TRUE, header = TRUE, legend = TRUE, dot_size = 3) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  meta$condition <- factor(meta$condition, levels = unique(meta$condition))
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  data_filtered <- data_filtered[complete.cases(data_filtered), ]
  
  transposed_expr <- t(data_filtered)
  
  tryCatch({
    zero_variance_columns <- which(apply(transposed_expr, 2, var) == 0)
    if (length(zero_variance_columns) > 0) {
      transposed_expr <- transposed_expr[, -zero_variance_columns, drop = FALSE]
    }
  }, error = function(e) {
    message("An error occurred while processing zero variance columns: ", e$message)
  })
  
  tryCatch({
    pca_result_transposed <- prcomp(transposed_expr, scale. = TRUE)
    explained_variance <- (pca_result_transposed$sdev^2) / sum(pca_result_transposed$sdev^2) * 100
    pca_scores <- as.data.frame(pca_result_transposed$x)
    pca_scores$sample <- rownames(pca_scores)
    
    pca_scores <- merge(pca_scores, meta, by.x = "sample", by.y = "sample")
    pca_scores$condition <- factor(pca_scores$condition, levels = levels(meta$condition))
    
    x_label <- sprintf("Principal Component 1 - %.2f%% variance", explained_variance[1])
    y_label <- sprintf("Principal Component 2 - %.2f%% variance", explained_variance[2])
    
    pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = dot_size, alpha = 0.7) +
      theme_minimal() +
      labs(x = x_label, y = y_label)
    
    if (header) {
      pca <- pca + labs(title = "PCA Plot")
    }
    
    pca <- pca + theme(legend.position = if (legend) "right" else "none")
    
    if (exists("plot_colors")) {
      pca <- pca + scale_color_manual(values = plot_colors)
    }
    
    return(pca)
  }, error = function(e) {
    message("An error occurred during PCA computation: ", e$message)
    return(NULL)
  })
}


abundance_plot <- function(data, meta, workflow = "Protein", color_package=TRUE) {
  data[data == 0] <- NA
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  annotated_columns <- meta$sample
  if (workflow == "Protein") {
    data_filtered <- data[, c("ProteinNames", annotated_columns), drop = FALSE]
    mean_intensities <- data.frame(Protein = data_filtered$ProteinNames)
  } else if (workflow == "Phosphosite") {
    data_filtered <- data[, c("PTM_Collapse_key", annotated_columns), drop = FALSE]
    mean_intensities <- data.frame(Protein = data_filtered$PTM_Collapse_key)
  }
  
  for (condition in unique(meta$condition)) {
    columns <- meta$sample[meta$condition == condition]
    if (workflow == "Protein") {
      condition_data <- data_filtered[, c("ProteinNames", columns), drop = FALSE]
    } else if (workflow == "Phosphosite") {
      condition_data <- data_filtered[, c("PTM_Collapse_key", columns), drop = FALSE]
    }
    condition_means <- rowMeans(condition_data[, columns, drop = FALSE], na.rm = TRUE)
    condition_means <- log10(condition_means + 1)
    mean_intensities[condition] <- condition_means
  }
  
  if (ncol(mean_intensities) > 2) {
    mean_intensities <- mean_intensities[rowSums(is.na(mean_intensities[, -1])) < ncol(mean_intensities) - 1, ]
  } else {
    mean_intensities <- mean_intensities[!is.na(mean_intensities[, 2]), ]
  }
  
  long_intensities <- mean_intensities %>%
    pivot_longer(-Protein, names_to = "Condition", values_to = "log10Intensity") %>%
    group_by(Condition) %>%
    mutate(Rank = rank(-log10Intensity, ties.method = "first")) %>%
    ungroup()
  
  unique_conditions <- unique(meta$condition)
  long_intensities$Condition <- factor(long_intensities$Condition, levels = unique_conditions)
  
  if (workflow == "Protein") {
    p <- ggplot(long_intensities, aes(x = Rank, y = log10Intensity, color = Condition)) +
      geom_line() +
      theme_minimal() +
      labs(title = "Abundance plot - all conditions",
           x = "Protein Rank",
           y = "log10 Protein Intensity",
           color = "Condition") +
      theme(legend.position = "right")
    if (exists("plot_colors")) {
      p <- p + scale_color_manual(values = plot_colors)
    }
  } else if (workflow == "Phosphosite") {
    p <- ggplot(long_intensities, aes(x = Rank, y = log10Intensity, color = Condition)) +
      geom_line() +
      theme_minimal() +
      labs(title = "Abundance plot - all conditions",
           x = "Phosphosite Rank",
           y = "log10 Phosphosite Intensity",
           color = "Condition") +
      theme(legend.position = "right")
    if (exists("plot_colors")) {
      p <- p + scale_color_manual(values = plot_colors)
    }
  }
  
  return(p)
}

interactive_abundance_plot <- function(data, meta, condition, workflow="Protein", search=NULL) {
  data[data == 0] <- NA
  annotated_columns <- meta$sample[meta$condition == condition]
  
  if (workflow=="Protein"){
    data_filtered <- data[, c("ProteinNames", annotated_columns), drop = FALSE]} else if (workflow=="Phosphosite"){
      data_filtered <- data[, c("PTM_Collapse_key", annotated_columns), drop = FALSE]}
  condition_means <- rowMeans(data_filtered[, annotated_columns, drop = FALSE], na.rm = TRUE)
  condition_means <- log10(condition_means + 1)
  
  if (workflow=="Protein"){
    mean_intensities <- data.frame(
      Feature = data_filtered$ProteinNames,
      log10Intensity = condition_means
    )} else if (workflow=="Phosphosite"){
      mean_intensities <- data.frame(
        Feature = data_filtered$PTM_Collapse_key,
        log10Intensity = condition_means
      )}
  mean_intensities <- mean_intensities[!is.na(mean_intensities$log10Intensity), ]
  
  mean_intensities <- mean_intensities %>%
    arrange(desc(log10Intensity)) %>%
    mutate(Rank = row_number())
  
  mean_intensities <- mean_intensities %>%
    mutate(Color = ifelse(Feature %in% search, "red", "blue"))
  
  p <- ggplot(mean_intensities, aes(x = Rank, y = log10Intensity, text = Feature)) +
    geom_point(aes(color = Color), size = 1) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(title = paste("Abundance plot - ", condition),
         x = paste(workflow, "Rank"),
         y = paste("log10", workflow, "Intensity")) +
    scale_color_identity() +
    theme(legend.position = "none")
  
  p_interactive <- ggplotly(p, tooltip = "text")
  return(p_interactive)
}

vline <- function(x) {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper", 
    x0 = x,
    x1 = x,
    line = list(color = "green", dash = "dot")
  )
}

hline <- function(y) {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = "green", dash="dot")
  )
}

filter_data <- function(data, meta, num, filterops = "per group") {
  meta$id <- sapply(meta$sample, extract_id)
  meta <- meta %>%
    group_by(condition) %>%
    mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
    ungroup()
  
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  if (filterops=="per group"){
    rows_to_keep <- apply(data_filtered, 1, function(row) {
      all(sapply(unique(meta$condition), function(condition) {
        condition_columns <- meta$sample[meta$condition == condition]
        sum(!is.na(row[condition_columns])) >= num
      }))
    })
  } else if (filterops =="in at least one group"){
    rows_to_keep <- apply(data_filtered, 1, function(row) {
      any(sapply(unique(meta$condition), function(condition) {
        condition_columns <- meta$sample[meta$condition == condition]
        sum(!is.na(row[condition_columns])) >= num
      }))
    })
  }
  
  filtered_data <- data[rows_to_keep, , drop = FALSE]
  return(filtered_data)
}

different_genes <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1) {
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  data_filtered <- data[, c("ProteinNames", annotated_columns1, annotated_columns2), drop = FALSE]
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, -1])) < (length(annotated_columns1) + length(annotated_columns2) - 2), ]
  
  log2fc <- apply(data_filtered[, -1], 1, function(row) {
    mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
    mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
    mean2 - mean1
  })
  
  pvals <- apply(data_filtered[, -1], 1, function(row) {
    values1 <- na.omit(as.numeric(row[annotated_columns1]))
    values2 <- na.omit(as.numeric(row[annotated_columns2]))
    if (length(values1) < 2 || length(values2) < 2) {
      return(NA)
    }
    t.test(values1, values2, var.equal = TRUE)$p.value
  })
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  data_filtered <- data_filtered[valid_idx, ]
  
  adj_pvals <- p.adjust(pvals, method = "BH")
  
  result_data <- data.frame(
    Protein = data_filtered$`ProteinNames`,
    log2FC = log2fc,
    pval = pvals,
    adj_pval = adj_pvals
  )
  
  if ("GeneNames" %in% colnames(data)) {
    result_data$Gene <- data$GeneNames[valid_idx]
  } else {
    result_data$Gene <- NA
  }
  
  sig_pval <- in_pval
  sig_log2fc <- in_log2fc
  
  upregulated_genes <- result_data %>%
    filter(adj_pval < sig_pval & log2FC > sig_log2fc) %>%
    pull(ifelse("Gene" %in% colnames(result_data), Gene, Protein))
  
  downregulated_genes <- result_data %>%
    filter(adj_pval < sig_pval & log2FC < -sig_log2fc) %>%
    pull(ifelse("Gene" %in% colnames(result_data), Gene, Protein))
  
  return(list(Upregulated = upregulated_genes, Downregulated = downregulated_genes))
}

enrichment_analysis <- function(gene_list, top_n = 10, min_num=20, max_num=300) {
  gene_list <- na.omit(gene_list)
  res <- gost(query = gene_list, organism = "hsapiens", sources = c("GO:CC", "GO:BP", "GO:MF"), user_threshold = 0.05, correction_method="g_SCS")

  if (is.null(res) || is.null(res$result)) {
    stop("No results returned from gost.")
  }
  
  res_df <- res$result %>%
    mutate(
      hits = intersection_size,
      hitsPerc = intersection_size / term_size * 100,
      pval = p_value
    ) %>%
    arrange(pval) %>%
    filter(term_size >= min_num & term_size <= max_num)
  
  if (nrow(res_df) == 0) {
    return(NULL)
  }
  
  ggplot(res_df %>% top_n(-top_n, wt = pval), 
         aes(x = hitsPerc, 
             y = reorder(term_name, -pval),
             colour = pval, 
             size = hits)) +
    geom_point() +
    expand_limits(x = 0) +
    labs(x = "Hits (%)", 
         y = "GO term", 
         colour = "p value", 
         size = "Count") +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") + 
    ggtitle("Gene Set Enrichment Analysis") 
}

corr_plot <- function(data, meta, method=FALSE, id=TRUE, full_range=FALSE) {
  conditions <- unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id==TRUE){
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  correlation_matrix <- cor(data_filtered, use = "pairwise.complete.obs", method = "pearson")
  hc <- hclust(as.dist(1 - correlation_matrix), method = "complete")
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))  
  
  col_range <- if (full_range) c(-1, 1) else range(correlation_matrix, na.rm = TRUE)
  
  if (method==FALSE){
    corrplot(correlation_matrix[hc$order, hc$order], 
             method = "color",          
             type = "full",            
             order = "original",         
             tl.col = "black",          
             tl.srt = 45,               
             tl.cex = 0.6,              
             number.cex = 0.6,          
             tl.pos = "lt",
             is.corr = FALSE,         
             col.lim = col_range)     
  } else {
    corrplot(correlation_matrix[hc$order, hc$order], 
             method = "ellipse",          
             type = "full",            
             order = "original",         
             tl.col = "black",          
             tl.srt = 45,               
             tl.cex = 0.6,              
             number.cex = 0.6,          
             tl.pos = "lt",
             is.corr = FALSE,         
             col.lim = col_range)     
  }
}

rename_cols <- function(df) {
  if ("PG.ProteinNames" %in% names(df)) {
    names(df)[names(df) == "PG.ProteinNames"] <- "ProteinNames"
  }
  if ("PG.ProteinGroups" %in% names(df)) {
    names(df)[names(df) == "PG.ProteinGroups"] <- "ProteinNames"
  }
  if ("Protein.IDs" %in% names(df)) {
    names(df)[names(df) == "Protein.IDs"] <- "ProteinNames"
  }
  if ("PG.Genes" %in% names(df)) {
    names(df)[names(df) == "PG.Genes"] <- "GeneNames"
  }
  if ("Protein.Names" %in% names(df)) {
    names(df)[names(df) == "Protein.Names"] <- "ProteinNames"
  }
  if ("Protein IDs" %in% names(df)) {
    names(df)[names(df) == "Protein IDs"] <- "ProteinNames"
  }
  if ("Genes" %in% names(df)) {
    names(df)[names(df) == "Genes"] <- "GeneNames"
  }
  if ("Gene names" %in% names(df)) {
    names(df)[names(df) == "Gene names"] <- "GeneNames"
  }
  if ("PEP.StrippedSequence" %in% names(df)) {
    names(df)[names(df) == "PEP.StrippedSequence"] <- "Stripped.Sequence"
  }
  if ("Sequence" %in% names(df)) {
    names(df)[names(df) == "Sequence"] <- "Stripped.Sequence"
  }
  if ("EG.ModifiedSequence" %in% names(df)) {
    names(df)[names(df) == "EG.ModifiedSequence"] <- "Modified.Sequence"
  }
  if ("R.FileName" %in% names(df)) {
    names(df)[names(df) == "R.FileName"] <- "File.Name"
  }
  
  matching_col <- grep("EG\\.TotalQuantity", names(df), value = TRUE)
  if (length(matching_col) > 0) {
    names(df)[names(df) == matching_col] <- "Precursor.Quantity"
  }
  
  if ("EG.iRTPredicted" %in% names(df)) {
    names(df)[names(df) == "EG.iRTPredicted"] <- "Predicted.RT"
  }
  if ("EG.iRTEmpirical" %in% names(df)) {
    names(df)[names(df) == "EG.iRTEmpirical"] <- "RT"
  }
  
  if ("Raw file" %in% names(df)) {
    df$Intensity <- as.numeric(df$Intensity)
    
    df <- aggregate(
      Intensity ~ GeneNames + `Protein names` + `Raw file`,
      data = df,
      FUN = function(x) sum(x, na.rm = TRUE)
    )
    
    df_wide <- reshape(
      df,
      idvar = c("GeneNames", "Protein names"),
      timevar = "Raw file",
      direction = "wide"
    )
    
    names(df_wide)[names(df_wide) == "Protein names"] <- "ProteinNames"
    df <- df_wide
  }
  
  if ("Reverse" %in% names(df)) {
    df <- df[!grepl("\\+", df$Reverse), ]
  }
  if ("Only identified by site" %in% names(df)) {
    df <- df[!grepl("\\+", df[["Only identified by site"]]), ]
  }
  
  return(df)
}


heatmap_plot <- function(data, meta, id=TRUE) {
  conditions = unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id)
  if (id==TRUE){
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  }else{
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered = data[, annotated_columns, drop = FALSE]
  data_col_names <-colnames(data_filtered)
  
  data_for_clustering <- t(apply(data_filtered, 1, function(row) {
    if (any(is.na(row))) {
      row[is.na(row)] <- mean(row, na.rm = TRUE)
    }
    return(row)
  }))
  
  if ("row" %in% "row") {
    data_for_clustering <- t(apply(data_for_clustering, 1, function(row) {
      scaled_row <- scale(row)
      if (any(is.nan(scaled_row)) || any(is.infinite(scaled_row))) {
        return(rep(0, length(scaled_row)))  # Replace problematic rows with zeros
      }
      return(scaled_row)
    }))
  }
  
  
  if ("row" %in% "row") {
    data_filtered <- t(apply(data_filtered, 1, function(row) {
      scaled_row <- scale(row)
      if (any(is.nan(scaled_row)) || any(is.infinite(scaled_row))) {
        return(rep(0, length(scaled_row)))  # Replace problematic rows with zeros
      }
      return(scaled_row)
    }))
    colnames(data_filtered) <- data_col_names
  }
  
  annotation_col <- data.frame(condition = meta$condition[match(data_col_names, meta$new_sample)])
  rownames(annotation_col) <- data_col_names
  
  palette <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
  na_color <- "grey"
  
  row_dist <- dist(data_for_clustering)
  row_clust <- hclust(row_dist)
  col_dist <- dist(t(data_for_clustering))
  col_clust <- hclust(col_dist)
  
  p=pheatmap(data_filtered,
             cluster_rows = row_clust, 
             cluster_cols = col_clust, 
             na_col = na_color, 
             color = palette, 
             annotation_col = annotation_col, 
             annotation_names_col = FALSE)
  return(p)
}

missing_value_plot_prec <- function(data2, meta, bin = 0, header = TRUE, text = TRUE, text_size = 3.88) {
  data2 <- data2 %>%
    pivot_wider(
      id_cols = File.Name,
      names_from = Precursor.Id,
      values_from = Precursor.Quantity,
      values_fill = NA
    )
  
  data2 <- as.data.frame(t(data2))
  colnames(data2) <- data2[1, ]
  data2 <- data2[-1, , drop = FALSE]
  
  colnames(data2) <- sapply(colnames(data2), extract_id)
  meta$sample <- sapply(meta$sample, extract_id)
  
  annotated_columns <- meta$sample
  data_filtered <- data2[, annotated_columns, drop = FALSE]
  
  na_count <- rowSums(is.na(data_filtered))
  
  if (bin > 0) {
    na_count[na_count > bin] <- paste0(">", bin)
  }
  
  miss_vals <- as.data.frame(table(na_count))
  miss_vals <- miss_vals[miss_vals$na_count != length(annotated_columns), ]
  
  if (bin > 0) {
    levels_vec <- c(as.character(0:bin), paste0(">", bin))
  } else {
    levels_vec <- as.character(0:max(as.numeric(as.character(miss_vals$na_count)), na.rm = TRUE))
  }
  miss_vals$na_count <- factor(miss_vals$na_count, levels = levels_vec)
  
  plot_title <- if (header) "Missing Value Plot - Precursor Level" else NULL
  
  p <- ggplot(data = miss_vals, aes(x = na_count, y = Freq)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(x = "Number of Missing Values", y = "Frequency", title = plot_title) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  if (text) {
    p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = text_size)
  }
  
  print(p)
}

missing_value_plot_pep <- function(data2, meta, bin = 0, header = TRUE, text = TRUE, text_size = 3.88) {
  plot_missing <- function(data_matrix, annotated_columns, plot_title) {
    data_matrix[data_matrix == 0] <- NA
    na_count <- rowSums(is.na(data_matrix))
    
    if (bin > 0) {
      na_count[na_count > bin] <- paste0(">", bin)
    }
    
    miss_vals <- as.data.frame(table(na_count))
    miss_vals <- miss_vals[miss_vals$na_count != length(annotated_columns), ]
    
    if (bin > 0) {
      levels_vec <- c(as.character(0:bin), paste0(">", bin))
    } else {
      levels_vec <- as.character(0:max(as.numeric(as.character(miss_vals$na_count)), na.rm = TRUE))
    }
    
    miss_vals$na_count <- factor(miss_vals$na_count, levels = levels_vec)
    
    p <- ggplot(data = miss_vals, aes(x = na_count, y = Freq)) +
      geom_bar(stat = "identity", fill = "blue") +
      labs(x = "Number of Missing Values", y = "Frequency", title = plot_title) +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    
    if (text) {
      p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = text_size)
    }
    
    print(p)
  }
  
  if ("File.Name" %in% names(data2)) {
    data2_names <- unique(data2$File.Name)
    for (name in data2_names) {
      meta$sample <- ifelse(meta$sample == name, name, meta$sample)
    }
    
    data2 <- data2 %>%
      pivot_wider(
        id_cols = File.Name,
        names_from = Stripped.Sequence,
        values_from = Precursor.Quantity,
        values_fill = NA,
        values_fn = list(Precursor.Quantity = max)
      )
    
    data2 <- as.data.frame(t(data2))
    colnames(data2) <- data2[1, ]
    data2 <- data2[-1, , drop = FALSE]
    colnames(data2) <- sapply(colnames(data2), extract_id_or_number)
    
    meta$sample <- sapply(meta$sample, extract_id_or_number)
    annotated_columns <- meta$sample
    data_filtered <- data2[, annotated_columns, drop = FALSE]
    
    plot_title <- if (header) "Missing Value Plot - Peptide Level" else NULL
    plot_missing(data_filtered, annotated_columns, plot_title)
    
  } else {
    annotated_columns <- meta$sample
    data2_filtered <- data2[, annotated_columns, drop = FALSE]
    plot_title <- if (header) "Missing value plot - peptide level" else NULL
    plot_missing(data2_filtered, annotated_columns, plot_title)
  }
}

coverage_plot_pep <- function(data2, meta, id = TRUE, header = TRUE, legend = TRUE, color_package = TRUE) {
  if ("File.Name" %in% names(data2)) {
    data2_names <- unique(data2$File.Name)
    for (name in data2_names) {
      meta$sample <- ifelse(meta$sample == name, name, meta$sample)
    }
    data2 <- data2 %>%
      pivot_wider(
        id_cols = File.Name,
        names_from = Stripped.Sequence,
        values_from = Precursor.Quantity,
        values_fill = list(Precursor.Quantity = NA),
        values_fn = list(Precursor.Quantity = max)
      )
    
    data2 <- as.data.frame(t(data2))
    colnames(data2) <- data2[1, ]
    data2 <- data2[-1, , drop = FALSE]
    colnames(data2) <- sapply(colnames(data2), extract_id_or_number)
    
    conditions <- unique(meta$condition)
    meta$sample <- as.character(meta$sample)
    meta$id <- sapply(meta$sample, extract_id_or_number)
    
    if (id == TRUE) {
      meta <- meta %>%
        group_by(condition) %>%
        mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
        ungroup()
    } else {
      meta <- meta %>%
        group_by(condition) %>%
        mutate(new_sample = paste0(condition, "_", row_number())) %>%
        ungroup()
    }
    
    rename_vector <- setNames(meta$new_sample, meta$id)
    colnames(data2) <- ifelse(colnames(data2) %in% meta$id, rename_vector[colnames(data2)], colnames(data2))
    annotated_columns <- meta$new_sample
    data_filtered <- data2[, annotated_columns, drop = FALSE]
    
    data_filtered[data_filtered == 0] <- NA
    data_filtered[!is.na(data_filtered)] <- 1
    data_filtered[is.na(data_filtered)] <- 0
    data_filtered[] <- lapply(data_filtered, as.numeric)
    
    column_sums <- colSums(data_filtered, na.rm = TRUE)
    sum_df <- data.frame(Sample = names(column_sums), Sum = column_sums)
    sum_df <- merge(sum_df, meta[, c("new_sample", "condition")], by.x = "Sample", by.y = "new_sample", all.x = TRUE)
    sum_df$Sample <- factor(sum_df$Sample, levels = meta$new_sample)
    
    coverage_plot <- ggplot(sum_df, aes(x = Sample, y = Sum, fill = condition)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_hline(yintercept = nrow(data_filtered), linetype = "dashed", color = "red")
    
    if (header) {
      coverage_plot <- coverage_plot +
        labs(title = "Peptides per sample", x = "Samples", y = "Sum of Precursor Quantities")
    } else {
      coverage_plot <- coverage_plot +
        labs(title = NULL, x = NULL, y = NULL)
    }
    
    if (!legend) {
      coverage_plot <- coverage_plot + theme(legend.position = "none")
    } else {
      coverage_plot <- coverage_plot + labs(fill = "Condition")
    }
    
    if (color_package && exists("plot_colors")) {
      if (length(plot_colors) < length(conditions)) {
        plot_colors <- rep(plot_colors, length.out = length(conditions))
      }
      coverage_plot <- coverage_plot + scale_fill_manual(values = plot_colors)
    }
    
    print(coverage_plot)
  } else {
    data2[data2 == 0] <- NA
    conditions = unique(meta$condition)
    meta$sample <- as.character(meta$sample)
    meta$id <- sapply(meta$sample, extract_id_or_number)
    
    if (id == TRUE) {
      meta <- meta %>%
        group_by(condition) %>%
        mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
        ungroup()
    } else {
      meta <- meta %>%
        group_by(condition) %>%
        mutate(new_sample = paste0(condition, "_", row_number())) %>%
        ungroup()
    }
    
    rename_vector <- setNames(meta$new_sample, meta$sample)
    colnames(data2) <- ifelse(colnames(data2) %in% meta$sample, rename_vector[colnames(data2)], colnames(data2))
    annotated_columns <- meta$new_sample
    data2_filtered = data2[, annotated_columns, drop = FALSE]
    
    data2_filtered[!is.na(data2_filtered)] <- 1
    data2_melted = melt(data2_filtered, variable.name = "Sample", value.name = "Value")
    data2_annotated = merge(data2_melted, meta, by.x = "Sample", by.y = "new_sample")
    data2_annotated$Sample <- factor(data2_annotated$Sample, levels = meta$new_sample)
    data2_annotated$condition <- factor(data2_annotated$condition, levels = conditions)
    
    p <- ggplot(data2_annotated, aes(x = Sample, y = Value, fill = condition)) +
      stat_summary(fun = sum, geom = "bar", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_hline(yintercept = length(data2$Stripped.Sequence), linetype = "dashed", color = "red")
    
    if (header) {
      p <- p + labs(title = "Peptides per sample", x = "Condition", y = "Number of peptides")
    } else {
      p <- p + labs(title = NULL, x = NULL, y = NULL)
    }
    
    if (!legend) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + labs(fill = "Condition")
    }
    
    if (color_package && exists("plot_colors")) {
      if (length(plot_colors) < length(conditions)) {
        plot_colors <- rep(plot_colors, length.out = length(conditions))
      }
      p <- p + scale_fill_manual(values = plot_colors)
    }
    
    print(p)
  }
}


RTvspredRT_plot <- function(data2, meta, method = "Hexbin Plot", add_line = FALSE, bins = 1000, header = TRUE) {
  x <- data2$Predicted.RT
  y <- data2$RT
  
  if (method == "Scatter Plot") {
    if (nrow(data2) > 100000) {
      set.seed(69)
      sample_indices <- sample(1:nrow(data2), 100000)
      data2 <- data2[sample_indices, ]
    }
    
    plot <- ggplot(data2, aes(x = Predicted.RT, y = RT)) +
      geom_point(alpha = 0.2, size = 0.1) +
      theme_minimal() +
      labs(x = "Predicted RT", y = "RT")
    
  } else if (method == "Density Plot") {
    plot <- ggplot(data2, aes(x = Predicted.RT, y = RT)) +
      stat_density_2d() +
      labs(x = "Predicted RT", y = "Actual RT") +
      theme_minimal()
    
  } else if (method == "Hexbin Plot") {
    plot <- ggplot(data2, aes(x = Predicted.RT, y = RT)) +
      geom_hex(bins = bins) +
      labs(x = "Predicted RT", y = "Actual RT") +
      theme_minimal()
  }
  
  if (add_line) {
    plot <- plot + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
  }
  
  if (header) {
    plot <- plot + ggtitle("Retention Time Plot")
  }
  
  print(plot)
}

clear_all_plots <- function() {
  while (dev.cur() > 1) {
    dev.off()
  }
}

modification_plot <- function(data2, meta, id = TRUE, header = TRUE) {
  data2_names <- unique(data2$File.Name)
  for (name in data2_names) {
    meta$sample <- ifelse(meta$sample == name, name, meta$sample)
  }
  
  data2 <- data2 %>%
    pivot_wider(
      id_cols = File.Name,
      names_from = Modified.Sequence,
      values_from = Precursor.Quantity,
      values_fill = list(Precursor.Quantity = NA),
      values_fn = list(Precursor.Quantity = max)
    )
  
  data2 <- as.data.frame(t(data2))
  colnames(data2) <- data2[1, ]
  data2 <- data2[-1, , drop = FALSE]
  colnames(data2) <- sapply(colnames(data2), extract_id_or_number)
  
  conditions <- unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id) {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$id)
  colnames(data2) <- ifelse(colnames(data2) %in% meta$id, rename_vector[colnames(data2)], colnames(data2))
  annotated_columns <- intersect(meta$new_sample, colnames(data2))
  data_filtered <- data2[, annotated_columns, drop = FALSE]
  data_filtered <- cbind(Modified.Sequence = rownames(data_filtered), data_filtered)
  rownames(data_filtered) <- 1:nrow(data_filtered)
  
  get_mod_count <- function(df, pattern) {
    df_mod <- df %>%
      filter(str_detect(Modified.Sequence, pattern)) %>%
      mutate(across(-Modified.Sequence, ~ ifelse(is.na(.), 0, 1)))
    colSums(df_mod[, -1])
  }
  
  col_sums_carb <- get_mod_count(data_filtered, "UniMod:4|Carbamidomethyl")
  col_sums_oxi <- get_mod_count(data_filtered, "UniMod:35|Oxidation")
  col_sums_ace <- get_mod_count(data_filtered, "UniMod:1|Acetyl")
  
  plot_data <- data.frame(
    Sample = rep(names(col_sums_carb), 3),
    Count = c(col_sums_carb, col_sums_oxi, col_sums_ace),
    Modification = rep(c("Carbamylation", "Oxidation", "Acetylation"), each = length(col_sums_carb))
  ) %>%
    filter(Count > 0)
  
  plot_data$Sample <- factor(plot_data$Sample, levels = meta$new_sample)
  
  plot <- ggplot(plot_data, aes(x = Sample, y = Count, fill = Modification)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey")
    ) +
    scale_fill_manual(
      values = c("Carbamylation" = "blue", "Oxidation" = "red", "Acetylation" = "green"),
      name = "Modification"
    ) +
    labs(
      x = "Sample",
      y = "Number of modified peptides",
      title = if (header) "Modifications per sample" else NULL
    )
  
  print(plot)
}

extract_id_or_number <- function(name) {
  id <- extract_id(name)
  if (is.na(id)) {
    number <- str_extract(name, "(?<=_)[0-9]+$")
    return(number)
  }
  return(id)
}

missed_cl_plot <- function(data2, meta, id = TRUE, text = TRUE, text_size = 3.88, header = TRUE) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  
  data2_names <- unique(data2$File.Name)
  for (name in data2_names) {
    meta$sample <- ifelse(meta$sample == name, name, meta$sample)
  }
  
  data2 <- data2 %>%
    pivot_wider(
      id_cols = File.Name,
      names_from = Stripped.Sequence,
      values_from = Precursor.Quantity,
      values_fill = list(Precursor.Quantity = NA),
      values_fn = list(Precursor.Quantity = max)
    )
  
  data2 <- as.data.frame(t(data2))
  colnames(data2) <- data2[1, ]
  data2 <- data2[-1, , drop = FALSE]
  colnames(data2) <- sapply(colnames(data2), extract_id_or_number)
  
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id == TRUE) {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$id)
  colnames(data2) <- ifelse(colnames(data2) %in% meta$id, rename_vector[colnames(data2)], colnames(data2))
  annotated_columns <- intersect(meta$new_sample, colnames(data2))
  data_filtered <- data2[, annotated_columns, drop = FALSE]
  data_filtered <- cbind(Stripped.Sequence = rownames(data_filtered), data_filtered)
  rownames(data_filtered) <- 1:nrow(data_filtered)
  
  count_all_RK <- function(sequence) {
    count_R <- str_count(sequence, "R")
    count_K <- str_count(sequence, "K")
    count <- count_R + count_K
    count <- pmax(0, count - 1)
    return(count)
  }
  
  rk_counts <- sapply(data_filtered$Stripped.Sequence, count_all_RK)
  
  for (col in colnames(data_filtered)[-1]) {
    data_filtered[[col]] <- ifelse(!is.na(data_filtered[[col]]), rk_counts, NA)
  }
  
  plot_data <- data_filtered %>%
    select(-Stripped.Sequence) %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Count") %>%
    filter(!is.na(Count)) %>%
    group_by(Sample, Count) %>%
    summarise(Occurrences = n(), .groups = 'drop') %>%
    ungroup() %>%
    group_by(Sample) %>%
    mutate(Total = sum(Occurrences), Percentage = (Occurrences / Total) * 100) %>%
    ungroup()
  
  plot_data$Sample <- factor(plot_data$Sample, levels = meta$new_sample)
  plot_data$Count <- factor(plot_data$Count, levels = sort(unique(plot_data$Count), decreasing = TRUE))
  
  missed_cleavage_labels <- plot_data %>%
    filter(Count == 1) %>%
    mutate(Label = sprintf("%.1f%%", Percentage))
  
  max_y_value <- max(plot_data$Percentage, na.rm = TRUE) * 1.05
  missed_cleavage_labels <- missed_cleavage_labels %>%
    mutate(Label_Y = max_y_value * 0.99)
  
  p <- ggplot(plot_data, aes(x = Sample, y = Percentage, fill = factor(Count))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = if (header) "Missed cleavages per sample" else NULL,
      x = "Sample",
      y = "Percentage of peptides with missed cleavage (%)",
      fill = "Number"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (text) {
    p <- p + geom_text(
      data = missed_cleavage_labels,
      aes(label = Label, y = Label_Y),
      vjust = 0,
      colour = "white",
      size = text_size
    )
  }
  
  return(p)
}

compare_prot_box <- function(data, meta, conditions, inputs, workflow="Protein") {
  condition_samples <- lapply(conditions, function(cond) {
    samples <- meta$sample[meta$condition == cond]
    if (length(samples) == 0) {
      warning(paste("No samples found for condition:", cond))
    }
    return(samples)
  })
  names(condition_samples) <- conditions
  samples_to_include <- unlist(condition_samples)
  
  if (workflow=="Protein"){
    data_filtered <- data[, c("ProteinNames", samples_to_include), drop = FALSE]
    data_protein <- data_filtered[data_filtered$ProteinNames == inputs, ]
    data_melted <- reshape2::melt(data_protein, id.vars = "ProteinNames", variable.name = "Sample", value.name = "Expression")
  } else if (workflow=="Phosphosite"){
    data_filtered <- data[, c("PTM_Collapse_key", samples_to_include), drop = FALSE]
    data_phos <- data_filtered[data_filtered$PTM_Collapse_key == inputs, ]
    data_melted <- reshape2::melt(data_phos, id.vars = "PTM_Collapse_key", variable.name = "Sample", value.name = "Expression")
  }
  data_melted$Condition <- NA
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    samples <- condition_samples[[i]]
    
    if (length(samples) > 0) {
      data_melted$Condition[data_melted$Sample %in% samples] <- condition
    }
  }
  data_melted <- na.omit(data_melted)
  
  valid_conditions <- conditions[sapply(conditions, function(cond) {
    sum(data_melted$Condition == cond) >= 3
  })]
  
  data_melted <- data_melted[data_melted$Condition %in% valid_conditions, ]
  
  if (length(valid_conditions) < 2) {
    stop("At least two conditions with at least three valid values each are required for comparison.")
  }
  
  data_melted$Condition <- factor(data_melted$Condition, levels = valid_conditions)
  
  if (workflow=="Protein"){
    boxplot <- ggplot(data_melted, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot() +
      stat_boxplot(geom = 'errorbar', width = 0.3) +
      labs(title = paste("Protein Expression -", inputs),
           x = "Condition",
           y = "Log2 intensity") + 
      theme_minimal() +
      theme(legend.position = "none")
  } else if (workflow=="Phosphosite"){
    boxplot <- ggplot(data_melted, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot() +
      stat_boxplot(geom = 'errorbar', width = 0.3) +
      labs(title = paste("Phosphosite Expression -", inputs),
           x = "Condition",
           y = "Log2 intensity") + 
      theme_minimal() +
      theme(legend.position = "none")
  }  
  
  return(boxplot)
}

compare_prot_line <- function(data, meta, conditions, inputs, id=TRUE, workflow="Protein") {
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  if (id == TRUE) {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  relevant_samples <- meta$new_sample
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  
  if (workflow=="Protein"){
    data_filtered <- data %>%
      filter(ProteinNames %in% inputs) %>%
      select(ProteinNames, one_of(relevant_samples))
    data_melted <- reshape2::melt(data_filtered, id.vars = "ProteinNames", variable.name = "Sample", value.name = "Value")
  } else if (workflow=="Phosphosite"){
    data_filtered <- data %>%
      filter(PTM_Collapse_key %in% inputs) %>%
      select(PTM_Collapse_key, one_of(relevant_samples))
    data_melted <- reshape2::melt(data_filtered, id.vars = "PTM_Collapse_key", variable.name = "Sample", value.name = "Value")
  }
  
  data_melted <- data_melted %>%
    left_join(meta %>% select(new_sample, condition), by = c("Sample" = "new_sample")) %>%
    na.omit() 
  
  ordered_samples <- unlist(lapply(conditions, function(cond) {
    meta$new_sample[meta$condition == cond] 
  }))
  
  data_melted$Sample <- factor(data_melted$Sample, levels = ordered_samples)
  if (workflow=="Protein"){
    data_melted$ProteinNames <- sapply(data_melted$ProteinNames, function(x) {
      strsplit(x, ";")[[1]][1]
    })}
  
  if (workflow=="Protein"){
    line_plot <- ggplot(data_melted, aes(x = Sample, y = Value, group = ProteinNames, color = ProteinNames)) +
      geom_line() +
      geom_point() +
      labs(title = "Protein Expression Across Samples",
           x = "Sample",
           y = "Log2 intensity") +
      guides(color = guide_legend(title = "Protein")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "right")
  } else if (workflow=="Phosphosite"){
    line_plot <- ggplot(data_melted, aes(x = Sample, y = Value, group = PTM_Collapse_key, color = PTM_Collapse_key)) +
      geom_line() +
      geom_point() +
      labs(title = "Phosphosite Expression Across Samples",
           x = "Sample",
           y = "Log2 intensity") +
      guides(color = guide_legend(title = "Phosphosite")) +
    theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "right")
  }
  
  return(line_plot)
}

RpeptideCollapse <- function(data, cutoff=0.75, collapse_level="PG"){
  colnames(data)[colnames(data) == "EG.TotalQuantity..Settings."] <- "EG.TotalQuantity (Settings)"
  collapse_data = as.data.frame(peptideCollapse(data, cutoff = cutoff, collapse_level = collapse_level))
  collapse_values = c()
  print("Advanced reordering!")
  for (z in 1:ncol(collapse_data)){
    for (i in 0:(nrow(collapse_data)-1)){
      value = collapse_data[i,z]
      collapse_values = c(collapse_values, value)
    }
  }
  collapse_summary = lapply(collapse_values, as.character)
  collapse_summary = data.frame(matrix(unlist(collapse_summary), nrow = nrow(collapse_data), byrow = F))
  colnames(collapse_summary) = colnames(collapse_data)
  collapse_summary <- collapse_summary %>%
    mutate(across(
      .cols = 1:(ncol(collapse_summary) - 5), 
      .fns = ~ as.numeric(ifelse(. == "nan", NA, .))
    ))
  print("Done")
  return(collapse_summary)
}

log2_transform_data <- function(data, meta) {
  annotated_columns = meta$sample
  data_filtered = data[, annotated_columns, drop = FALSE]
  data_filtered[data_filtered == 0] <- NA
  log2_data = log2(data_filtered)
  remaining_columns = setdiff(colnames(data), annotated_columns)
  combined_data = cbind(log2_data, data[, remaining_columns, drop = FALSE])
  return(combined_data)
}

inverseof_log2_transform_data <- function(data, meta) {
  annotated_columns = meta$sample
  log2_data = data[, annotated_columns, drop = FALSE]
  original_data = 2^log2_data
  remaining_columns = setdiff(colnames(data), annotated_columns)
  combined_data = cbind(original_data, data[, remaining_columns, drop = FALSE])
  return(combined_data)
}

impute_values <- function(data, meta, q=0.01, adj_std=1, ret=0, sample_wise=F, seed=69){
  set.seed(seed)
  data[data == 0] <- NA
  annotated_columns <- meta$sample
  all_columns <- colnames(data)
  remaining_columns <- setdiff(all_columns, annotated_columns)
  data_filtered <- data[, annotated_columns, drop = FALSE]
  combined_data <- unlist(data_filtered)
  missing_value_mask <- is.na(data_filtered)
  mask_df <- as.data.frame(missing_value_mask)
  
  data_without_na <- data_filtered[!rowSums(missing_value_mask), , drop = FALSE]
  data_with_na <- data_filtered[rowSums(missing_value_mask) > 0, , drop = FALSE]
  combined_data_without_na <- unlist(data_without_na)
  combined_data_with_na <- unlist(data_with_na)
  
  plot_data <- data.frame(
    value = c(combined_data_without_na, combined_data_with_na),
    category = rep(c("Without Missing values", "With Missing Values"), c(length(combined_data_without_na), length(combined_data_with_na)))
  )
  
  hist_plot1 <- ggplot(plot_data, aes(x = value, fill = category)) +
    geom_histogram(aes(y = ..density..), position = "identity", binwidth = 0.25, alpha = 0.5) +
    labs(title = "Overall distribution of data with and without missing values",
         x = "log2 Intensity",
         y = "Density") +
    theme_minimal()
  
  mean_val <- mean(combined_data, na.rm = TRUE)
  sd_val <- sd(combined_data, na.rm = TRUE)
  quantile_val <- quantile(combined_data, probs = q, na.rm = TRUE)
  
  hist_plot2 <- ggplot(data.frame(value = combined_data), aes(x = value)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.25, fill = "blue", alpha = 0.5) +
    stat_function(fun = dnorm, args = list(mean = mean_val, sd = sd_val), 
                  color = "red", size = 1) +
    labs(title = "Overall data distribution and norm fit",
         x = "log2 Intensity",
         y = "Density") +
    theme_minimal()
  
  if (sample_wise == T){
    for (col in annotated_columns) {
      col_mean <- mean(data_filtered[[col]], na.rm = TRUE)
      col_sd <- sd(data_filtered[[col]], na.rm = TRUE)
      num_na <- sum(is.na(data_filtered[[col]]))
      col_quantile_val <- quantile(data_filtered[[col]], probs = q, na.rm = TRUE)
      if (num_na > 0) {
        data_filtered[[col]][is.na(data_filtered[[col]])] <- rnorm(num_na, mean = col_quantile_val, sd = col_sd)
      }
    }
  } else if (sample_wise==F){
    for (col in annotated_columns) {
      num_na <- sum(is.na(data_filtered[[col]]))
      if (num_na > 0) {
        data_filtered[[col]][is.na(data_filtered[[col]])] <- rnorm(num_na, mean = quantile_val, sd = sd_val*adj_std)
      }
    }
  }
  
  imputed_data <- data_filtered[missing_value_mask]
  non_imputed_data <- data_filtered[!missing_value_mask]
  imputed_combined <- unlist(imputed_data)
  non_imputed_combined <- unlist(non_imputed_data)
  
  plot_data_imputed <- data.frame(
    value = c(non_imputed_combined, imputed_combined),
    category = rep(c("Non-Imputed Values", "Imputed Values"), 
                   c(length(non_imputed_combined), length(imputed_combined)))
  )
  
  hist_plot3 <- ggplot(plot_data_imputed, aes(x = value, fill = category)) +
    geom_histogram(aes(y = ..density..), position = "identity", binwidth = 0.25, alpha = 0.5) +
    labs(title = "Distribution of data after imputation",
         x = "log2 Intensity",
         y = "Density") +
    theme_minimal()
  
  data_imputed <- data
  data_imputed[, annotated_columns] <- data_filtered
  
  if (ret==0){
    return(data_imputed)
  } else if (ret==1){
    print(hist_plot1)
  }else if (ret==2){
    print(hist_plot2)
  }else if (ret==3){
    print(hist_plot3)
  }
}


read_data <- function(file) {
  ext <- tools::file_ext(file)
  
  df <- switch(ext,
               csv = read_csv(file, show_col_types = FALSE),
               tsv = read_tsv(file, show_col_types = FALSE),
               txt = read_delim(file, delim = "\t", show_col_types = FALSE),
               xlsx = read_excel(file),
               stop("Invalid file type")
  )
  
  col_spec <- if ("Reverse" %in% names(df)) {
    cols(Reverse = col_character())
  } else {
    cols()
  }
  
  df <- switch(ext,
               csv = read_csv(file, col_types = col_spec),
               tsv = read_tsv(file, col_types = col_spec),
               txt = read_delim(file, delim = "\t", col_types = col_spec),
               xlsx = read_excel(file),
               stop("Invalid file type")
  )
  
  df <- rename_cols(df)
  
  rownames(df) <- NULL
  
  return(df)
}

simple_phos_site_plot <- function(data, filter=0) {
  data = data[data$PTM_localization >= filter,]
  phosprot = length(unique(data$Protein_group))
  data$Seq = data$UPD_seq
  for (i in 1:length(data$Seq)) {
    data$Seq[i] <- gsub("[[:punct:]]", "", data$Seq[i])
    data$Seq[i] = toupper(data$Seq[i])
  }
  phospep = length(unique(data$Seq))
  phossite = length(unique(data$PTM_Collapse_key))
  
  df = data.frame(term = c("Phosphoproteins", "Phosphopeptides", "Phosphosites"), count = c(phosprot, phospep, phossite))
  
  ggplot(df, aes(x = reorder(term, -count), y = count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_label(aes(label = count), fill = "white", color = "black", size = 5, fontface = "bold") +
    coord_flip() +
    theme_bw() +
    easy_remove_x_axis() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

prepare_KSEA <- function(data, meta, condition1, condition2, pvalue = FALSE) {
  convert_to_centered <- function(seq) {
    star_index <- regexpr("\\*", seq)[1]
    if (star_index <= 1) return(NA)
    
    mod_residue <- toupper(substr(seq, star_index - 1, star_index - 1))
    clean_seq <- paste0(
      substr(seq, 1, star_index - 2),
      mod_residue,
      substr(seq, star_index + 1, nchar(seq))
    )
    
    mod_pos <- star_index - 1
    left <- substr(clean_seq, 1, mod_pos - 1)
    right <- substr(clean_seq, mod_pos + 1, nchar(clean_seq))
    
    left_tail <- if (nchar(left) >= 6) {
      substr(left, nchar(left) - 5, nchar(left))
    } else {
      paste0(paste0(rep("_", 6 - nchar(left)), collapse = ""), left)
    }
    
    right_tail <- if (nchar(right) >= 6) {
      substr(right, 1, 6)
    } else {
      paste0(right, paste0(rep("_", 6 - nchar(right)), collapse = ""))
    }
    
    paste0("_", left_tail, mod_residue, right_tail, "_")
  }
  
  volcano_result <- volcano_data_f(
    data = data,
    meta = meta,
    condition1 = condition1,
    condition2 = condition2,
    workflow = "Phosphosite",
    paired = "Unpaired"
  )
  
  sequence_info <- data %>% 
    dplyr::select(PTM_Collapse_key, UPD_seq) %>%
    dplyr::distinct()
  
  volcano_merged <- volcano_result %>%
    dplyr::rename(PTM_Collapse_key = Phossite) %>%
    dplyr::inner_join(sequence_info, by = "PTM_Collapse_key") %>%
    dplyr::filter(!is.na(log2FC))
  
  volcano_merged$Sequence <- sapply(volcano_merged$UPD_seq, convert_to_centered)
  volcano_merged <- volcano_merged %>% dplyr::filter(!is.na(Sequence))
  
  if (pvalue) {
    ksea_input <- volcano_merged %>%
      dplyr::select(Sequence, logFC = log2FC, pval = adj_pval)
  } else {
    ksea_input <- volcano_merged %>%
      dplyr::select(Sequence, logFC = log2FC)
  }
  
  return(ksea_input)
}

prepare_tree_data_es <- function(data){
  data <- data %>%
    select(kinase, dominant_enrichment_value_log2)
  return(data)
}

prepare_tree_data_pv <- function(data) {
  data <- data %>%
    select(kinase, dominant_adjusted_p_value_log10_abs) %>%
    mutate(dominant_adjusted_p_value_log10_abs = round(dominant_adjusted_p_value_log10_abs, 3))
  return(data)
}

phossite_coverage_plot <- function(data, meta, id=FALSE){
  conditions <- unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id == TRUE){
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n(", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  
  data_classI <- data[data$PTM_localization >= 0.75, annotated_columns, drop = FALSE]
  data_notclassI <- data[data$PTM_localization < 0.75, annotated_columns, drop = FALSE]
  
  count_classI <- colSums(!is.na(data_classI))
  count_notclassI <- colSums(!is.na(data_notclassI))
  
  plot_data <- data.frame(
    sample = rep(annotated_columns, 2),
    PTM_localization = rep(c("Class I", "Not Class I"), each = length(annotated_columns)),
    count = c(count_classI, count_notclassI)
  )
  
  plot_data$sample <- factor(plot_data$sample, levels = meta$new_sample)
  
  ggplot(plot_data, aes(x = sample, y = count, fill = PTM_localization)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(
      title = "Phosphosites per Sample",
      x = "Sample",
      y = "Number of Phosphosites",
      fill = "Phosphosite"
    )
}

heatmap_plot_nmv <- function(data, meta, id=TRUE) {
  conditions = unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id)
  if (id==TRUE){
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  }else{
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered = data[, annotated_columns, drop = FALSE]
  data_filtered <- na.omit(data_filtered)
  annotation_col <- data.frame(Condition = meta$condition)
  rownames(annotation_col) <- rename_vector
  data_z <- scale(data_filtered, center = TRUE, scale = T)
  
  p=pheatmap(data_z,
             annotation_col = annotation_col, 
             annotation_names_col = FALSE,
             show_rownames = FALSE)
  return(p)
}

volcano_data_f <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1, workflow="Protein", paired = "Unpaired") {
  data[data == 0] <- NA
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  if (workflow=="Protein"){
    data_filtered <- data[, c("ProteinNames", annotated_columns1, annotated_columns2), drop = FALSE]} else if (workflow=="Phosphosite"){
      data_filtered <- data[, c("PTM_Collapse_key", annotated_columns1, annotated_columns2), drop = FALSE]}
  
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, -1])) < (length(annotated_columns1) + length(annotated_columns2) - 2), ]
  
  if (paired=="Unpaired"){
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
      mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
      mean2-mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- na.omit(as.numeric(row[annotated_columns1]))
      values2 <- na.omit(as.numeric(row[annotated_columns2]))
      if (length(values1) < 2 || length(values2) < 2) {
        return(NA)
      }
      t.test(values1, values2, var.equal = TRUE)$p.value
    })
  } else if (paired=="Paired"){
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 1) {
        return(NA)
      }
      mean1 <- mean(paired_values$values1, na.rm = TRUE)
      mean2 <- mean(paired_values$values2, na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 2) {
        return(NA)
      }
      t.test(paired_values$values1, paired_values$values2, paired = TRUE)$p.value
    })
  }
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  data_filtered <- data_filtered[valid_idx, ]
  
  adj_pvals <- p.adjust(pvals, method = "BH")
  
  if (workflow=="Protein"){
    volcano_data <- data.frame(
      Protein = data_filtered$`ProteinNames`,
      log2FC = log2fc,
      pval = pvals,
      adj_pval = adj_pvals
    )} else if (workflow=="Phosphosite"){
      volcano_data <- data.frame(
        Phossite = data_filtered$`PTM_Collapse_key`,
        log2FC = log2fc,
        pval = pvals,
        adj_pval = adj_pvals
      )}
  
  sig_pval <- in_pval
  sig_log2fc <- in_log2fc
  
  volcano_data <- volcano_data %>%
    mutate(
      significance = case_when(
        adj_pval < sig_pval & log2FC > sig_log2fc ~ "Upregulated",
        adj_pval < sig_pval & log2FC < -sig_log2fc ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_adj_pval = -log10(adj_pval)
    )
  return(volcano_data)
}

kinase_volcano <- function(data, in_pval = 0.1){
  log10_pval_thresh <- -log10(in_pval)
  
  volcano_data <- data %>%
    mutate(
      significance = case_when(
        dominant_enrichment_value_log2 > 0 & dominant_adjusted_p_value_log10_abs > log10_pval_thresh ~ "Upregulated",
        dominant_enrichment_value_log2 < 0 & dominant_adjusted_p_value_log10_abs > log10_pval_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  color_mapping <- c("Downregulated" = "blue", "Not significant" = "gray", "Upregulated" = "red")
  
  plotly_plot <- plot_ly(volcano_data, 
                         x = ~dominant_enrichment_value_log2, 
                         y = ~dominant_adjusted_p_value_log10_abs, 
                         text = ~kinase,  
                         type = 'scatter', 
                         mode = 'markers',
                         color = ~significance,
                         colors = color_mapping,
                         marker = list(size = 5)) %>%
    layout(
      title = "Kinases",
      xaxis = list(
        title = "log2 (enrichment value)",
        showline = F, 
        showgrid = T, 
        zeroline = FALSE
      ),
      yaxis = list(
        title = '-log10 Adjusted p-value', 
        showline = F, 
        showgrid = T, 
        zeroline = FALSE
      ),
      hovermode = 'closest',
      shapes = list(
        hline(log10_pval_thresh)
      )
    )
  
  return(plotly_plot)
}

vis_coverage <- function(data, protein, chunk_size=30, db){
  #Data
  data = subset(data, ProteinNames == protein)
  found_seq = data$Stripped.Sequence
  found_seq = unique(found_seq)
  
  sequence <- db$V1[db$name == protein]
  
  #Sequence
  n = nchar(sequence)
  
  result_seqs <- list(strrep("-", n))
  
  for (pep in found_seq) {
    match_positions <- gregexpr(pep, sequence)[[1]]
    for (start_pos in match_positions) {
      if (start_pos != -1) {
        end_pos <- start_pos + nchar(pep) - 1
        
        placed <- FALSE
        for (i in seq_along(result_seqs)) {
          if (substring(result_seqs[[i]], start_pos, end_pos) == strrep("-", nchar(pep))) {
            substring(result_seqs[[i]], start_pos, end_pos) <- pep
            placed <- TRUE
            break
          }
        }
        
        if (!placed) {
          new_result_seq <- strrep("-", n)
          substring(new_result_seq, start_pos, end_pos) <- pep
          result_seqs <- append(result_seqs, list(new_result_seq))
        }
      }
    }
  }
  
  num_seq = ""
  count = 10
  for (i in 1:n){
    if (i %% count == 0) {
      for (i in 1:(10-nchar(as.character(count)))){
        num_seq = paste0(num_seq, " ")
      }
      num_seq = paste0(num_seq, count)
      count=count+10
    }
  }
  
  starts <- seq(1, n, by = 10)
  chunks <- substring(sequence, starts, pmin(starts + 10 - 1, n))
  num_chunks <- substring(num_seq, starts, pmin(starts + 10 - 1, n))
  result_chunks_list <- lapply(result_seqs, function(res_seq) {
    result_chunks <- substring(res_seq, starts, pmin(starts + 10 - 1, n))
    paste(result_chunks, collapse = " ")
  })
  
  sequence <- paste(chunks, collapse = " ")
  num_seq <- paste(num_chunks, collapse = " ")
  
  chunk_size = chunk_size + chunk_size*0.1
  n = nchar(sequence)
  starts <- seq(1, n, by = chunk_size)
  
  lines <- substring(sequence, starts, pmin(starts + chunk_size - 1, n))
  num_lines <- substring(num_seq, starts, pmin(starts + chunk_size - 1, n))
  
  result_lines_list <- lapply(result_chunks_list, function(res_chunk) {
    substring(res_chunk, starts, pmin(starts + chunk_size - 1, n))
  })
  
  for (i in seq_along(lines)) {
    cat(num_lines[i], "\n")
    for (result_lines in rev(result_lines_list)) {
      cat(result_lines[i], "\n")
    }
    cat(lines[i], "\n")
  }
}


model_3d <- function(data, protein, db){
  data <- data[data$ProteinNames == protein, ]
  uniprot_ids <- unique(data$PG.ProteinAccessions)
  
  uniprot_information <- fetch_uniprot(
    uniprot_ids = uniprot_ids,
    columns = c("sequence", "xref_pdb")
  )
  
  sequence <- db$V1[db$name == protein]
  
  pdb_id <- strsplit(uniprot_information$xref_pdb, ";")[[1]][1]
  
  view <- r3dmol() %>%
    m_add_model(data = m_fetch_pdb(pdb_id), format = "pdb") %>%  
    m_set_style(style = m_style_cartoon(color = "green")) %>%  
    m_add_surface(style = m_style_surface(color = "green", opacity = 0.4)) %>% 
    m_zoom_to() 
  
  found_seq <- unique(data$Stripped.Sequence)
  for (pep in found_seq) {
    positions <- gregexpr(pep, sequence)[[1]]
    
    for (start_pos in positions) {
      if (start_pos != -1) {
        end_pos <- start_pos + nchar(pep) - 1
        view <- view %>%
          m_set_style(
            sel = m_sel(resi = start_pos:end_pos),  
            style = m_style_cartoon(color = "red")   
          )
      }
    }
  }
  view
}

match_volc_and_org_data <- function(volc_data, org_data) {
  volc_data <- volc_data %>%
    rename(PTM_Collapse_key = Phossite)
  
  data_reduced <- org_data %>%
    select(Protein_group, Gene_group, UPD_seq, PTM_Collapse_key)
  
  data_reduced$UPD_seq <- toupper(gsub("[^A-Za-z0-9]", "", data_reduced$UPD_seq))
  
  matched_data <- volc_data %>%
    inner_join(data_reduced, by = "PTM_Collapse_key")
  
  matched_data <- matched_data %>%
    mutate(FC = 2^log2FC)
  
  matched_data <- matched_data %>%
    rename(
      Protein = Protein_group,
      Gene = Gene_group,
      Peptide = UPD_seq,
      p = pval
    ) %>%
    select(PTM_Collapse_key, Protein, Gene, Peptide, FC, p)
  
  matched_data <- matched_data %>%
    mutate(Residue.Both = str_extract(PTM_Collapse_key, "(?<=_)[^_]+(?=_[^_]+$)")) %>%
    select(-PTM_Collapse_key)
  
  matched_data <- matched_data %>%
    select(Protein, Gene, Peptide, Residue.Both, p, FC)
  
  return(matched_data)
}

kinact_kinase_activity <- function(matched_data, top_n = 0, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5) {
  KSEA_scores = KSEA.Scores(KSData, matched_data, NetworKIN=NetworKIN, NetworKIN.cutoff=NetworKIN.cutoff)
  
  KSEA_scores <- KSEA_scores %>%
    filter(m >= m.cutoff)
  
  if (top_n != 0) {
    if (top_n > 0) {
      KSEA_scores <- KSEA_scores %>% 
        arrange(desc(z.score)) %>% 
        slice(1:top_n)
    } else if (top_n < 0) {
      KSEA_scores <- KSEA_scores %>% 
        arrange(z.score) %>% 
        slice(1:abs(top_n))
    }
  }
  
  KSEA_scores <- KSEA_scores %>%
    mutate(Color = case_when(
      z.score > 0 & p.value < 0.05 ~ "red",    
      z.score < 0 & p.value < 0.05 ~ "blue", 
      TRUE ~ "grey"                         
    ))
  
  ggplot_kinase <- ggplot(KSEA_scores, aes(x = reorder(Kinase.Gene, z.score), y = z.score, fill = Color)) +
    geom_bar(stat = "identity") +  
    theme_minimal() + 
    labs(
      title = "Bar Plot of Kinase Scores",
      x = "Kinase Gene",
      y = "Z Score"
    ) +
    coord_flip() +  
    scale_fill_identity() +  
    theme(
      axis.text.x = element_text(size = 8), 
      axis.text.y = element_text(size = 8), 
      plot.title = element_text(size = 10),
      legend.position = "none" 
    )
  
  plotly_plot <- ggplotly(ggplot_kinase)
  
  return(plotly_plot)
}

downstream_phossite_volc <- function(matched_data, Kinase, NetworKIN = FALSE){
  new_data = KSEA.KS_table(KSData, matched_data, NetworKIN = NetworKIN)
  
  new_data= new_data %>%
    filter(Kinase.Gene == Kinase)
  
  merged_data <- inner_join(matched_data, new_data, 
                            by = c("Gene" = "Substrate.Gene", 
                                   "Residue.Both" = "Substrate.Mod"))
  
  merged_data <- merged_data %>%
    mutate(neg_log10_p = -log10(p),
           Protein_Residue_Both = paste(Protein, Residue.Both, sep = "_"),
           Color = case_when(
             log2FC > 1 & p < 0.05 ~ "red",  
             log2FC < -1 & p < 0.05 ~ "blue", 
             TRUE ~ "grey"                        
           )
    )
  
  ggplot_volcano <- ggplot(merged_data, aes(x = log2FC, y = neg_log10_p, text = Protein_Residue_Both)) +
    geom_point(aes(color = Color), alpha = 0.6) + 
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgreen") +  
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgreen") +  
    theme_minimal() +
    labs(
      title = Kinase,
      x = "log2 fold change",
      y = "-log10 p-value"
    ) +
    scale_color_identity() +
    theme(
      plot.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
  
  plotly_volcano <- ggplotly(ggplot_volcano, tooltip = "text")
  
  return(plotly_volcano)
}

init_colors <- function(color) {
  if (color == "Default") {
    if (exists("plot_colors", envir = .GlobalEnv)) {
      rm("plot_colors", envir = .GlobalEnv)
    }
  } else if (color == "Default16") {
    plot_colors <<- c(
      "#F8766D", "#E68613", "#CD9600", "#ABA300", "#7CAE00", "#0BB702",
      "#00BE67", "#00C19A", "#00BFC4", "#00B8E7", "#00A9FF", "#8494FF",
      "#C77CFF", "#ED68ED", "#FF61CC", "#FF68A1", "black"
    )
  } else if (color == "Warm/Cold") {
    plot_colors <<- c(
      "blue", "red", "green", "orange", "purple", "darkgoldenrod", "cyan",
      "gold", "orchid", "salmon", "navy", "brown", "steelblue", "coral",
      "skyblue", "peachpuff", "orchid"
    )
  } else if (color == "Black/Grey") {
    plot_colors <<- c("black", "grey")
  } else if (color == "Yue7") {
    plot_colors <<- c(
      "#2D5F85", "#5184B2", "#AAD4F8", "#F2F5FA",
      "#F1A7B5", "#D55276", "#AB3A54"
    )
  } else if (color == "Mario Document Input") {
    file_path <- file.path("www", "colors", "colors.txt")
    script_dir <- dirname(sys.frame(1)$ofile %||%
                            rstudioapi::getActiveDocumentContext()$path)
    full_path <- file.path(script_dir, file_path)
    plot_colors <<- readLines(full_path, warn = FALSE)
  }
}

coverage_plot_summary <- function(data, meta, id=TRUE, color_package=TRUE) {
  data[data == 0] <- NA
  conditions <- unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id == TRUE) {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n (", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  
  data_filtered[!is.na(data_filtered)] <- 1
  data_melted <- melt(data_filtered, variable.name = "Sample", value.name = "Value")
  data_annotated <- merge(data_melted, meta, by.x = "Sample", by.y = "new_sample")
  data_annotated$Sample <- factor(data_annotated$Sample, levels = meta$new_sample)
  data_annotated$condition <- factor(data_annotated$condition, levels = conditions)
  
  data_annotated_unique <- data_annotated %>%
    group_by(Sample) %>%
    summarize(Value = sum(Value, na.rm = TRUE),
              condition = first(condition))
  
  summary_data <- data_annotated_unique %>%
    group_by(condition) %>%
    summarize(mean_value = mean(Value, na.rm = TRUE),
              sd_value = sd(Value, na.rm = TRUE),
              .groups = "drop")
  
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(conditions)) {
      plot_colors <- rep(plot_colors, length.out = length(conditions))
    }
  }
  
  if ("ProteinNames" %in% colnames(data)) {
    p <- ggplot() +
      geom_bar(data = summary_data, aes(x = condition, y = mean_value, fill = condition),
               stat = "identity", position = "dodge", color = "black") +
      geom_errorbar(data = summary_data,
                    aes(x = condition, ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                    width = 0.2, position = position_dodge(0.9)) +
      geom_point(data = data_annotated_unique, aes(x = condition, y = Value),
                 position = position_jitter(width = 0.2), size = 2, alpha = 0.7) +
      theme_minimal() +
      labs(title = "Proteins per sample",
           x = "Condition", y = "Number of proteins", fill = "Condition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_hline(yintercept = length(data$ProteinNames), linetype = "dashed", color = "red")
  } else if ("PTM_Collapse_key" %in% colnames(data)) {
    p <- ggplot() +
      geom_bar(data = summary_data, aes(x = condition, y = mean_value, fill = condition),
               stat = "identity", position = "dodge", color = "black") +
      geom_errorbar(data = summary_data,
                    aes(x = condition, ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                    width = 0.2, position = position_dodge(0.9)) +
      geom_point(data = data_annotated_unique, aes(x = condition, y = Value),
                 position = position_jitter(width = 0.2), size = 2, alpha = 0.7) +
      theme_minimal() +
      labs(title = "Phosphosites per sample",
           x = "Condition", y = "Number of phosphosites", fill = "Condition") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_hline(yintercept = length(data$PTM_Collapse_key), linetype = "dashed", color = "red")
  }
  
  if (exists("plot_colors")) {
    p <- p + scale_fill_manual(values = plot_colors)
  }
  print(p)
}


phossite_coverage_plot_summary <- function(data, meta, id=FALSE){
  conditions <- unique(meta$condition)
  meta$sample <- as.character(meta$sample)
  meta$id <- sapply(meta$sample, extract_id_or_number)
  
  if (id == TRUE){
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number(), "\n(", id, ")")) %>%
      ungroup()
  } else {
    meta <- meta %>%
      group_by(condition) %>%
      mutate(new_sample = paste0(condition, "_", row_number())) %>%
      ungroup()
  }
  
  rename_vector <- setNames(meta$new_sample, meta$sample)
  colnames(data) <- ifelse(colnames(data) %in% meta$sample, rename_vector[colnames(data)], colnames(data))
  annotated_columns <- meta$new_sample
  
  data_classI <- data[data$PTM_localization >= 0.75, annotated_columns, drop = FALSE]
  data_notclassI <- data[data$PTM_localization < 0.75, annotated_columns, drop = FALSE]
  count_classI <- colSums(!is.na(data_classI))
  count_notclassI <- colSums(!is.na(data_notclassI))
  
  plot_data <- data.frame(
    sample = rep(annotated_columns, 2),
    PTM_localization = rep(c("Class I", "Not Class I"), each = length(annotated_columns)),
    count = c(count_classI, count_notclassI)
  )
  
  summarized_plot_data <- plot_data %>%
    mutate(condition = sub("_.*", "", sample)) %>%
    group_by(condition, PTM_localization) %>%    
    summarize(
      mean_count = mean(count, na.rm = TRUE),      
      sd_count = sd(count, na.rm = TRUE),          
      .groups = "drop"                           
    )
  
  plot_data <- plot_data %>%
    mutate(condition = sub("_.*", "", sample))
  
  ggplot() +
    geom_bar(
      data = summarized_plot_data,
      aes(x = condition, y = mean_count, fill = PTM_localization),
      stat = "identity", position = position_dodge(width = 0.8)
    ) +
    geom_errorbar(
      data = summarized_plot_data,
      aes(x = condition, ymin = mean_count - sd_count, ymax = mean_count + sd_count, group = PTM_localization),
      position = position_dodge(width = 0.8), width = 0.2
    ) +
    geom_point(
      data = plot_data,
      aes(x = condition, y = count, fill = PTM_localization),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = 2, color = "black", alpha = 0.7
    ) +
    theme_minimal() +
    labs(
      title = "Phosphosites per condition",
      x = "Condition",
      y = "Number of phosphosites",
      fill = "Phosphosite"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right" 
    )
}

volcano_plot <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1,
                         workflow = "Protein", paired = "Unpaired", uncorrected = FALSE) {
  
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  
  if (workflow == "Protein") {
    data_filtered <- data[, c("ProteinNames", annotated_columns1, annotated_columns2), drop = FALSE]
  } else if (workflow == "Phosphosite") {
    data_filtered <- data[, c("PTM_Collapse_key", annotated_columns1, annotated_columns2), drop = FALSE]
  }
  
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, -1])) < (length(annotated_columns1) + length(annotated_columns2) - 2), ]
  
  if (paired == "Unpaired") {
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
      mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- na.omit(as.numeric(row[annotated_columns1]))
      values2 <- na.omit(as.numeric(row[annotated_columns2]))
      if (length(values1) < 2 || length(values2) < 2) return(NA)
      t.test(values1, values2, var.equal = TRUE)$p.value
    })
  } else if (paired == "Paired") {
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 1) return(NA)
      mean2 <- mean(paired_values$values2, na.rm = TRUE)
      mean1 <- mean(paired_values$values1, na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 2) return(NA)
      t.test(paired_values$values1, paired_values$values2, paired = TRUE)$p.value
    })
  }
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  data_filtered <- data_filtered[valid_idx, ]
  
  adj_pvals <- if (uncorrected) pvals else p.adjust(pvals, method = "BH")
  
  label_column <- if (workflow == "Protein") "ProteinNames" else "PTM_Collapse_key"
  volcano_data <- data.frame(
    Feature = data_filtered[[label_column]],
    log2FC = log2fc,
    pval = pvals,
    adj_pval = adj_pvals
  )
  
  sig_pval <- in_pval
  sig_log2fc <- in_log2fc
  
  volcano_data <- volcano_data %>%
    mutate(
      significance = case_when(
        adj_pval < sig_pval & log2FC > sig_log2fc ~ "Upregulated",
        adj_pval < sig_pval & log2FC < -sig_log2fc ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_pval = -log10(adj_pval)
    )
  
  color_mapping <- c("Downregulated" = "blue", "Not significant" = "gray", "Upregulated" = "red")
  y_axis_label <- if (uncorrected) "-log10 p-value" else "-log10 adj. p-value"
  
  plot_ly(volcano_data,
          x = ~log2FC,
          y = ~neg_log10_pval,
          text = ~Feature,
          type = 'scatter',
          mode = 'markers',
          color = ~significance,
          colors = color_mapping,
          marker = list(size = 5)) %>%
    layout(
      title = paste0(condition1, " vs. ", condition2),
      xaxis = list(
        title = paste0("log2 fold change (", condition2, " - ", condition1, ")"),
        showline = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = y_axis_label,
        showline = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      hovermode = 'closest',
      shapes = list(
        vline(sig_log2fc),
        vline(-sig_log2fc),
        hline(-log10(sig_pval))
      )
    )
}

volcano_data_f <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1, workflow="Protein", paired = "Unpaired", uncorrected = FALSE) {
  data[data == 0] <- NA
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  
  if (workflow == "Protein") {
    data_filtered <- data[, c("ProteinNames", annotated_columns1, annotated_columns2), drop = FALSE]
  } else if (workflow == "Phosphosite") {
    data_filtered <- data[, c("PTM_Collapse_key", annotated_columns1, annotated_columns2), drop = FALSE]
  }
  
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, -1])) < (length(annotated_columns1) + length(annotated_columns2) - 2), ]
  
  if (paired == "Unpaired") {
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
      mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- na.omit(as.numeric(row[annotated_columns1]))
      values2 <- na.omit(as.numeric(row[annotated_columns2]))
      if (length(values1) < 2 || length(values2) < 2) {
        return(NA)
      }
      t.test(values1, values2, var.equal = TRUE)$p.value
    })
  } else if (paired == "Paired") {
    log2fc <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 1) {
        return(NA)
      }
      mean1 <- mean(paired_values$values1, na.rm = TRUE)
      mean2 <- mean(paired_values$values2, na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, -1], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 2) {
        return(NA)
      }
      t.test(paired_values$values1, paired_values$values2, paired = TRUE)$p.value
    })
  }
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  data_filtered <- data_filtered[valid_idx, ]
  
  if (uncorrected) {
    adj_pvals <- pvals
  } else {
    adj_pvals <- p.adjust(pvals, method = "BH")
  }
  
  if (workflow == "Protein") {
    volcano_data <- data.frame(
      Protein = data_filtered$`ProteinNames`,
      log2FC = log2fc,
      pval = pvals,
      adj_pval = adj_pvals
    )
  } else if (workflow == "Phosphosite") {
    volcano_data <- data.frame(
      Phossite = data_filtered$`PTM_Collapse_key`,
      log2FC = log2fc,
      pval = pvals,
      adj_pval = adj_pvals
    )
  }
  
  sig_pval <- in_pval
  sig_log2fc <- in_log2fc
  
  volcano_data <- volcano_data %>%
    mutate(
      significance = case_when(
        adj_pval < sig_pval & log2FC > sig_log2fc ~ "Upregulated",
        adj_pval < sig_pval & log2FC < -sig_log2fc ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_adj_pval = -log10(adj_pvals)
    )
  
  return(volcano_data)
}


tsne_plot <- function(data, meta, color_package = TRUE, perplexity = 30, max_iter = 1000, header = TRUE, legend = TRUE) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("The 'Rtsne' package is required for t-SNE. Please install it with install.packages('Rtsne').")
  }
  
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  meta$condition <- factor(meta$condition, levels = unique(meta$condition))
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  data_filtered <- data_filtered[complete.cases(data_filtered), ]
  
  transposed_expr <- t(data_filtered)
  
  tryCatch({
    zero_variance_columns <- which(apply(transposed_expr, 2, var) == 0)
    if (length(zero_variance_columns) > 0) {
      transposed_expr <- transposed_expr[, -zero_variance_columns, drop = FALSE]
    }
  }, error = function(e) {
    message("An error occurred while processing zero variance columns: ", e$message)
  })
  
  tryCatch({
    if (any(duplicated(transposed_expr))) {
      transposed_expr <- transposed_expr + matrix(rnorm(length(transposed_expr), mean = 0, sd = 1e-4), nrow = nrow(transposed_expr))
    }
    
    tsne_result <- Rtsne::Rtsne(
      transposed_expr,
      perplexity = perplexity,
      max_iter = max_iter,
      verbose = FALSE,
      check_duplicates = FALSE
    )
    
    tsne_scores <- as.data.frame(tsne_result$Y)
    colnames(tsne_scores) <- c("tSNE1", "tSNE2")
    tsne_scores$sample <- rownames(transposed_expr)
    
    tsne_scores <- merge(tsne_scores, meta, by.x = "sample", by.y = "sample")
    tsne_scores$condition <- factor(tsne_scores$condition, levels = levels(meta$condition))
    
    tsne <- ggplot(tsne_scores, aes(x = tSNE1, y = tSNE2, color = condition)) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal() +
      labs(
        x = "t-SNE 1",
        y = "t-SNE 2"
      )
    
    if (header) {
      tsne <- tsne + labs(title = "t-SNE Plot")
    }
    
    tsne <- tsne + theme(legend.position = if (legend) "right" else "none")
    
    if (exists("plot_colors")) {
      tsne <- tsne + scale_color_manual(values = plot_colors)
    }
    
    return(tsne)
  }, error = function(e) {
    message("An error occurred during t-SNE computation: ", e$message)
    return(NULL)
  })
}


umap_plot <- function(data, meta, color_package = TRUE, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", header = TRUE, legend = TRUE) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("The 'uwot' package is required for UMAP. Please install it with install.packages('uwot').")
  }
  
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  meta$condition <- factor(meta$condition, levels = unique(meta$condition))
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  data_filtered <- data_filtered[complete.cases(data_filtered), ]
  
  transposed_expr <- t(data_filtered)
  
  tryCatch({
    zero_variance_columns <- which(apply(transposed_expr, 2, var) == 0)
    if (length(zero_variance_columns) > 0) {
      transposed_expr <- transposed_expr[, -zero_variance_columns, drop = FALSE]
    }
  }, error = function(e) {
    message("An error occurred while processing zero variance columns: ", e$message)
  })
  
  tryCatch({
    umap_result <- uwot::umap(
      transposed_expr,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      metric = metric,
      scale = TRUE
    )
    
    umap_scores <- as.data.frame(umap_result)
    colnames(umap_scores) <- c("UMAP1", "UMAP2")
    umap_scores$sample <- rownames(transposed_expr)
    
    umap_scores <- merge(umap_scores, meta, by.x = "sample", by.y = "sample")
    umap_scores$condition <- factor(umap_scores$condition, levels = levels(meta$condition))
    
    umap <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, color = condition)) +
      geom_point(size = 3, alpha = 0.7) +  
      theme_minimal() +
      labs(
        x = "UMAP 1",
        y = "UMAP 2"
      )
    
    if (header) {
      umap <- umap + labs(title = "UMAP Plot")
    }
    
    umap <- umap + theme(legend.position = if (legend) "right" else "none")
    
    if (exists("plot_colors")) {
      umap <- umap + scale_color_manual(values = plot_colors)
    }
    
    return(umap)
  }, error = function(e) {
    message("An error occurred during UMAP computation: ", e$message)
    return(NULL)
  })
}


dim_func <- function(data, meta, method, color_package = TRUE, header = TRUE, legend = TRUE, dot_size=3) {
  plot_obj <- switch(method,
                     "PCA"  = pca_plot(data, meta, color_package, header = header, legend = legend, dot_size = dot_size),
                     "tSNE" = tsne_plot(data, meta, color_package, header = header, legend = legend),
                     "UMAP" = umap_plot(data, meta, color_package, header = header, legend = legend)
  )
  print(plot_obj)
}


customSpinnerWrapper <- function(plotExpr, workmode = "default") {
  if (workmode == "funny") {
    withSpinner(
      plotExpr,
      image = "loading/cat_loading.gif"
    )
  } else {
    withSpinner(
      plotExpr,
      type = 1,
      color = "#337ab7",
      size = 1
    )
  }
}

customSpinnerWrapper <- function(plotExpr, workmode = "default") {
  if (workmode == "funny") {
    gif_files <- list.files("www/loading", pattern = "\\.gif$", full.names = FALSE)
    
    if (length(gif_files) == 0) {
      warning("No GIFs found in www/loading/. Using default spinner.")
      return(
        withSpinner(
          plotExpr,
          type = 1,
          color = "#337ab7",
          size = 1
        )
      )
    }
    
    selected_gif <- sample(gif_files, 1)
    
    withSpinner(
      plotExpr,
      image = file.path("loading", selected_gif)
    )
  } else {
    withSpinner(
      plotExpr,
      type = 1,
      color = "#337ab7",
      size = 1
    )
  }
}

transform_meta <- function(data, meta, sort=TRUE){
  meta <- meta[, c("File Name", "Condition")]
  colnames(meta) <- c("sample", "condition")
  colnames_list <- colnames(data)
  meta$sample <- sapply(meta$sample, function(sample_str) {
    match <- colnames_list[grepl(sample_str, colnames_list)]
    if (length(match) > 0) match[1] else sample_str 
  })
  if(sort == TRUE){
    meta <- meta[order(meta$condition), ]
  }
  
  return(meta)
}

ctl <- function(x) {
  ifelse(x == "TRUE", TRUE,
         ifelse(x == "FALSE", FALSE, x))
}

read_sample_info <- function(file_path = "description.txt") {
  if (!file.exists(file_path)) return(NULL)
  
  file_text <- readChar(file_path, file.info(file_path)$size)
  
  if (substr(file_text, 1, 1) != "@") stop("File must start with an '@' symbol.")
  
  parts <- unlist(strsplit(file_text, "@", fixed = TRUE))
  parts <- trimws(parts[nzchar(trimws(parts))])
  
  if (length(parts) %% 2 != 0) stop("Number of headers and texts does not match.")
  
  headers <- parts[seq(1, length(parts), by = 2)]
  texts <- parts[seq(2, length(parts), by = 2)]
  
  df <- data.frame(header = headers, text = texts, stringsAsFactors = FALSE)
  return(df)
}

volcano_data_f <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1, workflow = "Protein", paired = "Unpaired", uncorrected = FALSE) {
  data[data == 0] <- NA
  
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  all_columns <- c(annotated_columns1, annotated_columns2)
  
  if (workflow == "Protein") {
    data_filtered <- data[, c("ProteinNames", all_columns), drop = FALSE]
    id_col <- "ProteinNames"
    id_col_name <- "Protein"
  } else if (workflow == "Phosphosite") {
    data_filtered <- data[, c("PTM_Collapse_key", all_columns), drop = FALSE]
    id_col <- "PTM_Collapse_key"
    id_col_name <- "Phossite"
  } else {
    stop("Invalid workflow specified.")
  }
  
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, all_columns])) < (length(all_columns) - 2), ]
  
  if (paired == "Unpaired") {
    log2fc <- apply(data_filtered[, all_columns], 1, function(row) {
      mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
      mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
      mean2 - mean1
    })
    
    pvals <- apply(data_filtered[, all_columns], 1, function(row) {
      values1 <- na.omit(as.numeric(row[annotated_columns1]))
      values2 <- na.omit(as.numeric(row[annotated_columns2]))
      if (length(values1) < 2 || length(values2) < 2) return(NA)
      t.test(values1, values2, var.equal = TRUE)$p.value
    })
    
  } else if (paired == "Paired") {
    log2fc <- apply(data_filtered[, all_columns], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 1) return(NA)
      mean(paired_values$values2) - mean(paired_values$values1)
    })
    
    pvals <- apply(data_filtered[, all_columns], 1, function(row) {
      values1 <- as.numeric(row[annotated_columns1])
      values2 <- as.numeric(row[annotated_columns2])
      paired_values <- na.omit(data.frame(values1, values2))
      if (nrow(paired_values) < 2) return(NA)
      t.test(paired_values$values2, paired_values$values1, paired = TRUE)$p.value
    })
  }
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  
  if (uncorrected) {
    adj_pvals <- pvals
  } else {
    adj_pvals <- p.adjust(pvals, method = "BH")
  }
  
  data_filtered <- data_filtered[valid_idx, ]
  id_column <- data_filtered[[id_col]]
  
  volcano_data <- data_filtered %>%
    mutate(
      !!id_col_name := id_column,
      log2FC = log2fc,
      pval = pvals,
      adj_pval = adj_pvals,
      neg_log10_adj_pval = -log10(adj_pval),
      significance = case_when(
        adj_pval < in_pval & log2FC > in_log2fc ~ "Upregulated",
        adj_pval < in_pval & log2FC < -in_log2fc ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  volcano_data <- volcano_data[, c(id_col_name, all_columns, "log2FC", "pval", "adj_pval", "neg_log10_adj_pval", "significance")]
  
  return(volcano_data)
}

first_digit_distribution <- function(data, meta) {
  sample_columns <- meta$sample
  data_values <- data[, sample_columns, drop = FALSE]
  
  data_values[data_values == 0] <- NA
  values_vector <- as.vector(unlist(data_values))
  values_vector <- values_vector[!is.na(values_vector) & values_vector > 0]
  
  first_digits <- as.numeric(substr(as.character(values_vector), 1, 1))
  first_digits <- first_digits[first_digits %in% 1:9]
  
  digit_freq <- as.data.frame(table(first_digits))
  digit_freq$Frequency <- digit_freq$Freq / sum(digit_freq$Freq)
  digit_freq$first_digits <- as.factor(digit_freq$first_digits)
  
  benford <- data.frame(
    first_digits = as.factor(1:9),
    Benford = log10(1 + 1 / (1:9))
  )
  
  plot_data <- merge(digit_freq, benford, by = "first_digits", all = TRUE)
  
  p <- ggplot(plot_data, aes(x = first_digits)) +
    geom_bar(aes(y = Frequency), stat = "identity", fill = "darkgreen") +
    geom_text(aes(y = Frequency, label = scales::percent(Frequency, accuracy = 0.1)),
              vjust = -0.5, size = 4) +
    geom_line(aes(y = Benford, group = 1), color = "red", linewidth = 1) +
    geom_point(aes(y = Benford), color = "red", size = 2) +
    labs(x = "First Digit", y = "Relative Frequency") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.1))) +
    theme_minimal()
  
  return(p)
}

data_pattern_structure <- function(data, meta) {
  sample_columns <- meta$sample
  data_values <- data[, sample_columns, drop = FALSE]
  all_values <- na.omit(unlist(data_values))
  value_freq <- table(all_values)
  freq_of_freq <- table(value_freq)
  
  freq_of_freq_percent <- 100 * freq_of_freq / sum(freq_of_freq)
  
  bar_midpoints <- barplot(freq_of_freq_percent,
                           xlab = "Duplicate Number Occurrences",
                           ylab = "Percentage (%)",
                           col = "skyblue",
                           ylim = c(0, max(freq_of_freq_percent) * 1.2))
  
  text(x = bar_midpoints,
       y = freq_of_freq_percent,
       labels = paste0(round(freq_of_freq_percent, 3), "%"),
       pos = 3, 
       cex = 0.8,
       col = "blue")
  
  invisible(freq_of_freq_percent)
}

read_export_count <- function(file_path) {
  if (file.exists(file_path)) {
    count <- suppressWarnings(as.integer(readLines(file_path, warn = FALSE)))
    if (is.na(count)) return(0)
    return(count)
  } else {
    return(0)
  }
}

increment_export_count <- function(file_path) {
  count <- read_export_count(file_path)
  count <- count + 1
  writeLines(as.character(count), file_path)
  return(count)
}

all_volcano_data <- function(data, meta, in_pval = 0.05, in_log2fc = 1, workflow = "Protein", paired = "Unpaired") {
  conditions <- unique(meta$condition)
  combs <- combn(conditions, 2, simplify = FALSE)
  
  wb <- createWorkbook()
  
  for (comp in combs) {
    cond1 <- comp[1]
    cond2 <- comp[2]
    
    volcano_df <- volcano_data_f(data, meta, cond1, cond2, in_pval, in_log2fc, workflow, paired)
    
    if (is.null(volcano_df) || nrow(volcano_df) == 0) next
    
    sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", paste0(cond1, "_vs_", cond2)), 1, 31)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, volcano_df)
  }
  
  tmpfile <- tempfile(fileext = ".xlsx")
  saveWorkbook(wb, tmpfile, overwrite = TRUE)
  return(tmpfile)
}

qqnorm_plot <- function(data, meta) {
  sample_columns <- meta$sample
  data_values <- data[, sample_columns, drop = FALSE]
  data_values[data_values == 0] <- NA
  
  values_vector <- as.vector(unlist(data_values))
  values_vector <- values_vector[!is.na(values_vector) & values_vector > 0]
  
  if (length(values_vector) > 100000) {
    set.seed(187) 
    values_vector <- sample(values_vector, 100000)
  }
  
  qqnorm(values_vector, main = "QQ Plot", pch = 19, col = "blue")
  qqline(values_vector, col = "red", lwd = 2)
}

pca_plot_interactive <- function(data, meta, color_package = TRUE, header = TRUE, legend = TRUE) {
  if (exists("plot_colors")) {
    if (length(plot_colors) < length(unique(meta$condition))) {
      plot_colors <- rep(plot_colors, length.out = length(unique(meta$condition)))
    }
  }
  
  meta$condition <- factor(meta$condition, levels = unique(meta$condition))
  annotated_columns <- meta$sample
  data_filtered <- data[, annotated_columns, drop = FALSE]
  data_filtered <- data_filtered[complete.cases(data_filtered), ]
  
  transposed_expr <- t(data_filtered)
  
  tryCatch({
    zero_variance_columns <- which(apply(transposed_expr, 2, var) == 0)
    if (length(zero_variance_columns) > 0) {
      transposed_expr <- transposed_expr[, -zero_variance_columns, drop = FALSE]
    }
  }, error = function(e) {
    message("An error occurred while processing zero variance columns: ", e$message)
  })
  
  tryCatch({
    pca_result_transposed <- prcomp(transposed_expr, scale. = TRUE)
    explained_variance <- (pca_result_transposed$sdev^2) / sum(pca_result_transposed$sdev^2) * 100
    pca_scores <- as.data.frame(pca_result_transposed$x)
    pca_scores$sample <- rownames(pca_scores)
    
    pca_scores <- merge(pca_scores, meta, by = "sample")
    pca_scores$condition <- factor(pca_scores$condition, levels = levels(meta$condition))
    
    x_label <- sprintf("Principal Component 1 - %.2f%% variance", explained_variance[1])
    y_label <- sprintf("Principal Component 2 - %.2f%% variance", explained_variance[2])
    
    plt <- plot_ly(
      data = pca_scores,
      x = ~PC1,
      y = ~PC2,
      type = 'scatter',
      mode = 'markers',
      color = ~condition,
      colors = if (exists("plot_colors")) plot_colors else NULL,
      text = ~paste("Sample:", sample, "<br>Condition:", condition),
      marker = list(size = 8, opacity = 0.8),
      hoverinfo = "text"
    )
    
    plt <- plt %>%
      layout(
        title = if (header) "PCA Plot" else NULL,
        zeroline = FALSE,   
        xaxis = list(title = x_label),
        yaxis = list(title = y_label),
        showlegend = legend
      )
    
    return(plt)
  }, error = function(e) {
    message("An error occurred during PCA computation: ", e$message)
    return(NULL)
  })
}

volcano_plot_sim <- function(data, meta, condition1, condition2, in_pval = 0.05, in_log2fc = 1,
                             workflow = "Protein", mod_var = 1, mod_n = 0) {
  
  annotated_columns1 <- meta$sample[meta$condition == condition1]
  annotated_columns2 <- meta$sample[meta$condition == condition2]
  
  if (workflow == "Protein") {
    data_filtered <- data[, c("ProteinNames", annotated_columns1, annotated_columns2), drop = FALSE]
  } else if (workflow == "Phosphosite") {
    data_filtered <- data[, c("PTM_Collapse_key", annotated_columns1, annotated_columns2), drop = FALSE]
  }
  
  data_filtered <- data_filtered[rowSums(is.na(data_filtered[, -1])) < (length(annotated_columns1) + length(annotated_columns2) - 2), ]
  
  log2fc <- apply(data_filtered[, -1], 1, function(row) {
    mean1 <- mean(as.numeric(row[annotated_columns1]), na.rm = TRUE)
    mean2 <- mean(as.numeric(row[annotated_columns2]), na.rm = TRUE)
    mean2 - mean1
  })
  
  pvals <- apply(data_filtered[, -1], 1, function(row) {
    x <- na.omit(as.numeric(row[annotated_columns1]))
    y <- na.omit(as.numeric(row[annotated_columns2]))
    
    n1_real <- length(x)
    n2_real <- length(y)
    if (n1_real < 2 || n2_real < 2) return(NA)
    
    mean1 <- mean(x)
    mean2 <- mean(y)
    var1 <- var(x)
    var2 <- var(y)
    
    n1 <- if (mod_n > 0) mod_n else n1_real
    n2 <- if (mod_n > 0) mod_n else n2_real
    
    pooled_var <- ((n1_real - 1) * var1 + (n2_real - 1) * var2) / (n1_real + n2_real - 2)
    pooled_var <- pooled_var * mod_var
    
    se <- sqrt(pooled_var * (1 / n1 + 1 / n2))
    t_stat <- (mean2 - mean1) / se
    df <- n1 + n2 - 2
    2 * pt(-abs(t_stat), df)
  })
  
  valid_idx <- !is.na(pvals)
  log2fc <- log2fc[valid_idx]
  pvals <- pvals[valid_idx]
  data_filtered <- data_filtered[valid_idx, ]
  
  adj_pvals <- p.adjust(pvals, method = "BH")
  
  label_column <- if (workflow == "Protein") "ProteinNames" else "PTM_Collapse_key"
  volcano_data <- data.frame(
    Feature = data_filtered[[label_column]],
    log2FC = log2fc,
    pval = pvals,
    adj_pval = adj_pvals
  )
  
  sig_pval <- in_pval
  sig_log2fc <- in_log2fc
  
  volcano_data <- volcano_data %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        adj_pval < sig_pval & log2FC > sig_log2fc ~ "Upregulated",
        adj_pval < sig_pval & log2FC < -sig_log2fc ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      neg_log10_pval = -log10(adj_pval)
    )
  
  color_mapping <- c("Downregulated" = "blue", "Not significant" = "gray", "Upregulated" = "red")
  
  plot_ly(volcano_data,
          x = ~log2FC,
          y = ~neg_log10_pval,
          text = ~Feature,
          type = 'scatter',
          mode = 'markers',
          color = ~significance,
          colors = color_mapping,
          marker = list(size = 5)) %>%
    layout(
      title = paste0(condition1, " vs. ", condition2),
      xaxis = list(
        title = paste0("log2 fold change (", condition2, " - ", condition1, ")"),
        showline = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "-log10 adj. p-value",
        showline = FALSE,
        showgrid = TRUE,
        zeroline = FALSE
      ),
      hovermode = 'closest',
      shapes = list(
        vline(sig_log2fc),
        vline(-sig_log2fc),
        hline(-log10(sig_pval))
      )
    )
}