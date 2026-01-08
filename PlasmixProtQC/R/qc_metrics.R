#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data --------------------------------
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample
  
  # Replace zero by NA -----------------------
  d[d == 0] <- NA
  
  # Statistics: number of features -----------
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)
  
  # Statistics: missing percentage -----------
  d_all_num <- nrow(d) * (ncol(d) - 1)
  d_missing_num <- length(which(is.na(d)))
  prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)
  
  # Check if replicates available ------------
  samples <- table(s)
  rep_samples <- samples[samples > 1]
  rep_num <- length(rep_samples)
  if (rep_num == 0) {
    stop("No replicates are available.")
  } else {
    # Calculating: absolute correlation ------
    d_mtx <- d[, 2:ncol(d)]
    d_cortest <- corr.test(d_mtx, method = "pearson", adjust = "fdr")
    d_pmtx <- d_cortest$p
    d_cormtx <- d_cortest$r
    d_cormtx[d_pmtx > 0.05] <- 0
    d_cordf <- melt(d_cormtx)
    d_cordf <- d_cordf[d_cordf$Var2 != d_cordf$Var1, ]
    d_cordf <- merge(d_cordf, m, by.x = "Var1", by.y = "library")
    d_cordf <- merge(d_cordf, m, by.x = "Var2", by.y = "library")
    d_cordf <- d_cordf[d_cordf$sample.x == d_cordf$sample.y, ]
    cor_value <- median(d_cordf$value)
    stat_acor <- round(cor_value, 3)
  }
  
  # Calculating: CV --------------------------
  d_long <- melt(d)
  # melt 可能会提示 "Using Feature as id variables"，这是正常的 info 消息
  
  d_long <- na.omit(d_long)
  if (length(d_long$value[d_long$value < 0])) {
    d_long$value <- 2^(d_long$value)
    message("There are negative values. Default to log2-transformed values.")
  }
  
  # 合并 metadata
  d_long <- merge(d_long, m, by.x = "variable", by.y = "library")
  
  # --- 修复报错的关键部分 ---
  # merge 后现在的 d_long 包含: variable, [Feature], value, sample, platform
  # 我们需要显式选择列，防止列数不匹配
  feature_col <- colnames(d)[1] # 自动获取第一列的列名 (e.g. "Feature")
  
  # 只保留计算 CV 需要的列
  d_long <- d_long[, c("variable", feature_col, "value", "sample")]
  
  # 现在可以安全重命名了
  colnames(d_long) <- c("library", "feature", "value", "sample")
  # ------------------------
  
  d_cv <- aggregate(
    value ~ feature + sample,
    data = d_long,
    FUN = function(x) sd(x) / mean(x)
  )
  stat_cv <- round(median(d_cv$value, na.rm = T) * 100, 3)
  
  # Output -----------------------------------
  stat_all <- c(stat_num, stat_missing, stat_acor, stat_cv)
  
  return(stat_all)
}

#' Get color mapping for samples
#' @param samples Vector of sample names
#' @return Named vector of colors
#' @export
get_sample_colors <- function(samples) {
  # --- 修改 1: 更新颜色配色方案 ---
  # M='#3171b8', Y='#53a949', P='#6d3390', X='#ffc65d', F='#ae231c'
  color_palette <- c(
    "M" = "#3171b8", # Blue
    "Y" = "#53a949", # Green
    "P" = "#6d3390", # Purple (Reference)
    "X" = "#ffc65d", # Yellow
    "F" = "#ae231c"  # Red
  )
  
  # Return colors only for existing samples
  # 如果遇到未知样本名，默认给灰色
  available_colors <- color_palette[samples]
  if(any(is.na(available_colors))) {
    available_colors[is.na(available_colors)] <- "gray50"
  }
  return(available_colors)
}

#' Calculating SNR value; Plotting a PCA panel
#' @param expr_dt A expression profile (at protein level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @import stats
#' @import utils
#' @importFrom rlang :=
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 ggsave
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # Load data --------------------------------------
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])
  
  # Replace NA by zero -----------------------------
  expr_df[is.na(expr_df)] <- 0
  
  # Label the grouping info ------------------------
  ids <- colnames(expr_df)
  group <- meta_dt$sample
  ids_group_mat <- data.table(id = ids, group = group)
  
  # PCA --------------------------------------------
  expr_df_t <- t(expr_df)
  
  # Remove constant/zero variance features before PCA
  col_vars <- apply(expr_df_t, 2, var, na.rm = TRUE)
  zero_var_cols <- which(col_vars == 0 | is.na(col_vars))
  
  n_features_original <- ncol(expr_df_t)
  n_features_removed <- length(zero_var_cols)
  
  if (n_features_removed > 0) {
    message(sprintf(
      "Removing %d features with zero variance before PCA analysis.",
      n_features_removed
    ))
    expr_df_t <- expr_df_t[, -zero_var_cols, drop = FALSE]
  }
  
  n_features_used <- ncol(expr_df_t)
  
  # Check if enough features remain
  if (n_features_used < 2) {
    stop("Not enough features with non-zero variance for PCA analysis. At least 2 features required.")
  }
  
  pca_prcomp <- prcomp(expr_df_t, retx = T, scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$sample_id <- rownames(pcs)
  pcs$sample <- meta_dt$sample
  
  # Calculating: SNR -------------------------------
  dt_perc_pcs <- data.table(
    PCX = 1:nrow(pcs),
    Percent = summary(pca_prcomp)$importance[2, ],
    AccumPercent = summary(pca_prcomp)$importance[3, ]
  )
  
  dt_dist <- data.table(
    id_a = rep(ids, each = length(ids)),
    id_b = rep(ids, time = length(ids))
  )
  
  dt_dist$group_a <- ids_group_mat[match(dt_dist$id_a, ids_group_mat$id)]$group
  dt_dist$group_b <- ids_group_mat[match(dt_dist$id_b, ids_group_mat$id)]$group
  
  dt_dist[, type := ifelse(id_a == id_b, "Same",
                           ifelse(group_a == group_b, "Intra", "Inter")
  )]
  
  dt_dist[, dist := (dt_perc_pcs[1]$Percent * (pcs[id_a, 1] - pcs[id_b, 1])^2 +
                       dt_perc_pcs[2]$Percent * (pcs[id_a, 2] - pcs[id_b, 2])^2)]
  
  dt_dist_stats <- dt_dist[, list(avg_dist = mean(dist)), by = list(type)]
  setkey(dt_dist_stats, type)
  signoise <- dt_dist_stats["Inter"]$avg_dist / dt_dist_stats["Intra"]$avg_dist
  signoise_db <- round(10 * log10(signoise), 3)
  
  # Plot -------------------------------------------
  p <- NULL
  if (plot) {
    unique_samples <- unique(meta_dt$sample)
    # 调用更新后的 get_sample_colors
    colors_custom <- get_sample_colors(unique_samples)
    
    text_custom_theme <- element_text(
      size = 16,
      face = "plain",
      color = "black",
      hjust = 0.5
    )
    scale_axis_x <- c(min(pcs$PC1), max(pcs$PC1))
    scale_axis_y <- c(min(pcs$PC2), max(pcs$PC2))
    
    pc1_prop <- summary(pca_prcomp)$importance[2, 1]
    pc2_prop <- summary(pca_prcomp)$importance[2, 2]
    text_axis_x <- sprintf("PC1(%.2f%%)", pc1_prop * 100)
    text_axis_y <- sprintf("PC2(%.2f%%)", pc2_prop * 100)
    limit_x <- c(1.1 * scale_axis_x[1], 1.1 * scale_axis_x[2])
    limit_y <- c(1.1 * scale_axis_y[1], 1.1 * scale_axis_y[2])
    
    p_title <- paste("SNR = ", signoise_db, sep = "")
    p_subtitle <- paste("(Proteins used in PCA = ", n_features_used,
                        "/", n_features_original, ")",
                        sep = ""
    )
    
    p <- ggplot(pcs, aes(x = .data$PC1, y = .data$PC2)) +
      geom_point(aes(color = sample), size = 8) +
      theme_few() +
      theme(
        plot.title = text_custom_theme,
        plot.subtitle = text_custom_theme,
        axis.title = text_custom_theme,
        axis.text = text_custom_theme,
        legend.title = text_custom_theme,
        legend.text = element_text(size = 16, color = "gray40")
      ) +
      labs(
        x = text_axis_x,
        y = text_axis_y,
        title = p_title,
        subtitle = p_subtitle
      ) +
      scale_color_manual(values = colors_custom) +
      scale_x_continuous(limits = limit_x) +
      scale_y_continuous(limits = limit_y) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }
  
  pc_num <- ncol(pcs)
  output <- data.table(pcs[, c((pc_num - 1):pc_num, 1:(pc_num - 2))])
  
  if (!is.null(output_dir)) {
    if (plot) {
      output_dir_final1 <- file.path(output_dir, "pca_plot.png")
      ggsave(output_dir_final1, p, width = 6, height = 5.5)
    }
    output_dir_final2 <- file.path(output_dir, "pca_table.tsv")
    write.table(output, output_dir_final2, sep = "\t", row.names = F)
  }
  
  return(list(table = output, SNR = signoise_db, snr_plot = p))
}

#' Analysis: differential expression
#' @param expr A expression table file (at peptide level)
#' @param group The grouping info
#' @import stats
#' @importFrom edgeR DGEList
#' @importFrom edgeR filterByExpr
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export
dep_analysis <- function(expr, group) {
  dge <- DGEList(counts = expr)
  design <- model.matrix(~group)
  
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot = F)
  fit <- lmFit(v, design)
  
  fit <- eBayes(fit)
  result <- topTable(fit, coef = ncol(design), sort.by = "logFC", number = Inf)
  
  # result$Sequence <- rownames(result)
  # result$Sequence.Number <- nrow(result)
  
  # --- 修改: 统一命名为 Feature ---
  result$Feature <- rownames(result)
  result$Feature.Number <- nrow(result)
  # ------------------------------
  
  result$Sample1 <- levels(group)[1]
  result$Sample2 <- levels(group)[2]
  result$Sample.Pair <- paste(levels(group)[2], levels(group)[1], sep = "/")
  
  return(result)
}

#' Calculating RC value; Plotting a scatterplot
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @param show_sample_pairs if True, samples in plot will be labeled.
#' @import stats
#' @import utils
#' @importFrom rlang .data
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggthemes theme_few
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggsave
#' @export
qc_cor <- function(expr_dt, meta_dt,
                   output_dir = NULL, plot = FALSE, show_sample_pairs = TRUE) {
  
  # 1. 确定当前数据的平台
  current_platform <- unique(meta_dt$platform)
  if (length(current_platform) > 1) {
    stop("Error: Multiple platforms detected in metadata. Please process one platform at a time.")
  }
  if (length(current_platform) == 0) {
    stop("Error: Platform info not found in metadata.")
  }
  
  # 2. 加载内置数据集
  utils::data("reference_dataset", package = "PlasmixProtQC", envir = environment())
  ref_dt <- reference_dataset
  
  if (!"platform" %in% colnames(ref_dt)) {
    stop("Internal reference dataset is missing 'platform' column.")
  }
  
  # 3. 根据平台过滤
  ref_dt <- ref_dt[ref_dt$platform == current_platform, ]
  
  if (nrow(ref_dt) == 0) {
    stop(paste0("No reference data found for platform: ", current_platform))
  }
  
  # 4. 数据预处理
  expr_ncol <- ncol(expr_dt)
  expr_df <- data.frame(expr_dt[, 2:expr_ncol], row.names = expr_dt[, 1])
  expr_matrix <- as.matrix(expr_df)
  expr_matrix <- expr_matrix[, meta_dt$library]
  
  # Replace NA by zero
  expr_matrix[is.na(expr_matrix)] <- 0
  
  # Check negative values
  if (length(expr_matrix[expr_matrix < 0])) {
    expr_matrix <- 2^(expr_matrix)
    message("There are negative values. Default to log2-transformed values.")
  }
  
  # Check grouping info
  samples <- as.character(unique(meta_dt$sample))
  if (length(samples) < 2) {
    stop("At least 2 different sample types are required for correlation analysis.")
  }
  
  # --- 修改 2: Reference sample logic (P as Denominator) ---
  # 优先检查 "P" 是否存在
  check_ref <- "P" %in% samples
  if (!check_ref) {
    # 如果没有 P，尝试回退到 D6 (旧逻辑兼容)，或者报错
    # 这里我们改为警告并使用第一个样本作为 Reference (通用逻辑)
    warning("Reference sample 'P' is not available. Using the first sample as reference.")
    reference_sample <- samples[1]
    samples <- c(reference_sample, samples[!samples %in% reference_sample])
  } else {
    reference_sample <- "P"
    samples <- c("P", samples[!samples %in% "P"])
  }
  
  pair_num <- length(samples) - 1
  
  # Differential expression analysis
  result_final <- c()
  for (j in 2:(pair_num + 1)) {
    # 这里的 Sample Pair 命名规则为 Test/Reference，即 X/P
    sample_pair_name <- paste(samples[j], reference_sample, sep = "/")
    
    ref_tmp <- ref_dt[ref_dt$sample_pair %in% sample_pair_name, ]
    
    col1 <- which(meta_dt$sample %in% samples[j])
    col2 <- which(meta_dt$sample %in% reference_sample)
    
    e_tmp <- expr_matrix[, c(col1, col2)]
    e_tmp <- e_tmp[apply(e_tmp, 1, function(x) length(which(x == 0)) < min(length(col1), length(col2))), ]
    
    expr_grouped <- e_tmp[rownames(e_tmp) %in% ref_tmp$feature, ]
    
    # 在做差异分析时，level 的顺序决定了分母。
    # levels = c(Ref, Test) -> logFC = Test - Ref = Test/Ref
    sample_pairs <- factor(
      x = rep(c(samples[j], reference_sample), times = c(length(col1), length(col2))),
      levels = c(reference_sample, samples[j]), # 第一个是分母 (P)
      ordered = T
    )
    
    if (nrow(expr_grouped) > 0) {
      result_tmp <- dep_analysis(expr = expr_grouped, group = sample_pairs)
      result_tmp <- result_tmp[result_tmp$adj.P.Val < 0.05, ]
      result_final <- rbind(result_final, result_tmp)
    }
  }
  
  # Calculating: RC
  if (is.null(result_final) || nrow(result_final) == 0) {
    return(list(DEPs = NULL, logfc = NULL, COR = NA, cor_plot = NULL))
  }
  
  result_final <- as.data.table(result_final)
  
  # --- 修改: 适配 dep_analysis 的新列名 Feature ---
  result_trim <- data.frame(
    feature = result_final$Feature,       # 这里改成了 Feature
    sample_pair = result_final$Sample.Pair,
    logFC.Test = result_final$logFC
  )
  
  result_trim$name <- paste(result_trim$feature, result_trim$sample_pair)
  ref_dt$name <- paste(ref_dt$feature, ref_dt$sample_pair)
  
  result_withref <- merge(result_trim, ref_dt, by = "name")
  
  df_test <- data.frame(
    Name = result_withref$name,
    Feature = result_withref$feature.x,
    Sample.Pair = result_withref$sample_pair.x,
    logFC.Test = result_withref$logFC.Test,
    logFC.Reference = result_withref$assigned_value
  )
  
  cor_value <- cor(x = df_test$logFC.Test, y = df_test$logFC.Reference)
  cor_value <- round(cor_value, 3)
  
  # Plotting
  p <- NULL
  if (plot) {
    text_custom_theme <- element_text(
      size = 16,
      face = "plain",
      color = "black",
      hjust = 0.5
    )
    
    scale_axis_r <- c(min(df_test$logFC.Reference), max(df_test$logFC.Reference))
    scale_axis_t <- c(min(df_test$logFC.Test), max(df_test$logFC.Test))
    limit <- max(abs(c(scale_axis_r, scale_axis_t)))
    limit_axis <- c(-limit, limit)
    
    plot_title <- paste("RC = ", cor_value, sep = "")
    plot_subtitle <- paste("(Number of features = ", nrow(df_test), ")", sep = "")
    
    # --- 修改 3: Plot Color Mapping with new Palette ---
    unique_comps <- unique(df_test$Sample.Pair)
    
    # 重新获取定义的颜色
    base_palette <- get_sample_colors(names(get_sample_colors(c("X")))) # Hacky way to get all palette, or redefine locally
    # 更简单的方法：直接硬编码我们知道的颜色，确保这里和 get_sample_colors 一致
    base_palette <- c(
      "M" = "#3171b8",
      "Y" = "#53a949",
      "P" = "#6d3390",
      "X" = "#ffc65d",
      "F" = "#ae231c"
    )
    
    pair_colors <- sapply(unique_comps, function(comp_name) {
      # 提取分子 (Test Sample)，例如 "X/P" -> "X"
      numerator_sample <- strsplit(as.character(comp_name), "/")[[1]][1]
      
      if (numerator_sample %in% names(base_palette)) {
        return(base_palette[[numerator_sample]])
      } else {
        return("gray")
      }
    })
    names(pair_colors) <- unique_comps
    
    p <- ggplot(df_test, aes(x = .data$logFC.Reference, y = .data$logFC.Test)) +
      theme_few() +
      theme(
        plot.title = text_custom_theme,
        plot.subtitle = text_custom_theme,
        axis.title = text_custom_theme,
        axis.text = text_custom_theme,
        legend.title = text_custom_theme,
        legend.text = element_text(size = 16, color = "gray40")
      ) +
      labs(
        y = "log2FC (Test Datasets)",
        x = "log2FC (Reference Datasets)",
        title = plot_title,
        subtitle = plot_subtitle
      ) +
      coord_fixed(xlim = limit_axis, ylim = limit_axis)
    
    if (show_sample_pairs == T) {
      p <- p +
        geom_point(aes(color = .data$Sample.Pair), size = 2.5, alpha = .5) +
        scale_color_manual(values = pair_colors, name = "Sample Pair")
    } else {
      p <- p + geom_point(color = "steelblue4", size = 2.5, alpha = .1)
    }
  }
  
  output_list <- list(
    DEPs = result_final,
    logfc = df_test,
    COR = cor_value,
    cor_plot = p
  )
  
  return(output_list)
}

#' Calculating Recall of nominal characteristics
#' @param expr_dt A expression table file (at peptide level)
#' @param meta_dt A metadata file (to map library to sample type)
#' @export
qc_recall <- function(expr_dt, meta_dt) {
  
  current_platform <- unique(meta_dt$platform)
  if (length(current_platform) > 1) {
    stop("Error: Multiple platforms detected in metadata.")
  }
  if (length(current_platform) == 0) {
    stop("Error: Platform info not found in metadata.")
  }
  
  utils::data("reference_dataset_quali", package = "PlasmixProtQC", envir = environment())
  ref_dt <- reference_dataset_quali
  
  if (!"platform" %in% colnames(ref_dt)) {
    stop("Internal reference dataset (quali) is missing 'platform' column.")
  }
  
  ref_dt <- ref_dt[ref_dt$platform == current_platform, ]
  
  if (nrow(ref_dt) == 0) {
    warning(paste0("No reference data found for platform: ", current_platform, ". Recall will be NA."))
    return(NA)
  }
  
  # 样本名清洗 (依然保留以防万一)
  ref_dt$sample_clean <- gsub("Quartet ", "", ref_dt$sample)
  
  sample_types <- unique(meta_dt$sample)
  ref_sample_types <- unique(ref_dt$sample_clean)
  common_types <- intersect(sample_types, ref_sample_types)
  
  if (length(common_types) == 0) {
    warning("No common sample types found between metadata and reference dataset.")
    return(NA)
  }
  
  group_recalls <- c()
  
  for (st in common_types) {
    ref_seqs <- unique(ref_dt$feature[ref_dt$sample_clean == st])
    n_ref <- length(ref_seqs)
    
    if (n_ref == 0) next
    
    libs <- meta_dt$library[meta_dt$sample == st]
    valid_libs <- libs[libs %in% colnames(expr_dt)]
    
    if (length(valid_libs) == 0) next
    
    lib_recalls <- c()
    for (lib in valid_libs) {
      vals <- expr_dt[[lib]]
      detected_idx <- which(!is.na(vals) & vals > 0)
      
      det_seqs <- expr_dt[[1]][detected_idx]
      n_det_ref <- length(intersect(det_seqs, ref_seqs))
      
      rec <- n_det_ref / n_ref
      lib_recalls <- c(lib_recalls, rec)
    }
    
    if (length(lib_recalls) > 0) {
      group_recalls <- c(group_recalls, mean(lib_recalls))
    }
  }
  
  if (length(group_recalls) == 0) return(NA)
  
  final_recall <- mean(group_recalls)
  return(final_recall)
}