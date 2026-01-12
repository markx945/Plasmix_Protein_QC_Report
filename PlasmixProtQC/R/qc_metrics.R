#' Statistics for basic information
#' @param expr_dt A expression profile
#' @param meta_dt A metadata file
#' @import stats
#' @importFrom psych corr.test
#' @importFrom reshape2 melt
#' @importFrom ggplot2 geom_abline
#' @export
qc_info <- function(expr_dt, meta_dt) {
  # Load data
  m <- meta_dt
  d <- expr_dt
  s <- meta_dt$sample
  
  # Replace zero by NA
  d[d == 0] <- NA
  
  # 1. Feature Count
  uniq_pro <- unique(d[, 1])
  stat_num <- length(uniq_pro)
  
  # 2. Missing %
  d_vals <- d[, 2:ncol(d)]
  d_all_num <- nrow(d_vals) * ncol(d_vals)
  d_missing_num <- length(which(is.na(d_vals)))
  prop_missing <- 0
  if (d_all_num > 0) prop_missing <- d_missing_num * 100 / d_all_num
  stat_missing <- round(prop_missing, 3)
  
  # 3. Correlation & CV
  samples <- table(s)
  rep_samples <- samples[samples > 1]
  
  stat_acor <- NA
  stat_cv <- NA
  
  if (length(rep_samples) > 0) {
    # Correlation (Log data correlation is fine)
    d_mtx <- d[, 2:ncol(d)]
    common_libs <- intersect(colnames(d_mtx), m$library)
    if(length(common_libs) >= 2) {
      d_mtx <- d_mtx[, common_libs]
      cor_res <- cor(d_mtx, use = "pairwise.complete.obs", method = "pearson")
      stat_acor <- round(median(cor_res[upper.tri(cor_res)], na.rm = TRUE), 3)
    }
    
    # CV Calculation
    d_long <- melt(d, id.vars = colnames(d)[1]) 
    d_long <- na.omit(d_long)
    d_long <- merge(d_long, m, by.x = "variable", by.y = "library")
    
    # 强制假设输入数据是 Log2 的，先转回 Linear 再算 CV
    d_long$value_linear <- 2^(d_long$value)
    
    d_cv <- aggregate(
      value_linear ~ Feature + sample,
      data = d_long,
      FUN = function(x) sd(x) / mean(x)
    )
    stat_cv <- round(median(d_cv$value_linear, na.rm = T) * 100, 3)
  }
  
  return(c(stat_num, stat_missing, stat_acor, stat_cv))
}

#' Get color mapping
#' @export
get_sample_colors <- function(samples) {
  cols <- c("M"="#3171b8", "Y"="#53a949", "P"="#6d3390", "X"="#ffc65d", "F"="#ae231c")
  res <- cols[samples]
  res[is.na(res)] <- "gray50"
  return(res)
}

#' Calculating SNR
#' @import stats
#' @import utils
#' @importFrom data.table data.table setkey
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme labs guides guide_legend ggsave scale_x_continuous scale_y_continuous
#' @importFrom ggthemes theme_few
#' @export
qc_snr <- function(expr_dt, meta_dt, output_dir = NULL, plot = TRUE) {
  # 数据准备
  expr_df <- data.frame(expr_dt[, 2:ncol(expr_dt)], row.names = expr_dt[, 1])
  expr_df[is.na(expr_df)] <- 0 
  
  # 对齐 Metadata
  libs <- colnames(expr_df)
  meta_sub <- meta_dt[match(libs, meta_dt$library), ]
  group <- meta_sub$sample
  
  # --- 修改 1: PCA 缩放逻辑调整 ---
  # 检查平台信息
  current_platform <- unique(meta_sub$platform)
  
  # 默认不 Scale (SomaScan/Olink 都是相对定量单位)
  do_scale <- FALSE
  
  # 只有 DIA 平台才开启 Scale
  if(length(current_platform) > 0 && any(grepl("DIA", current_platform, ignore.case = TRUE))) {
    do_scale <- TRUE
  }
  
  # PCA
  expr_t <- t(expr_df)
  vars <- apply(expr_t, 2, var)
  expr_t <- expr_t[, vars > 0 & !is.na(vars), drop = FALSE]
  
  if(ncol(expr_t) < 2) stop("Too few features for PCA.")
  
  # 应用 Scale 参数
  pca <- prcomp(expr_t, scale. = do_scale)
  
  pcs <- as.data.frame(pca$x)
  pcs$sample <- group
  
  # 计算 SNR
  prop <- summary(pca)$importance[2, 1:2] 
  
  n <- nrow(pcs)
  dist_val <- matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      d2 <- prop[1]*(pcs[i,1]-pcs[j,1])^2 + prop[2]*(pcs[i,2]-pcs[j,2])^2
      dist_val[i,j] <- d2
    }
  }
  
  intra_dists <- c()
  inter_dists <- c()
  
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      if(group[i] == group[j]) {
        intra_dists <- c(intra_dists, dist_val[i,j])
      } else {
        inter_dists <- c(inter_dists, dist_val[i,j])
      }
    }
  }
  
  snr <- mean(inter_dists) / mean(intra_dists)
  snr_db <- round(10 * log10(snr), 2)
  
  # Plot
  p <- NULL
  if(plot) {
    cols <- get_sample_colors(unique(group))
    p <- ggplot(pcs, aes(x = PC1, y = PC2, color = sample)) +
      geom_point(size = 6) +
      scale_color_manual(values = cols) +
      theme_few() +
      labs(title = paste0("SNR = ", snr_db),
           x = sprintf("PC1 (%.1f%%)", prop[1]*100),
           y = sprintf("PC2 (%.1f%%)", prop[2]*100))
  }
  
  if(!is.null(output_dir)) {
    if(plot) ggsave(file.path(output_dir, "pca_plot.png"), p, width=6, height=5)
    write.table(pcs, file.path(output_dir, "pca_table.tsv"), sep="\t", quote=F)
  }
  
  return(list(SNR = snr_db, snr_plot = p))
}

#' Calculating RC
#' @importFrom utils data
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic labs coord_fixed
#' @export
qc_cor <- function(expr_dt, meta_dt, output_dir = NULL, plot = FALSE, show_sample_pairs = TRUE) {
  
  # 1. 准备参考数据
  plat <- unique(meta_dt$platform)
  
  # 尝试加载数据
  tryCatch({
    utils::data("reference_dataset", package = "PlasmixProtQC", envir = environment())
  }, error = function(e) {
    utils::data("reference_dataset", package = "protqc", envir = environment())
  })
  
  ref_dt <- reference_dataset[reference_dataset$platform == plat, ]
  if(nrow(ref_dt) == 0) stop("No reference data.")
  
  # 2. 准备表达矩阵 (Feature x Sample)
  mat <- as.matrix(expr_dt[, -1])
  rownames(mat) <- expr_dt[[1]]
  mat <- mat[, meta_dt$library] 
  
  # [修正 1] 删除这里的 mat[is.na(mat)] <- 0，保留 NA 状态！
  # mat[is.na(mat)] <- 0
  
  
  # --- 修改 2: RC 计算逻辑 (还原Log2 -> 算比值) ---
  
  # 2.1 将 Log2 数据还原为线性信号值 (2^x)
  # 无论什么平台，前提假设是用户输入了 Log2 转换后的数据
  mat_linear <- 2^mat
  
  # 3. 确定 Reference Sample (P)
  samps <- unique(meta_dt$sample)
  ref_samp <- "P"
  if(!"P" %in% samps) {
    ref_samp <- samps[1]
    warning("No 'P' sample. Using ", ref_samp)
  }
  test_samps <- setdiff(samps, ref_samp)
  
  # 4. 循环计算 Ratio
  res_all <- data.frame()
  
  for(ts in test_samps) {
    # 提取列索引
    cols_test <- which(meta_dt$sample == ts)
    cols_ref <- which(meta_dt$sample == ref_samp)
    
    # 提取线性数据子集
    dat_ref <- mat_linear[, cols_ref, drop=FALSE]
    dat_test <- mat_linear[, cols_test, drop=FALSE]
    
    # [修正 2] 过滤逻辑升级：Reference 组必须有有效值
    # 计算每行非 NA 的数量
    count_ref <- rowSums(!is.na(dat_ref))
    count_test <- rowSums(!is.na(dat_test))
    
    # 规则：Ref 组至少有 1 个有效值，Test 组至少有 1 个有效值
    # 如果 Reference 是全 NA，Ratio = Test / NaN = NaN，这种必须丢弃
    keep_rows <- (count_ref > 0) & (count_test > 0)
    
    # dat_ref <- dat_ref[keep_rows, , drop=FALSE]
    # dat_test <- dat_test[keep_rows, , drop=FALSE]
    # 
    # # 过滤掉 全0/全NA 的行
    # keep_rows <- (rowSums(dat_ref > 0, na.rm=TRUE) > 0) & (rowSums(dat_test > 0, na.rm=TRUE) > 0)
    # 
    # dat_ref <- dat_ref[keep_rows, , drop=FALSE]
    # dat_test <- dat_test[keep_rows, , drop=FALSE]
    
    current_features <- rownames(dat_ref)
    
    # 与参考集 Feature 取交集
    pair_name <- paste(ts, ref_samp, sep="/")
    feats_ref <- ref_dt$feature[ref_dt$sample_pair == pair_name]
    common_feats <- intersect(current_features, feats_ref)
    
    if(length(common_feats) < 3) next
    
    # 筛选特征
    dat_ref <- dat_ref[common_feats, , drop=FALSE]
    dat_test <- dat_test[common_feats, , drop=FALSE]
    
    # *** 核心计算: 原始信号值相比 (Mean Test / Mean Ref) ***
    # 计算行平均值 (Arithmetic Mean of Linear Signal)
    mean_ref <- rowMeans(dat_ref, na.rm = TRUE)
    mean_test <- rowMeans(dat_test, na.rm = TRUE)
    
    # 直接相除得到 Ratio
    ratio_val <- mean_test / mean_ref
    
    # 构建临时结果
    df_tmp <- data.frame(
      Feature = common_feats,
      Sample.Pair = pair_name,
      # 依然沿用 logFC.Test 这个列名以保持兼容性，但实际存储的是 Linear Ratio
      logFC.Test = ratio_val 
    )
    
    res_all <- rbind(res_all, df_tmp)
  }
  
  if(nrow(res_all) == 0) return(list(COR = NA, cor_plot = NULL))
  
  # 5. 合并参考值
  res_all$key <- paste(res_all$Feature, res_all$Sample.Pair)
  ref_dt$key <- paste(ref_dt$feature, ref_dt$sample_pair)
  
  merged <- merge(res_all, ref_dt, by="key")
  
  # # --- 升级版离群点检测 (Log-Space Outlier Detection) ---
  # # 1. 转换到 Log2 空间进行检测 (这样分布是对称的，IQR 更准确)
  # #    加极小值防止 log(0)
  # log2_vals <- log2(merged$logFC.Test + 1e-6)
  # 
  # # 2. 计算 Log 空间的 IQR
  # Q1 <- quantile(log2_vals, 0.25, na.rm = TRUE)
  # Q3 <- quantile(log2_vals, 0.75, na.rm = TRUE)
  # IQR_val <- Q3 - Q1
  # 
  # # 3. 设定 Log 空间的阈值
  # # 这里系数 1.5 是标准，如果想更严格可以设为 1.2
  # upper_log <- Q3 + 1.5 * IQR_val
  # lower_log <- Q1 - 1.5 * IQR_val
  # 
  # # 4. 找出离群点
  # is_outlier <- (log2_vals > upper_log) | (log2_vals < lower_log)
  # 
  # # *额外硬阈值*：Ratio > 100 或 < 0.01 无论如何都算离群 (生物学上极罕见)
  # is_extreme <- (merged$logFC.Test > 100) | (merged$logFC.Test < 0.01)
  # final_mask <- is_outlier | is_extreme
  # 
  # n_outliers <- sum(final_mask, na.rm = TRUE)
  # 
  # if (n_outliers > 0) {
  #   message(sprintf("Detected %d outliers in RC calculation (Log-IQR method). Removing them.", n_outliers))
  #   
  #   # 打印前几个被剔除的极端值，方便用户 debug
  #   removed_examples <- head(merged[final_mask, ], 3)
  #   message("Top removed outliers:")
  #   print(removed_examples[, c("Feature", "Sample.Pair", "logFC.Test")])
  #   
  #   merged_clean <- merged[!final_mask, ]
  # } else {
  #   merged_clean <- merged
  # }
  # ----------------------------------------------------
  
  # 计算相关性 (Test Ratio vs Reference Ratio)
  cor_val <- cor(merged$logFC.Test, merged$assigned_value, use="complete.obs")
  cor_val <- round(cor_val, 3)
  
  # Plot
  p <- NULL
  if(plot) {
    # [修复] 安全的颜色映射逻辑，防止下标越界
    unique_pairs <- unique(merged$Sample.Pair)
    
    # 1. 提取所有涉及的分子样本 (Test Samples)
    numerators <- unique(sapply(strsplit(as.character(unique_pairs), "/"), `[`, 1))
    
    # 2. 获取这些样本的标准颜色
    sample_cols_base <- get_sample_colors(numerators) # named vector
    
    # 3. 构建 Pair -> Color 的映射向量
    pair_colors_vec <- c()
    for(pair in unique_pairs) {
      num <- strsplit(as.character(pair), "/")[[1]][1]
      # 安全查找：如果 num 存在于 sample_cols_base 的名字中
      if(num %in% names(sample_cols_base)) {
        pair_colors_vec[pair] <- sample_cols_base[num]
      } else {
        pair_colors_vec[pair] <- "gray"
      }
    }
    
    p <- ggplot(merged, aes(x = assigned_value, y = logFC.Test)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(aes(color = Sample.Pair), alpha = 0.6, size = 2) +
      scale_color_manual(values = pair_colors_vec) +
      theme_few() +
      coord_fixed(xlim = c(0, 5), ylim = c(0, 5)) +
      scale_x_continuous(breaks = 0:5) +
      scale_y_continuous(breaks = 0:5) +
      labs(title = paste0("PCC = ", cor_val),
           x = "Reference Value",
           y = "Test Value")
  }
  
  return(list(COR = cor_val, cor_plot = p, logfc = merged))
}

#' Recall
#' @export
qc_recall <- function(expr_dt, meta_dt) {
  plat <- unique(meta_dt$platform)
  utils::data("reference_dataset_quali", package = "PlasmixProtQC", envir = environment())
  ref <- reference_dataset_quali[reference_dataset_quali$platform == plat, ]
  if(nrow(ref) == 0) return(NA)
  
  samps <- unique(meta_dt$sample)
  common <- intersect(samps, gsub("Quartet ", "", unique(ref$sample)))
  
  if(length(common) == 0) return(NA)
  
  vals <- c()
  for(s in common) {
    feats_ref <- ref$feature[gsub("Quartet ", "", ref$sample) == s]
    libs <- meta_dt$library[meta_dt$sample == s]
    
    for(l in libs) {
      det <- expr_dt[[1]][!is.na(expr_dt[[l]]) & expr_dt[[l]] > 0]
      recall <- length(intersect(det, feats_ref)) / length(feats_ref)
      vals <- c(vals, recall)
    }
  }
  return(mean(vals))
}