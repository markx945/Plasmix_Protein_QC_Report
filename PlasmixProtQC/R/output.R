#' Generating a table of conclusion
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @param output_dir A directory of the output file(s)
#' @param plot if True, a plot will be output.
#' @importFrom utils write.table
#' @export
qc_conclusion <- function(exp_path, meta_path, output_dir = NULL, plot = TRUE) {
  
  # 1. Load the input data
  data_list <- input_data(exp_path, meta_path)
  
  expr_dt <- data_list$expr_dt
  meta    <- data_list$metadata
  
  # 2. Run the QC pipelines
  allmetrics_results <- qc_allmetrics(expr_dt, meta, output_dir, plot)
  
  # 提取结果表
  output_table <- allmetrics_results$output_table
  
  # 3. Save & Output
  if (!is.null(output_dir)) {
    output_dir_final <- file.path(output_dir, "conclusion_table.tsv")
    write.table(output_table, output_dir_final, sep = "\t", row.names = F, quote = F)
  }
  
  # --- 修复逻辑：扁平化输出结构，方便 Report 调用 ---
  final_list <- list(
    # 核心表格 (Report 用这个名字)
    qc_metrics_table = output_table,
    
    # 核心图片 (Report 直接调用 qc_result$snr_plot)
    snr_plot = allmetrics_results$snr_results$snr_plot,
    cor_plot = allmetrics_results$cor_results$cor_plot, # 即 RC plot
    
    # 保留完整细节以备查
    detail_results = allmetrics_results,
    
    # 方便调试
    input_meta = meta
  )
  
  return(final_list)
}