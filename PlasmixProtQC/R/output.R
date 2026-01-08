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
  
  # 只有 expr_dt 和 metadata
  expr_dt <- data_list$expr_dt
  meta    <- data_list$metadata
  
  # 2. Run the QC pipelines
  # 只传两个参数
  allmetrics_results <- qc_allmetrics(expr_dt, meta, output_dir, plot)
  
  # 提取结果表
  output_table <- allmetrics_results$output_table
  
  # 3. Save & Output
  if (!is.null(output_dir)) {
    output_dir_final <- file.path(output_dir, "conclusion_table.tsv")
    write.table(output_table, output_dir_final, sep = "\t", row.names = F, quote = F)
  }
  
  final_list <- list(
    results = allmetrics_results,
    conclusion = output_table
  )
  
  return(final_list)
}