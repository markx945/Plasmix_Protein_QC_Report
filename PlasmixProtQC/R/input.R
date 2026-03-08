#' Input files
#' @param exp_path A file path of the expression table file
#' @param meta_path A file path of the metadata file
#' @import stats
#' @import utils
#' @importFrom data.table fread
#' @importFrom dplyr %>%
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @export
input_data <- function(exp_path, meta_path) {
  # Load data ------------------------------------------------
  # fread 能够自动识别 .txt, .csv, .tsv 的分隔符
  expr <- fread(exp_path)
  meta <- fread(meta_path)

  expr <- as.data.frame(expr)
  meta <- meta %>%
    rename_with(tolower)

  # Check metadata format ------------------------------------
  # 检查必要的列是否存在
  req_cols <- c("library", "sample", "platform")

  if ("name" %in% colnames(meta)) {
    meta <- meta %>% dplyr::rename(library = name)
  }

  if (!all(c("library", "sample", "platform") %in% colnames(meta))) {
    stop('The columns named "library", "sample", and "platform" are required in metadata.')
  }

  meta_final <- meta %>% select(library, sample, platform)

  if (any(duplicated(meta_final$library))) {
    stop("Duplicated library IDs in metadata.")
  }

  # Check expression data format -----------------------------

  # --- 新增逻辑: 兼容旧格式 (Type + Feature) ---
  # 如果前两列正好是 Type 和 Feature，则移除 Type 列
  if (colnames(expr)[1] == "Type" && colnames(expr)[2] == "Feature") {
    message("Detected legacy format with 'Type' column. Removing 'Type' and using 'Feature' as identifier.")
    expr <- expr[, -1] # 删除第一列
  }
  # -------------------------------------------

  # 1. 检查是否有重复列名
  if (length(which(duplicated(colnames(expr))))) {
    stop("Duplicated column names in data.")
  }

  # 2. 确定 Feature 列
  # 强制取第一列作为 Feature，并重命名为 "Feature"
  feature_col_name <- colnames(expr)[1]

  # 防止第一列不是 "Feature" 但也不是 "Type" 的情况 (比如叫 "Protein")
  if (feature_col_name != "Feature") {
    message(sprintf("Renaming the first column '%s' to 'Feature'.", feature_col_name))
    colnames(expr)[1] <- "Feature"
  }

  # 3. 检查 Library 列匹配
  # 排除第一列后，剩余列名应该是 Library ID
  expr_libs <- colnames(expr)[2:ncol(expr)]

  # 找出 Metadata 中有但在表达矩阵中没有的样本
  missing_libs <- setdiff(meta_final$library, expr_libs)
  if (length(missing_libs) > 0) {
    warning(sprintf("%d libraries in metadata are missing from expression data.", length(missing_libs)))
  }

  # 只保留 Metadata 中存在的列 (取交集)
  valid_libs <- intersect(expr_libs, meta_final$library)

  if (length(valid_libs) == 0) {
    stop("No common library IDs found between expression data and metadata.")
  }

  # 重组表达矩阵：第一列 Feature + 匹配的 Library 列
  expr_final <- expr[, c("Feature", valid_libs)]

  # 返回统一的 List
  data_list <- list(
    "expr_dt" = expr_final,
    "metadata" = meta_final
  )

  return(data_list)
}
