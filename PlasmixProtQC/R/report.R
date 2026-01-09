#' @title Generate Protein Quality Report (Chinese Version)
#' @description Generates a report following the ShenKang Plasmix template structure.
#' @param qc_result list, the output from qc_conclusion()
#' @param report_template character, path to a .docx file to use as style template (optional)
#' @param report_dir character, directory to save the report
#' @param report_name character, filename of the report
#' @param batch_name character, optional, manual batch name. If NULL, tries to infer from metadata.
#' @importFrom dplyr %>%
#' @importFrom flextable flextable theme_box align width bold bg color body_add_flextable set_header_labels
#' @importFrom officer read_docx body_add_par body_add_gg body_add_break body_add_fpar fpar fp_text prop_section
#' @export
generate_protein_report <- function(qc_result,
                                    report_template = NULL,
                                    report_dir = NULL,
                                    report_name = NULL,
                                    batch_name = NULL) {
  
  # --- 1. 参数检查与初始化 ---
  if (is.null(qc_result)) stop("qc_result is required.")
  
  if (is.null(report_dir)) {
    report_dir <- file.path(getwd(), "output")
    dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (is.null(report_name)) {
    report_name <- paste0("Plasmix_QC_Report_", format(Sys.Date(), "%Y%m%d"), ".docx")
  }
  output_file <- file.path(report_dir, report_name)
  
  # --- 2. 提取并处理数据 ---
  raw_table <- qc_result$qc_metrics_table
  meta      <- qc_result$input_meta
  
  # 尝试获取批次名 (优先用参数，其次查 metadata 的 study_id，最后用默认值)
  if (is.null(batch_name)) {
    if ("study_id" %in% colnames(meta) && length(unique(meta$study_id)) == 1) {
      batch_name <- unique(meta$study_id)[1]
    } else {
      batch_name <- paste0(format(Sys.Date(), "%Y%m%d"), "批次")
    }
  }
  
  # 提取数值 (根据英文 Metrics 名提取)
  # 注意：这里要匹配 qc_pipelines.R 中定义的英文名
  n_prot_val <- raw_table$Value[grep("Number of features", raw_table$Quality_Metrics)]
  recall_val <- raw_table$Value[grep("Recall", raw_table$Quality_Metrics)]
  snr_val    <- raw_table$Value[grep("Signal-to-Noise Ratio", raw_table$Quality_Metrics)]
  rc_val     <- raw_table$Value[grep("Relative Correlation", raw_table$Quality_Metrics)]
  
  # 转换为数值
  val_prot   <- as.numeric(n_prot_val)
  val_recall <- as.numeric(recall_val)
  val_snr    <- as.numeric(snr_val)
  val_rc     <- as.numeric(rc_val)
  
  # 处理 NA
  if(length(val_prot)==0) val_prot <- NA
  if(length(val_recall)==0) val_recall <- NA
  if(length(val_snr)==0) val_snr <- NA
  if(length(val_rc)==0) val_rc <- NA
  
  # --- 3. 质量判定逻辑 (用于最后一列) ---
  # 标准: Protein >= 1000, Recall >= 0.90, SNR >= 5, RC >= 0.80
  pass_prot   <- !is.na(val_prot) && val_prot >= 1000
  pass_recall <- !is.na(val_recall) && val_recall >= 0.90
  pass_snr    <- !is.na(val_snr) && val_snr >= 5
  pass_rc     <- !is.na(val_rc) && val_rc >= 0.80
  
  is_all_pass <- pass_prot && pass_recall && pass_snr && pass_rc
  overall_quality <- ifelse(is_all_pass, "Pass", "Fail")
  
  # 格式化数值用于显示 (保留两位小数)
  txt_prot   <- format(val_prot, big.mark = ",")
  txt_recall <- sprintf("%.2f", val_recall)
  txt_snr    <- sprintf("%.2f", val_snr)
  txt_rc     <- sprintf("%.2f", val_rc)
  
  # --- 4. 构建表格数据 ---
  # 定义列名
  cols <- c("批次", "蛋白鉴定数", "标称特性召回率", "信噪比", "相对相关性", "整体质量")
  
  # 第一行：推荐标准
  row_std <- c("推荐质量标准", "≥1,000", "≥0.90", "≥5", "≥0.80", "全部通过")
  
  # 第二行：实际数据
  row_dat <- c(batch_name, txt_prot, txt_recall, txt_snr, txt_rc, overall_quality)
  
  df_rep <- data.frame(rbind(row_std, row_dat), stringsAsFactors = FALSE)
  colnames(df_rep) <- cols
  
  # 创建 Flextable
  ft <- flextable(df_rep) %>%
    theme_box() %>%
    align(align = "center", part = "all") %>%
    width(width = 1.2) %>%
    # 表头加粗背景灰
    bold(part = "header") %>%
    bg(part = "header", bg = "#EFEFEF") %>%
    # 第一行(标准)字体颜色变灰
    color(i = 1, color = "gray40") %>%
    # 第二行(数据)加粗
    bold(i = 2) %>%
    # 根据结果给"整体质量"上色
    color(i = 2, j = 6, color = ifelse(overall_quality == "Pass", "green", "red"))
  
  # --- 5. 文档文本内容 (来自提供的模板) ---
  txt_summary <- "本报告基于多项组学关键质量控制指标，总结了 Plasmix血浆参考物质所生成数据的质量情况。质量控制流程从用户输入血浆蛋白组表达矩阵开始，分别计算每批次的蛋白鉴定数（Number of identified proteins）、标称特性召回率 (Recall of nominal characteristics)、 信噪比 (Signal-to-Noise Ratio, SNR)、与参考数据集的相对相关性 (Relative Correlation with Reference Datasets, RC)及整体质量判断。"
  
  txt_def_title <- "质量控制指标"
  txt_def_1 <- "蛋白鉴定数（Number of identified proteins）：该指标反映在蛋白组检测中成功鉴定并映射到基因符号的蛋白数量。通常期望获得尽可能多的蛋白特征，以支持后续的生物学分析。"
  txt_def_2 <- "标称特性召回率（Recall of nominal characteristics）：定义为标称特性蛋白在测试数据集中被成功检测到的比例。"
  txt_def_3 <- "信噪比（Signal-to-Noise Ratio, SNR）：用于刻画某一检测平台、实验室或批次区分不同生物样本组之间内在生物学差异（“信号”）与同一样本组技术重复变异（“噪声”）的能力。SNR 越高，表明区分多组样本差异的能力越强。"
  txt_def_4 <- "与参考数据集的相对相关性（Relative Correlation with Reference Datasets, RC）：定义为在给定样本对之间，测试数据集中比值型表达水平与对应比值型参考数据集之间的 Pearson 相关系数，用于表征比值表达谱在数值层面的整体一致性趋势。为提高分析可靠性，在进行比值表达分析前，首先对每个样本组的技术重复取均值。差异倍数（fold change）采用 log2 转换。"
  
  txt_refs <- c(
    "1. Zheng Y, et al. Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials. Nature Biotechnology, 2024.",
    "2. Liu Y, et al. Harmonizing plasma proteomics data with the sample-to-reference ratio approach. BioRxiv (2026)",
    "3. 上海临床队列组学检测工作指引（征求意见稿）, 2025/11/26"
  )
  
  txt_disclaimer <- "本数据质量报告仅针对所评估的特定数据集提供分析结果，仅供信息参考之用。尽管已尽最大努力确保分析结果的准确性和可靠性，但本报告按“现状（AS IS）”提供，不附带任何形式的明示或暗示担保。报告作者及发布方不对基于本报告内容所采取的任何行动承担责任。本报告中的结论不应被视为对任何产品或流程质量的最终判定，也不应用于关键应用场景、商业决策或法规合规用途，除非经过专业核查和独立验证。对于分析结果的正确性、准确性、可靠性或适用性，不作任何明示或暗示的保证。"
  
  # --- 6. 生成 Word 文档 ---
  
  # 如果提供了模板，就加载模板；否则新建空白文档
  if (!is.null(report_template) && file.exists(report_template)) {
    doc <- read_docx(path = report_template)
  } else {
    doc <- read_docx()
  }
  
  # 文档构建顺序：
  # Title -> Summary -> Table -> Definitions -> Refs -> Disclaimer -> Plots
  
  # 标题
  doc <- doc %>%
    body_add_par("Plasmix血浆蛋白组质量报告", style = "heading 1") %>%
    body_add_par("", style = "Normal") # 空行
  
  # 摘要
  doc <- doc %>%
    body_add_par("摘要", style = "heading 2") %>%
    body_add_par(txt_summary, style = "Normal") %>%
    body_add_par("", style = "Normal")
  
  # 表格
  doc <- doc %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")
  
  # 定义
  doc <- doc %>%
    body_add_par(txt_def_title, style = "heading 2") %>%
    body_add_par(txt_def_1, style = "Normal") %>%
    body_add_par(txt_def_2, style = "Normal") %>%
    body_add_par(txt_def_3, style = "Normal") %>%
    body_add_par(txt_def_4, style = "Normal") %>%
    body_add_par("", style = "Normal")
  
  # 参考文献
  doc <- doc %>%
    body_add_par("参考文献", style = "heading 2")
  for(ref in txt_refs) {
    doc <- doc %>% body_add_par(ref, style = "Normal")
  }
  doc <- doc %>% body_add_par("", style = "Normal")
  
  # 免责声明
  doc <- doc %>%
    body_add_par("免责声明", style = "heading 2") %>%
    body_add_par(txt_disclaimer, style = "Normal") %>%
    body_add_break() # 分页
  
  # 图表部分
  # Signal-to-Noise Ratio
  doc <- doc %>%
    body_add_par("Signal-to-Noise Ratio", style = "heading 1") %>%
    body_add_gg(value = qc_result$snr_plot, style = "centered", width = 6, height = 5) %>%
    body_add_par("", style = "Normal")
  
  # Correlation
  doc <- doc %>%
    body_add_par("Correlation with Reference Datasets", style = "heading 1") %>%
    body_add_gg(value = qc_result$cor_plot, style = "centered", width = 6, height = 5) %>%
    body_add_par("", style = "Normal")
  
  # --- 7. 保存 ---
  print(doc, target = output_file)
  message(paste("Report generated successfully at:", output_file))
}