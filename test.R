# 1. 加载 PlasmixProtQC 包
library(PlasmixProtQC)

# 2. 定义文件路径
# 请确保这两个文件在你的工作目录，或者填写绝对路径
expr_file <- system.file("extdata", "plasmix_somascan_test_expr.csv", package = "PlasmixProtQC")
meta_file <- system.file("extdata", "plasmix_somascan_test_meta.csv", package = "PlasmixProtQC")

# olink
# expr_file <- system.file("extdata", "plasmix_olink_test_expr.csv", package = "PlasmixProtQC")
# meta_file <- system.file("extdata", "plasmix_olink_test_meta.csv", package = "PlasmixProtQC")

# DIA
# expr_file <- system.file("extdata", "plasmix_dia_test_expr.csv", package = "PlasmixProtQC")
# meta_file <- system.file("extdata", "plasmix_dia_test_meta.csv", package = "PlasmixProtQC")

template <- system.file("extdata", "Plasmix_template.docx", package = "PlasmixProtQC")

qc_res <- qc_conclusion(exp_path = expr_file, meta_path = meta_file)

# 4. 查看结果
# message("\n=== 分析完成 ===")
# print("结论表:")
# print(qc_res$conclusion)

# 如果想直接在 RStudio 里看图，可以运行：
# print(qc_res$results$snr_results$snr_plot)  # 查看 SNR PCA 图
# print(qc_res$results$cor_results$cor_plot)  # 查看 RC 相关性图

# 5. 生成报告
# report_template 可选，如果你想复用那个 docx 里的字体样式，可以传入路径
generate_protein_report(
  qc_result = qc_res,
  report_template = template,
  # report_dir = "./output_dir", # 可选，指定输出目录
  # report_name = "我的测试报告.docx", # 可选，指定报告文件名
  # batch_name = "20260108测试批次" # 可选，指定表格里的批次名
)


# # 1. 提取完整的结果表
# res_table <- qc_res$detail_results$cor_results$logfc
#
# # 2. 筛选 logFC.Test (实际是 Ratio) 大于 100 的行
# # 这里的列名虽然叫 logFC.Test，但根据我们修改后的逻辑，里面存的是 Linear Ratio
# outliers_100 <- res_table[res_table$logFC.Test > 100, ]
#
# # 3. 按数值从大到小排序，查看前几行
# outliers_sorted <- outliers_100[order(outliers_100$logFC.Test, decreasing = TRUE), ]
#
# print(">>> Ratio > 100 的数据统计：")
# print(nrow(outliers_sorted)) # 看看有多少个
#
# print(">>> Ratio 最大的前 10 个：")
# print(head(outliers_sorted, 10))
