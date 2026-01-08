devtools::load_all()

devtools::document()

# 2. 定义文件路径
# 请确保这两个文件在你的工作目录，或者填写绝对路径
expr_file <- "./test_data/plasmix_somascan_test_expr.txt"
meta_file <- "./test_data/plasmix_somascan_test_meta.txt"

output_dir <- './test_output/'

tryCatch({
  qc_res <- qc_conclusion(
    exp_path = expr_file,
    meta_path = meta_file,
    output_dir = out_dir
  )
  
  # 4. 查看结果
  message("\n=== 分析完成 ===")
  print("结论表:")
  print(qc_res$conclusion)
  
  # 如果想直接在 RStudio 里看图，可以运行：
  # print(qc_res$results$snr_results$snr_plot)  # 查看 SNR PCA 图
  # print(qc_res$results$cor_results$cor_plot)  # 查看 RC 相关性图
  
}, error = function(e) {
  message("运行报错: ", e$message)
})

qc_res
