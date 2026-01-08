# test_simulation.R
library(tidyverse)
library(usethis)
library(devtools)

# ==============================================================================
# 1. 模拟并保存内置参考数据 (Reference Data) -> data/
# ==============================================================================
message(">>> 步骤1: 生成模拟参考数据集 (Reference)...")

# 设定平台名称
sim_platform <- "DIA"

# 1.1 设定特征池 (模拟 100 个 Uniprot ID)
# 假设这是参考标准品中“绝对存在”的 100 个蛋白
all_ref_features <- sprintf("P%05d", 1:100) # 生成 P00001 - P00100

samples <- c("X", "Y", "M", "F", "P")

# 1.2 生成定性参考 (Reference Quali)
# 这 100 个蛋白在所有样本中都应该被检测到
ref_quali <- expand.grid(sample = samples, feature = all_ref_features) %>%
  mutate(platform = sim_platform) %>%
  select(platform, sample, feature)

# 1.3 生成定量参考 (Reference Quant)
# 设定理论上的差异倍数 (LogFC)，以 P 为分母
# X/P = 2, Y/P = -2, M/P = 1, F/P = -1
ref_quant_list <- list()
pairs <- c("X/P", "Y/P", "M/P", "F/P")
expected_logfc <- c(2, -2, 1, -1)

# 设置随机种子保证可重复性
set.seed(123)

for(i in seq_along(pairs)) {
  # 为这 100 个蛋白生成参考值 (理论值 + 极小波动)
  vals <- rnorm(n = 100, mean = expected_logfc[i], sd = 0.05)
  df <- data.frame(
    platform = sim_platform,
    sample_pair = pairs[i],
    feature = all_ref_features,
    assigned_value = vals
  )
  ref_quant_list[[i]] <- df
}
ref_quant <- bind_rows(ref_quant_list)

# 保存到包的数据目录
usethis::use_data(reference_dataset_quali = ref_quali, 
                  reference_dataset = ref_quant, 
                  overwrite = TRUE)

message("   参考数据已生成并保存到 data/。")

# ==============================================================================
# 2. 重新加载包
# ==============================================================================
message(">>> 步骤2: 重新加载包以应用新数据...")
devtools::load_all()

# ==============================================================================
# 3. 模拟用户输入数据 (User Input Data)
# ==============================================================================
message(">>> 步骤3: 生成模拟用户输入数据...")

# 3.1 决定“检测到”了哪些蛋白 (关键步骤！)
# 场景：参考集有100个，我们模拟只检测到了其中的 90 个 (Recall = 0.9)
# 另外检测到了 20 个参考集以外的蛋白 (噪音，不影响 Recall 分子，只增加分母)

# A. 从 Reference 中挑选 90 个 (True Positives)
detected_tp <- sample(all_ref_features, 90)

# B. 生成 20 个干扰蛋白 (False Positives / Noise)
detected_noise <- sprintf("N%05d", 1:20)

# C. 合并作为表达矩阵的行名
input_features <- c(detected_tp, detected_noise)
message(sprintf("   模拟特征数: %d (其中 %d 个来自参考集, %d 个为噪音)", 
                length(input_features), length(detected_tp), length(detected_noise)))

# 3.2 生成 Metadata
reps <- 3
meta_df <- expand.grid(sample = samples, rep = 1:reps) %>%
  mutate(library = paste(sample, rep, sep = "_"),
         platform = sim_platform) %>% # 必须匹配 Reference 的平台
  select(library, sample, platform)

meta_df
# 3.3 生成表达矩阵
# 初始化矩阵 (行=特征, 列=样本)
expr_matrix <- matrix(nrow = length(input_features), ncol = nrow(meta_df))
rownames(expr_matrix) <- input_features
colnames(expr_matrix) <- meta_df$library

# 填充矩阵数值
base_expr <- 20 # 基础表达量

for(col in 1:ncol(expr_matrix)) {
  samp_type <- meta_df$sample[col]
  
  for(row in 1:nrow(expr_matrix)) {
    feat <- input_features[row]
    
    # 1. 确定理论差异 (LogFC)
    offset <- 0
    
    if(feat %in% detected_tp) {
      # 如果是参考蛋白，去 reference_dataset 里查表找 LogFC
      if(samp_type != "P") {
        pair_name <- paste0(samp_type, "/P")
        # 查找对应 LogFC
        # 注意：这里为了速度简化了查找，实际可以用 merge
        target_val <- ref_quant$assigned_value[ref_quant$feature == feat & ref_quant$sample_pair == pair_name]
        if(length(target_val) > 0) offset <- target_val
      }
    } else {
      # 如果是噪音蛋白，随机给一个偏移，与样本类型无关
      offset <- rnorm(1, 0, 1) 
    }
    
    # 2. 添加噪音 (影响 SNR)
    # 给参考蛋白加小噪音，给噪音蛋白加大噪音
    if(feat %in% detected_tp) {
      noise <- rnorm(1, mean = 0, sd = 0.2) 
    } else {
      noise <- rnorm(1, mean = 0, sd = 0.5)
    }
    
    expr_matrix[row, col] <- base_expr + offset + noise
  }
}

# 转换为数据框，第一列名为 Feature
expr_df <- as.data.frame(expr_matrix)
expr_df <- tibble::rownames_to_column(expr_df, var = "Feature")

# 保存为 CSV
write.csv(expr_df, "test_expression.csv", row.names = FALSE)
write.csv(meta_df, "test_metadata.csv", row.names = FALSE)

message("   测试文件已生成: test_expression.csv, test_metadata.csv")

# ==============================================================================
# 4. 运行 QC 流程
# ==============================================================================
message(">>> 步骤4: 运行 QC 分析...")

out_dir <- "test_output"
if(!dir.exists(out_dir)) dir.create(out_dir)

# 运行
qc_res <- qc_conclusion(
  exp_path = "test_expression.csv",
  meta_path = "test_metadata.csv",
  output_dir = out_dir
)

# 打印结果
print("--- 最终 QC 结果表 ---")
print(qc_res$conclusion)

message("\n预期结果:")
message("1. Recall 应该接近 0.90 (因为我们只挑选了 90% 的参考蛋白)")
message("2. SNR 应该较高 (>10)")
message("3. RC 应该接近 1.0 (因为数据是基于参考值生成的)")
message("测试完成！")