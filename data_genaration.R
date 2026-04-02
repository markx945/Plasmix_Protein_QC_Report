# 安装必要的包 (如果尚未安装)
# install.packages(c("usethis", "readr", "dplyr"))

library(usethis)
library(readr)
library(dplyr)

# 1. 读取数据 (请确保文件名路径正确)
# 读取标称特性 (对应 reference_dataset_quali)
raw_quali <- read_csv("../标称特性.csv") # 你的文件名
# 读取特性量值 (对应 reference_dataset)
raw_quant <- read_csv("../特性量值.csv") # 你的文件名

raw_quali

# 2. 清洗与标准化 reference_dataset_quali (Recall用)
# 目标列名: platform, sample, peptide_sequence
# 假设你的CSV里列名可能是中文，这里我们需要重命名为英文
# 请根据你实际CSV的列名修改下面的 rename 逻辑
reference_dataset_quali <- raw_quali %>%
  # 假设 CSV 列名是: 平台, 样本, 肽段序列
  # rename(platform = `平台`, sample = `样本`, peptide_sequence = `肽段序列`) 
  # 如果你的列名已经是英文，确保包含以下列：
  select(platform, sample, feature) %>%
  mutate(sample = gsub("Quartet ", "", sample)) %>% # 清洗样本名，去掉可能的 "Quartet " 前缀
  distinct() # 去重

reference_dataset_quali

# 3. 清洗与标准化 reference_dataset (RC用)
# 目标列名: platform, sample_pair, peptide_sequence, logFC
reference_dataset <- raw_quant %>%
  # 假设 CSV 列名是: 平台, 样本对, 肽段序列, assigned_value
  # rename(platform = `平台`, sample_pair = `样本对`, feature = `肽段序列`, assigned_value = `assigned_value`)
  # 如果你的列名已经是英文，确保包含以下列：
  select(platform, sample_pair, feature, assigned_value) %>%
  distinct()

# 4. 保存为包的内置数据
# 这会将数据保存到 data/ 文件夹下，分别名为 reference_dataset.rda 和 reference_dataset_quali.rda
usethis::use_data(reference_dataset_quali, reference_dataset, overwrite = TRUE)

message("数据集已生成并保存到 data/ 目录")
