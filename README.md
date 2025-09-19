# 肠道微生物组分析流程 (Gut Microbiome Analysis Pipeline)

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![QIIME2](https://img.shields.io/badge/QIIME2-2023.9-green.svg)](https://qiime2.org/)
[![PICRUSt2](https://img.shields.io/badge/PICRUSt2-2.5.2-orange.svg)](https://github.com/picrust/picrust2)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

一个完整的肠道微生物组分析流程，从FASTQ原始数据到可视化HTML报告，提供健康评估、疾病风险预测和个性化建议。

## 🎯 主要功能

### 1. 数据预处理
- **FASTQ质控**：自动质量过滤和去噪
- **ASV聚类**：使用DADA2生成高精度的ASV表
- **物种注释**：基于SILVA 138数据库的准确分类
- **功能注释**：使用PICRUSt2预测微生物功能
- **批量处理**：支持多样本并行分析

### 2. 核心分析模块
- **多样性分析**：Alpha多样性（Shannon、Simpson、Chao1）和Beta多样性
- **肠型判定**：基于属水平的3种肠型分类（拟杆菌型、普氏菌型、瘤胃球菌型）
- **菌群评估**：43种核心细菌的健康评分
- **疾病风险**：14种疾病的风险评估
- **生物年龄**：基于菌群组成的肠道年龄预测
- **功能预测**：维生素合成、短链脂肪酸产生能力分析

### 3. 智能报告
- **HTML可视化报告**：交互式图表展示
- **中文注释**：完整的中文解读和建议
- **个性化建议**：饮食、益生菌、生活方式建议

## 📁 项目结构

```
16s_bighealth_pipe/
│
├── 📂 scripts/
│   ├── 📂 preprocessing/           # 数据预处理
│   │   └── backend_process.sh     # FASTQ批处理与功能注释主脚本
│   │
│   ├── 📂 analysis/                # 分析模块
│   │   ├── 1_basic_analysis.py    # 基础分析（多样性、B/F比值）
│   │   ├── 2_enterotype.py        # 肠型分析
│   │   ├── 3_bacteria_eval.py     # 菌群健康评估
│   │   ├── 4_disease_risk.py      # 疾病风险评估  
│   │   ├── 5_age_predict.py       # 生物年龄预测
│   │   ├── 6_functional_prediction.py  # 功能预测
│   │   └── 7_cn_annotation.py     # 中文注释生成
│   │
│   └── 📂 report/                  # 报告生成
│       ├── generate_report.py     # HTML报告生成器
│       ├── report_template.html   # 报告模板
│       ├── report_styles.css      # 样式文件
│       └── report_scripts.js      # 交互脚本
│
├── 📂 database/                    # 数据库文件
│   ├── silva_138_99_16S_338f_806r_classifier.qza  # SILVA物种数据库
│   ├── disease_associations.json  # 疾病关联数据
│   ├── core_bacteria_annotations.json  # 核心菌中文注释
│   ├── core_ec_translations.json      # EC酶中文翻译
│   ├── core_pathway_translations.json # 通路中文翻译
│   └── 📂 process_bak/            # 数据库生成脚本
│
├── 📄 pipeline.sh                  # 主控制脚本
├── 📄 requirements.txt            # Python依赖
├── 📄 README.md                   # 本文档
└── 📄 quick_start.md              # 快速开始指南
```

## 🚀 快速开始

### 环境要求

- **操作系统**: Linux/macOS (推荐Ubuntu 20.04+)
- **Python**: 3.8+
- **QIIME2**: 2023.9 (用于FASTQ预处理和物种注释)
- **PICRUSt2**: 2.5.2 (用于功能预测)
- **内存**: 最少8GB，推荐16GB+
- **存储**: 至少20GB可用空间

### 安装步骤

1. **克隆项目**
```bash
git clone https://github.com/chenziqing0111/16s_bighealth_pipe.git
cd 16s_bighealth_pipe
```

2. **安装Python依赖**
```bash
pip install -r requirements.txt
```

3. **安装QIIME2** (必需)
```bash
# 使用conda安装
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2023.9-py38-linux-conda.yml
conda env create -n qiime2-amplicon-2023.9 --file qiime2-amplicon-2023.9-py38-linux-conda.yml
conda activate qiime2-amplicon-2023.9
```

4. **安装PICRUSt2** (用于功能预测)
```bash
# 在QIIME2环境中安装
conda activate qiime2-amplicon-2023.9
conda install -c bioconda -c conda-forge picrust2=2.5.2
```

5. **下载必需数据库**
```bash
# SILVA数据库已包含在项目中
# 如需更新，可从以下地址下载：
# wget https://data.qiime2.org/2023.9/common/silva-138-99-nb-classifier.qza
# mv silva-138-99-nb-classifier.qza database/

# 下载gutMDisorder疾病关联数据库
wget http://bio-annotation.cn/gutMDisorder/download/gutMDisorder.csv
mv gutMDisorder.csv database/
```

### 使用方法

#### 方式1: 从FASTQ文件开始（完整流程）
```bash
# 基础用法
./pipeline.sh /path/to/fastq_dir /path/to/output -d database/silva_138_99_16S_338f_806r_classifier.qza

# 并行处理多个样本
./pipeline.sh /path/to/fastq_dir /path/to/output -d database/silva_138_99_16S_338f_806r_classifier.qza -j 8

# 带元数据的分析（包含年龄、性别等信息）
./pipeline.sh /path/to/fastq_dir /path/to/output -m metadata.tsv -d database/silva_138_99_16S_338f_806r_classifier.qza
```

#### 方式2: 从ASV表开始（跳过预处理）
```bash
# 如果你已有ASV表（例如来自其他软件）
./pipeline.sh merged_asv_table.tsv /path/to/output --skip-preprocessing

# 批量分析多个样本
./pipeline.sh merged_asv_table.tsv /path/to/output --skip-preprocessing -j 4
```

### 输入文件格式

#### FASTQ文件命名规范
```
sample001_R1.fastq.gz  # 正向reads
sample001_R2.fastq.gz  # 反向reads
```

#### ASV表格式（TSV）
```
ASV_ID	sample001	sample002	Taxon	Confidence
ASV_1	1234	5678	Bacteroides	0.99
ASV_2	234	456	Prevotella	0.95
```

#### 元数据文件格式（可选）
```
sample_id	age	gender	group
sample001	35	M	healthy
sample002	42	F	patient
```

## 📊 输出结果

### 预处理目录结构
```
preprocessing/
├── merged_asv_taxonomy_table.tsv    # 合并的ASV表（含物种注释）
├── merged_functional_annotation.tsv # 合并的功能注释表(KO)
├── functional_ec_annotation.tsv     # EC酶注释表
├── functional_pathway_annotation.tsv # 通路注释表
├── functional_annotation_stats.json # 功能注释统计
├── functional_summary.json          # 功能摘要
├── processing_stats.json            # 处理统计
├── sample_list.txt                  # 样本列表
├── functional_prediction/           # PICRUSt2功能预测结果
│   ├── KO_metagenome_out/          # KO预测
│   ├── EC_metagenome_out/          # EC预测
│   └── pathways_out/               # 通路预测
├── single_samples/                  # 各样本处理中间文件
└── logs/                            # 处理日志
```

### 完整输出目录
```
output/
├── preprocessing/              # 预处理与功能注释结果（见上）
├── analysis_results/          # 分析结果
│   ├── sample001/
│   │   ├── diversity/        # 多样性分析
│   │   ├── enterotype/       # 肠型判定
│   │   ├── bacteria_scores/  # 菌群评分
│   │   ├── disease_risk/     # 疾病风险
│   │   ├── age_prediction/   # 年龄预测
│   │   ├── functional_prediction/  # 功能预测
│   │   └── cn_annotations.json     # 中文注释
│   └── sample002/
├── reports/                   # HTML报告
│   ├── sample001_report.html
│   └── sample002_report.html
└── logs/                      # 运行日志
```

### 关键输出文件

1. **HTML报告** (`reports/sample_id_report.html`)
   - 综合健康评分和等级
   - 多样性指标可视化
   - 菌群组成饼图和柱状图
   - 疾病风险雷达图
   - 功能预测结果
   - 个性化健康建议

2. **JSON结果** (每个分析模块一个)
   - `basic_analysis.json`: 多样性和B/F比值
   - `enterotype_analysis.json`: 肠型分类结果
   - `bacteria_evaluation.json`: 菌群健康评分
   - `disease_risk_assessment.json`: 14种疾病风险
   - `age_prediction.json`: 生物年龄预测
   - `functional_prediction.json`: 功能预测结果
   - `cn_annotations.json`: 中文注释汇总

## 🔬 分析模块详解

### 1. 多样性分析
- **Shannon指数**: 评估物种多样性和均匀度
- **Simpson指数**: 物种优势度
- **Chao1指数**: 物种丰富度估计
- **B/F比值**: 拟杆菌门/厚壁菌门比值（正常范围0.84-4.94）

### 2. 肠型分类
基于属水平组成的K-means聚类（K=3）：
- **拟杆菌型** (Bacteroides-dominant): 高蛋白饮食相关
- **普氏菌型** (Prevotella-dominant): 高纤维饮食相关
- **瘤胃球菌型** (Ruminococcus-dominant): 混合型

### 3. 菌群健康评估
评估43种核心细菌：
- **13种有益菌**: 双歧杆菌、乳酸杆菌、粪杆菌等
- **15种条件致病菌**: 韦荣氏菌、克雷伯菌、肠球菌等
- **15种有害菌**: 志贺氏菌、沙门氏菌、艰难梭菌等

### 4. 疾病风险评估
基于gutMDisorder数据库，覆盖14种常见疾病：
- **消化系统**: IBD、IBS、便秘、息肉、结直肠癌
- **代谢疾病**: 糖尿病、代谢综合征、肝病
- **心血管**: 冠心病、高血压
- **神经精神**: 抑郁、阿尔茨海默病
- **其他**: 痛风、湿疹

### 5. 功能预测 (PICRUSt2)
- **维生素合成**: 12种维生素（B1-B12、K等）的合成能力
- **短链脂肪酸**: 丁酸、丙酸、乙酸产生能力
- **代谢通路**: KEGG通路富集分析
- **酶功能**: EC酶分类与功能预测

## ⚙️ 高级配置

### 自定义参数
```bash
# 修改并行任务数
./pipeline.sh input output -j 16  # 使用16个并行任务

# 修改线程数
./pipeline.sh input output -t 8   # 每个任务使用8线程

# 跳过特定步骤
./pipeline.sh input output --skip-taxonomy     # 跳过物种注释
./pipeline.sh input output --skip-functional   # 跳过功能预测
```

### 数据库说明

1. **SILVA数据库** (`silva_138_99_16S_338f_806r_classifier.qza`)
   - 用于16S rRNA基因物种注释
   - 已针对338F-806R引物优化

2. **疾病关联数据库** (`disease_associations.json`)
   - 基于gutMDisorder数据库处理
   - 包含14种疾病的菌群关联信息
   - 原始数据和处理脚本在`database/process_bak/`

3. **中文注释数据库**
   - `core_bacteria_annotations.json`: 核心菌群中文名称和描述
   - `core_ec_translations.json`: EC酶中文翻译
   - `core_pathway_translations.json`: KEGG通路中文翻译

## 🐛 故障排除

### 常见问题

1. **内存不足错误**
   - 减少并行任务数: `-j 2`
   - 分批处理样本

2. **QIIME2未找到**
   - 确保已激活conda环境: `conda activate qiime2-amplicon-2023.9`
   - 检查环境变量设置

3. **PICRUSt2功能预测失败**
   - 确保已安装PICRUSt2: `conda install -c bioconda picrust2`
   - 检查ASV表格式是否正确
   - 确保有代表序列文件

4. **物种注释失败**
   - 检查SILVA数据库文件是否存在
   - 确保文件格式为.qza
   - 验证数据库文件完整性

5. **报告生成失败**
   - 检查所有分析模块是否完成
   - 查看日志文件获取详细错误信息
   - 确保中文字体支持

### 获取帮助
```bash
# 查看主脚本帮助
./pipeline.sh --help

# 查看预处理脚本帮助
bash scripts/preprocessing/backend_process.sh --help

# 查看各模块帮助
python scripts/analysis/1_basic_analysis.py --help
```

## 📚 参考文献

1. Turnbaugh PJ, et al. A core gut microbiome in obese and lean twins. Nature. 2009
2. Qin J, et al. A human gut microbial gene catalogue established by metagenomic sequencing. Nature. 2010
3. Human Microbiome Project Consortium. Structure, function and diversity of the healthy human microbiome. Nature. 2012
4. Arumugam M, et al. Enterotypes of the human gut microbiome. Nature. 2011
5. Galkin F, et al. Human Gut Microbiome Aging Clock Based on Taxonomic Profiling and Deep Learning. iScience. 2020
6. Douglas GM, et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol. 2020

## 📄 许可证

本项目采用MIT许可证。详见[LICENSE](LICENSE)文件。

## 👥 贡献

欢迎提交Issue和Pull Request！

项目地址：[https://github.com/chenziqing0111/16s_bighealth_pipe](https://github.com/chenziqing0111/16s_bighealth_pipe)

## 📧 联系方式

如有问题或建议，请通过以下方式联系：
- GitHub Issues: [项目Issues页面](https://github.com/chenziqing0111/16s_bighealth_pipe/issues)

## 🙏 致谢

- QIIME2开发团队
- PICRUSt2开发团队
- SILVA数据库维护者
- gutMDisorder数据库作者
- 所有开源贡献者

---
**最后更新**: 2024年12月
