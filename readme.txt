gut_microbiome_pipeline/
│
├── 📂 scripts/
│   ├── 📂 preprocessing/
│   │   ├── backend_process.sh       # FASTQ→ASV表（多样本）
│   │   └── merge_tables.py          
│   │
│   ├── 📂 analysis/                  # 都从ASV表开始
│   │   ├── 1_basic_analysis.py      
│   │   ├── 2_enterotype.py          
│   │   ├── 3_bacteria_eval.py       
│   │   ├── 4_disease_risk.py        
│   │   └── 5_age_predict.py         
│   │
│   └── 📂 report/
│       ├── generate_report.py       # 生成报告
│       └── 📂 templates/            
│           ├── report_template.html
│           ├── report_styles.css
│           └── report_scripts.js
│
├── 📂 database/
│   └── ...（数据库文件）
│
├── 📄 pipeline.sh                    # 智能处理脚本（单个/批量自动识别）
├── 📄 requirements.txt
└── 📄 README.md


output/
├── preprocessing/              # FASTQ预处理结果
│   ├── merged_asv_taxonomy_table.tsv
│   └── sample_list.txt
├── analysis_results/          # 每个样本的分析结果
│   ├── sample001/
│   │   ├── diversity/
│   │   ├── enterotype/
│   │   ├── bacteria_scores/
│   │   ├── disease_risk/
│   │   └── age_prediction/
│   └── sample002/
├── reports/                   # HTML报告
│   ├── sample001_report.html
│   └── sample002_report.html
└── logs/                      # 日志文件
    ├── pipeline_20241220_143022.log
    ├── sample001.log
    └── sample002.log


🎯 完整分析流程（基于现成资源）
📊 Phase 1: 基础分析（你已有代码支持）
输入：ASV表 + 物种注释表 + 样本信息
↓
1. Alpha多样性
   - Shannon、Simpson、Chao1指数
   - 多样性评分（正常/偏低/偏高）

2. Beta多样性  
   - PCoA/PCA分析
   - 聚类分析

3. 物种组成
   - 门/纲/目/科/属/种各水平相对丰度
   - Top10优势菌展示
   
4. B/F比值计算
   - 拟杆菌门/厚壁菌门比值
   - 评估：正常范围0.84-4.94
🔬 Phase 2: 肠型分析
基于属水平组成
↓
K-means聚类（K=3）
↓
判断肠型：
- 拟杆菌型（Bacteroides主导）
- 普氏菌型（Prevotella主导）  
- 瘤胃球菌型（Ruminococcus主导）
🦠 Phase 3: 特定菌群评估
从ASV表提取目标菌的相对丰度
↓
分三类统计：

1. 核心有益菌（13个属/种）
   - 双歧杆菌、乳酸杆菌、粪杆菌等
   - 计算：总体有益菌评分

2. 条件致病菌（15个属/种）
   - 韦荣氏菌、克雷伯菌、肠球菌等
   - 计算：条件致病菌负荷

3. 核心有害菌（15个属/种）
   - 志贺氏菌、沙门氏菌、艰难梭菌等
   - 计算：有害菌风险评分
🧬 Phase 4: 功能预测（PICRUSt2）
输入：ASV表 + 代表序列
↓
PICRUSt2分析
↓
输出：
1. 维生素合成能力（12种维生素）
2. 氨基酸代谢能力（20种氨基酸）
3. 短链脂肪酸产生能力
4. 天然色素合成（类胡萝卜素、花青素等）
🏥 Phase 5: 疾病风险评估（基于数据库）
使用gutMDisorder数据库：
下载数据库CSV文件
↓
对每种疾病：
1. 提取该疾病相关的所有菌
2. 计算公式：
   风险分 = Σ(上调菌×丰度×权重) - Σ(下调菌×丰度×权重)
↓
标准化到0-100分
↓
分级：
- 0-40: 低风险
- 40-70: 中风险  
- 70-100: 高风险
覆盖的疾病（按你报告的14种）：

消化系统：IBD、IBS、便秘、息肉、结直肠癌
代谢疾病：糖尿病、代谢综合征、肝病
心血管：冠心病、高血压
神经精神：抑郁、阿尔茨海默
其他：痛风、湿疹

🎂 Phase 6: 生物年龄预测
使用Galkin模型（2020）：
GitHub克隆代码
↓
输入：39个核心属的相对丰度
↓
深度学习模型预测
↓
输出：生物年龄（岁）
↓
计算：肠道年龄差异 = 生物年龄 - 实际年龄
🍎 Phase 7: 个性化建议生成
基于以上所有分析结果
↓
1. 饮食建议
   - 根据B/F比值 → 调整膳食纤维/蛋白质
   - 根据有益菌水平 → 推荐益生元食物
   - 根据疾病风险 → 个性化忌口建议

2. 益生菌补充建议  
   - 缺乏的有益菌 → 对应益生菌产品
   - 标准：双歧杆菌<0.183% 建议补充

3. 生活方式建议
   - 基于整体评分和疾病风险
📈 Phase 8: 报告生成
汇总所有分析结果
↓
生成HTML报告，包含：
1. 菌群得分总览（雷达图）
2. 详细检测数据（表格）
3. 风险评估（柱状图）
4. 个性化建议（文字）

✅ 所需资源清单
数据库：

✅ SILVA 138（物种注释）- 你已有
📥 gutMDisorder（疾病关联）- 需下载
📥 MicrobiomeAgeClock（年龄模型）- GitHub

软件工具：

✅ QIIME2（基础分析）- 你已有
📥 PICRUSt2（功能预测）- 需安装
✅ Python（数据处理）- 你已有

参考范围：

✅ 各菌正常范围（从文献/数据库获取）
✅ 疾病风险阈值（gutMDisorder提供）
✅ B/F比值范围（0.84-4.94）


