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
