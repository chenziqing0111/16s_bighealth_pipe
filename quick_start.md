1.环境
ssh -p 7002 user010@43.138.80.104

bh(快捷键，直达工作目录：/data/user010/lvxz_dev/16s_bighealth_pipe)

conda activate /home/user010/miniconda3/envs/qiime2-2024.2



2.测试
# 全流程
bash pipeline.sh ./data ./backend_output_1 -d ./database/silva_138_99_16S_338f_806r_classifier.qza -j 6 -t 6

# 已有预处理数据，仅进行下游分析及生成报告
bash pipeline.sh backend_output/preprocessing backend_output --skip-preprocessing -j 1

3.目录直达
# 项目结构
主脚本：/data/user010/lvxz_dev/16s_bighealth_pipe/pipeline.sh
模块脚本目录：/data/user010/lvxz_dev/16s_bighealth_pipe/scripts
数据库目录：/data/user010/lvxz_dev/16s_bighealth_pipe/database

# 测试分析
数据目录：/data/user010/lvxz_dev/16s_bighealth_pipe/data
预处理目录：/data/user010/lvxz_dev/16s_bighealth_pipe/backend_output/preprocessing
下游分析目录：/data/user010/lvxz_dev/16s_bighealth_pipe/backend_output/analysis_results
报告目录：/data/user010/lvxz_dev/16s_bighealth_pipe/backend_output/reports
