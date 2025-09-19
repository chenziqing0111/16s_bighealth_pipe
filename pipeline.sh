#!/bin/bash
#=============================================================================
# 肠道菌群分析流程 - 智能主控脚本
# 功能：
# 1. 自动识别输入类型（FASTQ目录 或 ASV表）
# 2. 自动识别样本数（单个/批量）
# 3. 并行处理多个样本
# 4. 生成完整分析报告
#=============================================================================

set -e  # 出错即停止

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 默认参数
THREADS=4
PARALLEL_JOBS=4

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 输出函数
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# 使用说明
usage() {
    cat << EOF
肠道菌群分析流程 - 使用说明

用法: 
  $0 INPUT OUTPUT [OPTIONS]

参数:
  INPUT               输入路径（FASTQ目录 或 ASV表文件）
  OUTPUT              输出目录

可选参数:
  -m, --metadata      元数据文件（包含样本信息）
  -d, --database      SILVA数据库文件（.qza格式）
  -t, --threads       每个任务的线程数（默认: $THREADS）
  -j, --jobs          并行任务数（默认: $PARALLEL_JOBS）
  --skip-preprocessing 跳过预处理（输入必须是ASV表）
  --skip-taxonomy     跳过物种注释
  -h, --help          显示帮助

示例:
  # 从FASTQ开始完整分析
  $0 /path/to/fastq_dir /path/to/output -d silva.qza
  
  # 从ASV表开始分析
  $0 merged_asv_table.tsv /path/to/output --skip-preprocessing
  
  # 带元数据的批量分析
  $0 /path/to/fastq_dir /path/to/output -m metadata.tsv -j 8

输出:
  OUTPUT/
  ├── preprocessing/          # 预处理结果
  ├── analysis_results/       # 分析结果
  │   └── sample_id/         # 每个样本的结果
  ├── reports/               # HTML报告
  └── logs/                  # 日志文件

EOF
    exit 0
}

# 检查依赖
check_dependencies() {
    log "检查环境依赖..."
    
    # 检查Python
    if ! command -v python3 &> /dev/null; then
        error "未找到Python3"
        exit 1
    fi
    
    # 检查QIIME2（如果需要预处理）
    if [ "$SKIP_PREPROCESSING" != "true" ]; then
        if ! command -v qiime &> /dev/null; then
            warn "未找到QIIME2，将跳过预处理步骤"
            SKIP_PREPROCESSING="true"
        fi
    fi
    
    # 检查必需的Python包
    python3 -c "import pandas, numpy, json" 2>/dev/null || {
        error "缺少必需的Python包，请运行: pip install pandas numpy"
        exit 1
    }
    
    # 检查可选的Python包
    python3 -c "import sklearn" 2>/dev/null || {
        warn "未安装scikit-learn，部分功能可能受限"
    }
    
    info "✓ 环境检查完成"
}

# 预处理FASTQ文件
run_preprocessing() {
    local input_dir="$1"
    local output_dir="$2"
    
    log "Step 1: 预处理FASTQ文件"
    
    # 检查预处理脚本
    if [ ! -f "$SCRIPT_DIR/scripts/preprocessing/backend_process.sh" ]; then
        error "找不到预处理脚本: backend_process.sh"
        exit 1
    fi
    
    # 构建预处理命令
    local cmd="bash $SCRIPT_DIR/scripts/preprocessing/backend_process.sh"
    cmd="$cmd --from-fastq $input_dir"
    cmd="$cmd --output $output_dir/preprocessing"
    cmd="$cmd --threads $THREADS"
    cmd="$cmd --parallel $PARALLEL_JOBS"
    
    # 添加数据库参数
    if [ -n "$SILVA_DB" ] && [ -f "$SILVA_DB" ]; then
        cmd="$cmd --database $SILVA_DB"
    elif [ "$SKIP_TAXONOMY" != "true" ]; then
        warn "未指定SILVA数据库，将跳过物种注释"
    fi
    
    # 执行预处理
    info "执行命令: $cmd"
    eval $cmd 2>&1 | tee -a "$LOG_FILE"
    
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        error "预处理失败"
        exit 1
    fi
    
    # 检查输出文件
    ASV_TABLE="$output_dir/preprocessing/merged_asv_taxonomy_table.tsv"
    if [ ! -f "$ASV_TABLE" ]; then
        error "预处理未生成ASV表: $ASV_TABLE"
        exit 1
    fi
    
    info "✓ 预处理完成，生成ASV表: $ASV_TABLE"
}

# 处理单个样本的分析
process_sample() {
    local sample_id="$1"
    local asv_table="$2"
    local output_base="$3"
    local metadata_file="$4"
    
    log "处理样本: $sample_id"
    
    # 创建样本输出目录
    local sample_dir="$output_base/analysis_results/$sample_id"
    mkdir -p "$sample_dir"/{diversity,enterotype,bacteria_scores,disease_risk,age_prediction}
    
    # 创建样本日志
    local sample_log="$output_base/logs/${sample_id}.log"
    
    # 提取单样本ASV数据
    log "  提取样本数据..."
    python3 << PYTHON 2>> "$sample_log"
import pandas as pd
import sys

try:
    # 读取ASV表
    df = pd.read_csv('$asv_table', sep='\t', index_col=0)
    
    # 识别分类列
    tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    tax_cols = [col for col in tax_cols if col in df.columns]
    
    # 检查样本是否存在
    if '$sample_id' not in df.columns:
        print(f"错误: 样本 $sample_id 不在ASV表中", file=sys.stderr)
        sys.exit(1)
    
    # 提取样本数据和分类信息
    sample_data = df[['$sample_id'] + tax_cols]
    
    # 过滤零值行（可选，保留所有数据更准确）
    # sample_data = sample_data[sample_data['$sample_id'] > 0]
    
    # 保存
    sample_data.to_csv('$sample_dir/sample_asv.tsv', sep='\t')
    print(f"  ✓ 提取 {len(sample_data)} 个ASVs")
    
except Exception as e:
    print(f"错误: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON

    if [ $? -ne 0 ]; then
        error "  样本数据提取失败"
        return 1
    fi
    
    # 并行运行所有分析
    log "  运行分析模块..."
    
    {
        # 1. 基础分析
        python3 "$SCRIPT_DIR/scripts/analysis/1_basic_analysis.py" \
            --input "$sample_dir/sample_asv.tsv" \
            --output "$sample_dir/diversity" \
            >> "$sample_log" 2>&1 &
        
        # 2. 肠型分析
        python3 "$SCRIPT_DIR/scripts/analysis/2_enterotype.py" \
            --input "$sample_dir/sample_asv.tsv" \
            --output "$sample_dir/enterotype" \
            >> "$sample_log" 2>&1 &
        
        # 3. 菌群评估
        if [ -f "$SCRIPT_DIR/database/normal_ranges.json" ]; then
            python3 "$SCRIPT_DIR/scripts/analysis/3_bacteria_eval.py" \
                --input "$sample_dir/sample_asv.tsv" \
                --ranges "$SCRIPT_DIR/database/normal_ranges.json" \
                --output "$sample_dir/bacteria_scores" \
                >> "$sample_log" 2>&1 &
        else
            python3 "$SCRIPT_DIR/scripts/analysis/3_bacteria_eval.py" \
                --input "$sample_dir/sample_asv.tsv" \
                --output "$sample_dir/bacteria_scores" \
                >> "$sample_log" 2>&1 &
        fi
        
        # 4. 疾病风险
        if [ -f "$SCRIPT_DIR/database/gutMDisorder.csv" ]; then
            python3 "$SCRIPT_DIR/scripts/analysis/4_disease_risk.py" \
                --input "$sample_dir/sample_asv.tsv" \
                --database "$SCRIPT_DIR/database/gutMDisorder.csv" \
                --output "$sample_dir/disease_risk" \
                >> "$sample_log" 2>&1 &
        else
            python3 "$SCRIPT_DIR/scripts/analysis/4_disease_risk.py" \
                --input "$sample_dir/sample_asv.tsv" \
                --output "$sample_dir/disease_risk" \
                >> "$sample_log" 2>&1 &
        fi
        
        # 5. 年龄预测
        # 如果有元数据，尝试获取年龄
        local age_param=""
        if [ -n "$metadata_file" ] && [ -f "$metadata_file" ]; then
            age=$(python3 -c "
import pandas as pd
try:
    df = pd.read_csv('$metadata_file', sep='\t')
    # 尝试不同的列名
    for col in ['age', 'Age', 'AGE', '年龄']:
        if col in df.columns:
            row = df[df.iloc[:,0] == '$sample_id']
            if not row.empty:
                print(int(row[col].iloc[0]))
                break
except:
    pass
" 2>/dev/null)
            
            if [ -n "$age" ]; then
                age_param="--age $age"
            fi
        fi
        
        if [ -f "$SCRIPT_DIR/database/age_markers.json" ]; then
            python3 "$SCRIPT_DIR/scripts/analysis/5_age_predict.py" \
                --input "$sample_dir/sample_asv.tsv" \
                --markers "$SCRIPT_DIR/database/age_markers.json" \
                $age_param \
                --output "$sample_dir/age_prediction" \
                >> "$sample_log" 2>&1 &
        else
            python3 "$SCRIPT_DIR/scripts/analysis/5_age_predict.py" \
                --input "$sample_dir/sample_asv.tsv" \
                $age_param \
                --output "$sample_dir/age_prediction" \
                >> "$sample_log" 2>&1 &
        fi
        
        # 6. 功能预测
        python3 "$SCRIPT_DIR/scripts/analysis/6_functional_prediction.py" \
            --input "$sample_dir/sample_asv.tsv" \
            --output "$sample_dir/functional_prediction" \
            >> "$sample_log" 2>&1 &
                
        # 等待所有分析完成
        wait
    }
       
    # 7. 中文注释（新增）
    log "  运行中文注释..."
    
    # 检查注释脚本是否存在
    if [ -f "$SCRIPT_DIR/scripts/analysis/7_cn_annotation.py" ]; then
        # 检查中文注释数据库
        local database_dir="${SCRIPT_DIR}/database"
        
        if [ -f "$database_dir/core_bacteria_annotations.json" ] && \
           [ -f "$database_dir/core_pathway_translations.json" ] && \
           [ -f "$database_dir/core_ec_translations.json" ]; then
            
            # 运行中文注释
            python3 "$SCRIPT_DIR/scripts/analysis/7_cn_annotation.py" \
                --sample-dir "$sample_dir" \
                --database "$database_dir" \
                >> "$sample_log" 2>&1
            
            if [ $? -eq 0 ]; then
                log "  ✓ 中文注释完成"
            else
                warn "  中文注释执行出错，将使用英文信息"
            fi
        else
            warn "  未找到完整的中文注释数据库"
        fi
    fi
    
    # 检查分析结果
    local success=true
    for module in diversity enterotype bacteria_scores disease_risk age_prediction functional_prediction; do
        if [ ! -f "$sample_dir/$module/"*.json ]; then
            warn "  模块 $module 未生成结果"
            success=false
        fi
    done
    
    if [ "$success" = true ]; then
        log "  ✓ 所有分析完成"
        
        # 生成报告
        log "  生成HTML报告..."
        
        if [ -f "$SCRIPT_DIR/scripts/report/generate_report.py" ]; then
            local report_cmd="python3 $SCRIPT_DIR/scripts/report/generate_report.py"
            report_cmd="$report_cmd --sample-dir $sample_dir"  # 改这里：使用 --sample-dir
            report_cmd="$report_cmd --output $output_base/reports/${sample_id}_report.html"
            
            # 注意：可能不需要 --metadata 参数了，先注释掉
            # if [ -n "$metadata_file" ] && [ -f "$metadata_file" ]; then
            #     report_cmd="$report_cmd --metadata $metadata_file"
            # fi
            
            eval $report_cmd >> "$sample_log" 2>&1
            
            if [ $? -eq 0 ]; then
                info "  ✓ 报告生成: $output_base/reports/${sample_id}_report.html"
            else
                warn "  报告生成失败，请查看日志: $sample_log"
            fi
        else
            warn "  未找到报告生成脚本"
        fi
    else
        error "  部分分析失败，请查看日志: $sample_log"
    fi
    
    log "  样本 $sample_id 处理完成"
}

# 主函数
main() {
    # 解析参数
    INPUT=""
    OUTPUT_DIR=""
    METADATA=""
    SILVA_DB=""
    SKIP_PREPROCESSING=""
    SKIP_TAXONOMY=""
    
    # 至少需要2个参数
    if [ $# -lt 2 ]; then
        usage
    fi
    
    # 获取必需参数
    INPUT="$1"
    OUTPUT_DIR="$2"
    shift 2
    
    # 解析可选参数
    while [ $# -gt 0 ]; do
        case "$1" in
            -m|--metadata)
                METADATA="$2"
                shift 2
                ;;
            -d|--database)
                SILVA_DB="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -j|--jobs)
                PARALLEL_JOBS="$2"
                shift 2
                ;;
            --skip-preprocessing)
                SKIP_PREPROCESSING="true"
                shift
                ;;
            --skip-taxonomy)
                SKIP_TAXONOMY="true"
                shift
                ;;
            -h|--help)
                usage
                ;;
            *)
                warn "未知参数: $1"
                shift
                ;;
        esac
    done
    
    # 验证输入
    if [ ! -e "$INPUT" ]; then
        error "输入路径不存在: $INPUT"
        exit 1
    fi
    
    # 创建输出目录结构
    mkdir -p "$OUTPUT_DIR"/{preprocessing,analysis_results,reports,logs}
    
    # 设置日志文件
    LOG_FILE="$OUTPUT_DIR/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
    exec 1> >(tee -a "$LOG_FILE")
    exec 2>&1
    
    # 打印配置
    echo "========================================="
    echo "肠道菌群分析流程"
    echo "========================================="
    echo "开始时间: $(date)"
    echo "输入路径: $INPUT"
    echo "输出目录: $OUTPUT_DIR"
    echo "线程数: $THREADS"
    echo "并行任务: $PARALLEL_JOBS"
    [ -n "$METADATA" ] && echo "元数据: $METADATA"
    [ -n "$SILVA_DB" ] && echo "数据库: $SILVA_DB"
    echo "========================================="
    
    # 检查依赖
    check_dependencies
    
    # 判断输入类型并处理
    # 判断输入类型并处理
    if [ -d "$INPUT" ]; then
        # 输入是目录
        
        # 检查是否是预处理输出目录（包含merged_asv_taxonomy_table.tsv）
        if [ -f "$INPUT/merged_asv_taxonomy_table.tsv" ]; then
            log "检测到预处理输出目录，使用其中的ASV表..."
            ASV_TABLE="$INPUT/merged_asv_taxonomy_table.tsv"
            PREPROCESSING_DIR="$INPUT"  # 保存预处理目录路径
            info "✓ 找到ASV表: $ASV_TABLE"
            
            # 导出预处理目录供后续使用
            export PREPROCESSING_DIR
            
        elif [ -f "$INPUT/preprocessing/merged_asv_taxonomy_table.tsv" ]; then
            # 可能是完整的输出目录
            log "检测到完整输出目录，使用预处理结果..."
            ASV_TABLE="$INPUT/preprocessing/merged_asv_taxonomy_table.tsv"
            PREPROCESSING_DIR="$INPUT/preprocessing"
            info "✓ 找到ASV表: $ASV_TABLE"
            
            # 导出预处理目录供后续使用
            export PREPROCESSING_DIR
            
        elif [ "$SKIP_PREPROCESSING" = "true" ]; then
            error "目录中未找到ASV表，且指定了--skip-preprocessing"
            exit 1
            
        else
            # 是FASTQ目录，需要预处理
            log "检测到FASTQ目录，开始预处理..."
            run_preprocessing "$INPUT" "$OUTPUT_DIR"
            PREPROCESSING_DIR="$OUTPUT_DIR/preprocessing"
            export PREPROCESSING_DIR
    fi
        
    elif [ -f "$INPUT" ]; then
        # 输入是文件
        log "检测到文件输入，检查格式..."
        
        # 检查是否是ASV表
        if head -1 "$INPUT" | grep -q -E "ASV|OTU|Feature|Taxon|Confidence"; then
            ASV_TABLE="$INPUT"
            info "✓ 识别为ASV表"
        else
            error "无法识别的文件格式"
            exit 1
        fi
    else
        error "无效的输入: $INPUT"
        exit 1
    fi
    
    # 获取样本列表
    log "Step 2: 识别样本列表"
    
    SAMPLES=$(python3 << PYTHON
import pandas as pd
import sys

try:
    df = pd.read_csv('$ASV_TABLE', sep='\t', index_col=0)
    
    # 排除分类列
    tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in df.columns if col not in tax_cols]
    
    for sample in sample_cols:
        print(sample)
        
except Exception as e:
    print(f"错误: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON
)
    
    if [ -z "$SAMPLES" ]; then
        error "未找到有效样本"
        exit 1
    fi
    
    # 统计样本数
    SAMPLE_COUNT=$(echo "$SAMPLES" | wc -l)
    info "找到 $SAMPLE_COUNT 个样本"
    
    # 显示前5个样本
    echo "$SAMPLES" | head -5 | while read sample; do
        echo "  - $sample"
    done
    [ $SAMPLE_COUNT -gt 5 ] && echo "  ..."
    
    # 处理样本
    log "Step 3: 分析样本"
    
    # 导出函数和变量供并行使用
    export -f process_sample log error warn info
    export SCRIPT_DIR OUTPUT_DIR
    
    if [ $SAMPLE_COUNT -eq 1 ]; then
        # 单样本处理
        info "单样本模式"
        SAMPLE_ID=$(echo "$SAMPLES" | head -1)
        process_sample "$SAMPLE_ID" "$ASV_TABLE" "$OUTPUT_DIR" "$METADATA"
        
    else
        # 多样本并行处理
        info "批量模式: $SAMPLE_COUNT 个样本, $PARALLEL_JOBS 并行"
        
        # 使用GNU parallel或xargs
        if command -v parallel &> /dev/null; then
            echo "$SAMPLES" | parallel -j $PARALLEL_JOBS \
                process_sample {} "$ASV_TABLE" "$OUTPUT_DIR" "$METADATA"
        else
            echo "$SAMPLES" | xargs -P $PARALLEL_JOBS -I {} bash -c \
                'process_sample "$@"' _ {} "$ASV_TABLE" "$OUTPUT_DIR" "$METADATA"
        fi
    fi
    
    # 统计结果
    echo "========================================="
    log "Step 4: 分析完成统计"
    
    # 统计成功的报告
    SUCCESS_COUNT=$(ls -1 "$OUTPUT_DIR/reports/"*.html 2>/dev/null | wc -l)
    FAILED_COUNT=$((SAMPLE_COUNT - SUCCESS_COUNT))
    
    info "成功生成报告: $SUCCESS_COUNT/$SAMPLE_COUNT"
    [ $FAILED_COUNT -gt 0 ] && warn "失败样本: $FAILED_COUNT"
    
    # 列出生成的报告
    if [ $SUCCESS_COUNT -gt 0 ]; then
        echo ""
        echo "生成的报告:"
        ls -1 "$OUTPUT_DIR/reports/"*.html | head -10 | while read report; do
            echo "  - $(basename "$report")"
        done
        [ $SUCCESS_COUNT -gt 10 ] && echo "  ..."
    fi
    
    echo ""
    echo "========================================="
    echo "完成时间: $(date)"
    echo "输出目录: $OUTPUT_DIR"
    echo "日志文件: $LOG_FILE"
    echo "========================================="
    
    # 返回状态
    [ $FAILED_COUNT -eq 0 ] && exit 0 || exit 1
}

# 执行主函数
main "$@"
