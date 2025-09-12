#!/bin/bash
#=============================================================================
# 后端基础分析流程 - 从FASTQ到ASV表（修复版）
# 修复问题：
# 1. SILVA数据库路径处理
# 2. rep_seqs.qza文件保留
# 3. 物种注释确保执行
# 4. 合并表格正确生成
#=============================================================================

# 默认参数
THREADS=4
PARALLEL_JOBS=1
SILVA_DB=""  # 不设置默认值，必须用户提供或从环境变量获取
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 输出函数
log() { 
    echo "[$(date +'%H:%M:%S')] $1"
}

error_exit() {
    echo "[ERROR] $1" >&2
    exit 1
}

warn() {
    echo "[WARN] $1" >&2
}

# 使用说明
usage() {
    cat << EOF
后端基础分析流程 - 从FASTQ生成ASV表

用法: 
  $0 --from-fastq DIR -o OUTPUT -d DATABASE [OPTIONS]

必需参数:
  --from-fastq DIR    FASTQ文件目录
  -o, --output DIR    输出目录
  -d, --database FILE SILVA数据库文件（.qza格式）

可选参数:
  -t, --threads N     每个样本的线程数 (默认: $THREADS)
  -j, --parallel N    并行处理的样本数 (默认: $PARALLEL_JOBS)
  --skip-taxonomy     跳过物种注释
  -h, --help          显示帮助

输出文件:
  OUTPUT/merged_asv_taxonomy_table.tsv  - ASV丰度表+物种注释
  OUTPUT/sample_list.txt                - 样本列表
  OUTPUT/processing_stats.json          - 处理统计

示例:
  # 基础处理
  $0 --from-fastq raw_data -o output -d silva_classifier.qza
  
  # 并行处理
  $0 --from-fastq raw_data -o output -d silva_classifier.qza -j 4

注意：
  SILVA数据库文件必须提供，或设置环境变量 SILVA_DB

EOF
    exit 0
}

# 生成样本列表
generate_sample_list() {
    local fastq_dir="$1"
    local output_dir="$2"
    
    log "扫描FASTQ文件..."
    
    python3 << PYTHON
import os
import re
from pathlib import Path

fastq_dir = Path("$fastq_dir")
output_file = Path("$output_dir") / "sample_list.txt"

def parse_fastq_name(filename):
    """解析FASTQ文件名"""
    patterns = [
        r'(.+?)[._]R([12])[._]',
        r'(.+?)[._]([12])\.',
        r'(.+?)[._](forward|reverse)',
        r'(.+?)[._](fwd|rev)',
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            sample = match.group(1)
            sample = re.sub(r'\.(338F_806R|16S|V[34]-V[34]|raw)', '', sample)
            direction = match.group(2)
            if direction in ['1', 'forward', 'fwd']: 
                return sample, 'R1'
            if direction in ['2', 'reverse', 'rev']: 
                return sample, 'R2'
    return None, None

def infer_group(sample_name):
    """推断分组"""
    group = re.sub(r'[_-]\d+$', '', sample_name)
    return group if group and not group.isdigit() else sample_name

# 查找所有FASTQ文件
fastq_files = list(fastq_dir.glob('**/*.f*q.gz')) + list(fastq_dir.glob('**/*.f*q'))
print(f"  找到 {len(fastq_files)} 个FASTQ文件")

if not fastq_files:
    print("  错误: 没有找到FASTQ文件")
    sys.exit(1)

# 配对文件
samples = {}
for file in fastq_files:
    sample, direction = parse_fastq_name(file.name)
    if sample and direction:
        if sample not in samples:
            samples[sample] = {}
        samples[sample][direction] = str(file.absolute())

# 筛选完整配对
complete = {s: f for s, f in samples.items() if 'R1' in f and 'R2' in f}
print(f"  识别到 {len(complete)} 对完整配对")

if not complete:
    print("  错误: 没有找到完整的配对文件")
    sys.exit(1)

# 写入样本列表
with open(output_file, 'w') as f:
    groups = {}
    for sample in complete:
        group = infer_group(sample)
        groups.setdefault(group, []).append(sample)
    
    for group in sorted(groups):
        for sample in sorted(groups[group]):
            f.write(f"{sample}\t{complete[sample]['R1']}\t{complete[sample]['R2']}\t{group}\n")

print(f"  样本列表保存到: {output_file}")

# 显示分组统计
for group in sorted(groups):
    print(f"    组 {group}: {len(groups[group])} 个样本")
PYTHON
}

# 处理单个样本
process_single_sample() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    local group="$4"
    local outdir="$OUTPUT_DIR/single_samples/$sample"
    
    # 确保临时目录环境变量在子进程中生效
    export TMPDIR="$OUTPUT_DIR/tmp"
    export TEMP="$OUTPUT_DIR/tmp"
    export TMP="$OUTPUT_DIR/tmp"
    export TEMPDIR="$OUTPUT_DIR/tmp"
    
    # Python也需要这些
    export PYTHONTMPDIR="$OUTPUT_DIR/tmp"
    
    mkdir -p "$outdir/exported"
    
    log "处理样本: $sample (组: $group)"
    
    # Step 1: 导入数据
    echo -e "sample-id,absolute-filepath,direction\n$sample,$r1,forward\n$sample,$r2,reverse" > "$outdir/manifest.csv"
    
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$outdir/manifest.csv" \
        --output-path "$outdir/data.qza" \
        --input-format PairedEndFastqManifestPhred33 2>&1 | \
        grep -E "Imported|saved" || true
    
    # Step 2: 去除引物
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$outdir/data.qza" \
        --p-front-f ACTCCTACGGGAGGCAGCA \
        --p-front-r GGACTACHVGGGTWTCTAAT \
        --o-trimmed-sequences "$outdir/trimmed.qza" \
        --quiet
    
    # Step 3: DADA2降噪
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$outdir/trimmed.qza" \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-max-ee-f 5 \
        --p-max-ee-r 5 \
        --p-n-threads $THREADS \
        --o-table "$outdir/table.qza" \
        --o-representative-sequences "$outdir/rep_seqs.qza" \
        --o-denoising-stats "$outdir/stats.qza" \
        --quiet
    
    # Step 4: 导出基础结果
    for item in table rep_seqs stats; do
        qiime tools export \
            --input-path "$outdir/${item}.qza" \
            --output-path "$outdir/exported" 2>&1 | \
            grep -E "Exported|saved" || true
    done
    
    # 转换BIOM为TSV
    biom convert \
        -i "$outdir/exported/feature-table.biom" \
        -o "$outdir/asv_table.tsv" \
        --to-tsv
    
    # Step 5: 物种分类（检查全局和环境变量中的SILVA_DB）
    local silva_path="${SILVA_DB}"
    
    # 如果是相对路径，转换为绝对路径
    if [ -n "$silva_path" ] && [[ "$silva_path" != /* ]]; then
        silva_path="$(cd "$(dirname "$silva_path")" && pwd)/$(basename "$silva_path")"
    fi
    
    if [ "$SKIP_TAXONOMY" != "true" ] && [ -n "$silva_path" ] && [ -f "$silva_path" ]; then
        log "  运行物种注释..."
        qiime feature-classifier classify-sklearn \
            --i-classifier "$silva_path" \
            --i-reads "$outdir/rep_seqs.qza" \
            --o-classification "$outdir/taxonomy.qza" \
            --p-n-jobs $THREADS \
            --quiet
        
        qiime tools export \
            --input-path "$outdir/taxonomy.qza" \
            --output-path "$outdir/exported" 2>&1 | \
            grep -E "Exported|saved" || true
            
        log "  物种注释完成"
    else
        if [ "$SKIP_TAXONOMY" = "true" ]; then
            log "  跳过物种注释（--skip-taxonomy）"
        elif [ -z "$silva_path" ]; then
            log "  警告: SILVA数据库路径未设置"
        else
            log "  警告: SILVA数据库不存在: $silva_path"
        fi
    fi
    
    # 清理临时文件（保留rep_seqs.qza用于后续可能的分析）
    rm -f "$outdir"/data.qza "$outdir"/trimmed.qza "$outdir"/manifest.csv 2>/dev/null || true
    # 注意：保留 table.qza, rep_seqs.qza, stats.qza, taxonomy.qza
    
    log "样本 $sample 完成"
}

# 合并所有ASV表和物种注释
merge_asv_tables() {
    log "合并ASV表和物种注释..."
    
    python3 << PYTHON
import pandas as pd
import json
from pathlib import Path

output_dir = Path("$OUTPUT_DIR")
sample_list_file = output_dir / "sample_list.txt"

# 读取样本列表
samples_df = pd.read_csv(sample_list_file, sep='\t', names=['sample_id', 'r1', 'r2', 'group'])

# 收集所有ASV表
all_asvs = set()
asv_data = {}
taxonomy_data = None  # 只需要读取一次
stats_data = []

log_msg = []

for _, row in samples_df.iterrows():
    sample_id = row['sample_id']
    sample_dir = output_dir / 'single_samples' / sample_id
    
    # 读取ASV表
    asv_file = sample_dir / 'asv_table.tsv'
    if asv_file.exists():
        try:
            df = pd.read_csv(asv_file, sep='\t', comment='#', index_col=0)
            # 跳过可能的标题行
            if len(df) > 0 and (df.index[0] == 'OTU ID' or str(df.index[0]).startswith('#')):
                df = pd.read_csv(asv_file, sep='\t', skiprows=1, index_col=0)
            
            all_asvs.update(df.index)
            asv_data[sample_id] = df[sample_id] if sample_id in df.columns else df.iloc[:, 0]
        except Exception as e:
            log_msg.append(f"  警告: 无法读取 {sample_id} 的ASV表: {e}")
    
    # 读取物种注释（只需要从一个样本读取）
    if taxonomy_data is None:
        tax_file = sample_dir / 'exported' / 'taxonomy.tsv'
        if tax_file.exists():
            try:
                tax_df = pd.read_csv(tax_file, sep='\t')
                if 'Feature ID' in tax_df.columns:
                    tax_df.set_index('Feature ID', inplace=True)
                elif tax_df.columns[0] not in ['Taxon', 'Confidence']:
                    tax_df.set_index(tax_df.columns[0], inplace=True)
                taxonomy_data = tax_df
                log_msg.append(f"  从 {sample_id} 读取物种注释: {len(taxonomy_data)} 个ASVs")
            except Exception as e:
                log_msg.append(f"  警告: 无法读取物种注释: {e}")
    
    # 读取统计信息
    stats_file = sample_dir / 'exported' / 'stats.tsv'
    if stats_file.exists():
        try:
            with open(stats_file, 'r') as f:
                lines = [line for line in f.readlines() if not line.startswith('#')]
            
            if len(lines) > 1:
                import io
                stats_df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
                if len(stats_df) > 0:
                    stats_row = stats_df.iloc[0].to_dict()
                    stats_row['sample_id'] = sample_id
                    stats_row['group'] = row['group']
                    stats_data.append(stats_row)
        except:
            pass

# 打印日志
for msg in log_msg:
    print(msg)

# 创建合并的ASV表
print(f"  合并 {len(asv_data)} 个样本的ASV数据...")
merged_asv = pd.DataFrame(index=sorted(all_asvs))
for sample_id, data in asv_data.items():
    merged_asv[sample_id] = data

merged_asv = merged_asv.fillna(0).astype(int)

# 准备最终的合并表
merged_table = merged_asv.copy()

# 如果有物种注释，添加到表格
if taxonomy_data is not None and len(taxonomy_data) > 0:
    print(f"  添加物种注释...")
    
    # 确保索引是字符串类型
    merged_table.index = merged_table.index.astype(str)
    taxonomy_data.index = taxonomy_data.index.astype(str)
    
    # 找到共同的ASVs
    common_asvs = list(set(merged_table.index) & set(taxonomy_data.index))
    print(f"    共同ASVs: {len(common_asvs)} / {len(merged_table)}")
    
    # 添加注释列
    if 'Taxon' in taxonomy_data.columns:
        merged_table['Taxon'] = taxonomy_data['Taxon']
        merged_table['Taxon'] = merged_table['Taxon'].fillna('Unclassified')
    
    if 'Confidence' in taxonomy_data.columns:
        merged_table['Confidence'] = taxonomy_data['Confidence']
        merged_table['Confidence'] = merged_table['Confidence'].fillna(0.0)
    
    has_taxonomy = True
else:
    # 没有物种注释，添加空列
    print("  警告: 没有找到物种注释，添加空注释列")
    merged_table['Taxon'] = 'Unclassified'
    merged_table['Confidence'] = 0.0
    has_taxonomy = False

# 保存合并的表格（总是保存为同一个文件名）
output_file = output_dir / 'merged_asv_taxonomy_table.tsv'
merged_table.to_csv(output_file, sep='\t')
print(f"  合并表保存到: {output_file}")
print(f"    - 总ASV数: {len(merged_asv)}")
print(f"    - 总样本数: {len(asv_data)}")
print(f"    - 总reads数: {merged_asv.sum().sum():,}")
print(f"    - 包含物种注释: {'是' if has_taxonomy else '否'}")

# 保存处理统计
stats_summary = {
    'total_samples': len(samples_df),
    'total_asvs': len(merged_asv),
    'total_reads': int(merged_asv.sum().sum()),
    'samples': stats_data,
    'groups': samples_df['group'].value_counts().to_dict(),
    'has_taxonomy': has_taxonomy
}

output_stats = output_dir / 'processing_stats.json'
with open(output_stats, 'w') as f:
    json.dump(stats_summary, f, indent=2)
print(f"  处理统计保存到: {output_stats}")
PYTHON
}

# 并行处理
run_parallel_processing() {
    local sample_list="$1"
    local max_jobs="$2"
    local job_count=0
    
    log "开始并行处理（最大 $max_jobs 个任务）"
    
    while IFS=$'\t' read -r sample r1 r2 group; do
        [ -z "$sample" ] && continue
        
        # 等待空闲槽位
        while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
            sleep 2
        done
        
        # 启动后台任务
        {
            process_single_sample "$sample" "$r1" "$r2" "$group"
        } > "$OUTPUT_DIR/logs/${sample}.log" 2>&1 &
        
        job_count=$((job_count + 1))
        log "  启动任务 $job_count: $sample (PID: $!)"
        
    done < "$sample_list"
    
    # 等待所有任务完成
    log "等待所有任务完成..."
    wait
    
    log "所有 $job_count 个任务已完成"
}

# 主函数
main() {
    # 初始化变量
    FASTQ_DIR=""
    OUTPUT_DIR=""
    SKIP_TAXONOMY=""
    
    # 解析参数
    while [ $# -gt 0 ]; do
        case "$1" in
            --from-fastq)
                FASTQ_DIR="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -j|--parallel)
                PARALLEL_JOBS="$2"
                shift 2
                ;;
            -d|--database)
                SILVA_DB="$2"
                shift 2
                ;;
            --skip-taxonomy)
                SKIP_TAXONOMY="true"
                shift
                ;;
            -h|--help)
                usage
                ;;
            *)
                shift
                ;;
        esac
    done
    
    # 验证参数
    if [ -z "$FASTQ_DIR" ]; then
        error_exit "请指定FASTQ目录 (--from-fastq)"
    fi
    
    if [ -z "$OUTPUT_DIR" ]; then
        error_exit "请指定输出目录 (-o)"
    fi
    
    if [ ! -d "$FASTQ_DIR" ]; then
        error_exit "FASTQ目录不存在: $FASTQ_DIR"
    fi
    
    # 处理SILVA数据库
    if [ "$SKIP_TAXONOMY" != "true" ]; then
        # 如果没有通过参数提供，尝试从环境变量获取
        if [ -z "$SILVA_DB" ]; then
            SILVA_DB="${SILVA_DATABASE:-}"
        fi
        
        if [ -z "$SILVA_DB" ]; then
            warn "未指定SILVA数据库，将跳过物种注释"
            warn "使用 -d 参数指定数据库文件，或设置环境变量 SILVA_DATABASE"
            SKIP_TAXONOMY="true"
        else
            # 转换为绝对路径
            if [[ "$SILVA_DB" != /* ]]; then
                SILVA_DB="$(cd "$(dirname "$SILVA_DB")" 2>/dev/null && pwd)/$(basename "$SILVA_DB")"
            fi
            
            if [ ! -f "$SILVA_DB" ]; then
                warn "SILVA数据库文件不存在: $SILVA_DB"
                warn "将跳过物种注释"
                SKIP_TAXONOMY="true"
            fi
        fi
    fi
    
    # 显示配置
    log "========================================="
    log "后端基础分析流程"
    log "========================================="
    log "配置:"
    log "  FASTQ目录: $FASTQ_DIR"
    log "  输出目录: $OUTPUT_DIR"
    log "  每样本线程数: $THREADS"
    log "  并行任务数: $PARALLEL_JOBS"
    if [ "$SKIP_TAXONOMY" = "true" ]; then
        log "  物种注释: 跳过"
    else
        log "  SILVA数据库: $SILVA_DB"
    fi
    log "========================================="
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"/{single_samples,logs}
    
    # 导出变量（确保子进程能访问）
    export OUTPUT_DIR THREADS SILVA_DB SKIP_TAXONOMY
    export -f log warn error_exit process_single_sample
    
    # Step 1: 生成样本列表
    log "Step 1: 扫描和配对FASTQ文件"
    generate_sample_list "$FASTQ_DIR" "$OUTPUT_DIR"
    SAMPLE_LIST="$OUTPUT_DIR/sample_list.txt"
    
    # 统计样本数
    TOTAL_SAMPLES=$(wc -l < "$SAMPLE_LIST")
    log "总样本数: $TOTAL_SAMPLES"
    
    # Step 2: 处理单个样本
    log "Step 2: 处理单个样本"
    
    if [ "$PARALLEL_JOBS" -eq 1 ]; then
        # 串行处理
        log "串行处理模式"
        while IFS=$'\t' read -r sample r1 r2 group; do
            [ -z "$sample" ] && continue
            process_single_sample "$sample" "$r1" "$r2" "$group"
        done < "$SAMPLE_LIST"
    else
        # 并行处理
        run_parallel_processing "$SAMPLE_LIST" "$PARALLEL_JOBS"
    fi
    
    # Step 3: 合并结果
    log "Step 3: 合并所有样本结果"
    merge_asv_tables
    
    log "========================================="
    log "后端处理完成!"
    log "输出文件:"
    log "  - 合并表: $OUTPUT_DIR/merged_asv_taxonomy_table.tsv"
    log "  - 样本列表: $OUTPUT_DIR/sample_list.txt"
    log "  - 处理统计: $OUTPUT_DIR/processing_stats.json"
    log "========================================="
    
    # 清理临时目录
    if [ -d "$OUTPUT_DIR/tmp" ]; then
        log "清理临时文件..."
        rm -rf "$OUTPUT_DIR/tmp"
    fi
}

# 执行主函数
main "$@"
