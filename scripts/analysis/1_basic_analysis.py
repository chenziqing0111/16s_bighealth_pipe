#!/usr/bin/env python3
"""
基础分析：多样性、物种组成、B/F比值
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

class BasicAnalyzer:
    def __init__(self, asv_table_path):
        """初始化分析器"""
        self.asv_table = self._load_asv_table(asv_table_path)
        self.sample_id = self._get_sample_id()
        self.results = {}
        
    def _load_asv_table(self, path):
        """加载ASV表"""
        df = pd.read_csv(path, sep='\t', index_col=0)
        return df
    
    def _get_sample_id(self):
        """获取样本ID（第一个非分类列）"""
        tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        sample_cols = [col for col in self.asv_table.columns if col not in tax_cols]
        return sample_cols[0] if sample_cols else None
    
    def calculate_alpha_diversity(self):
        """计算Alpha多样性"""
        if not self.sample_id:
            return
            
        counts = self.asv_table[self.sample_id].values
        counts = counts[counts > 0]  # 去除零值
        
        if len(counts) == 0:
            return
        
        # Shannon指数
        proportions = counts / counts.sum()
        shannon = -np.sum(proportions * np.log(proportions))
        
        # Simpson指数
        simpson = 1 - np.sum(proportions ** 2)
        
        # Chao1估计
        observed_asvs = len(counts)
        singletons = np.sum(counts == 1)
        doubletons = np.sum(counts == 2)
        
        if doubletons > 0:
            chao1 = observed_asvs + (singletons ** 2) / (2 * doubletons)
        else:
            chao1 = observed_asvs + singletons * (singletons - 1) / 2 if singletons > 1 else observed_asvs
        
        # Pielou均匀度
        evenness = shannon / np.log(observed_asvs) if observed_asvs > 1 else 1
        
        self.results['alpha_diversity'] = {
            'shannon': float(shannon),
            'simpson': float(simpson),
            'chao1': float(chao1),
            'observed_asvs': int(observed_asvs),
            'evenness': float(evenness),
            'total_reads': int(self.asv_table[self.sample_id].sum())
        }
        
        # 评估多样性水平
        if 350 <= observed_asvs <= 770:
            diversity_status = '正常'
        elif observed_asvs < 350:
            diversity_status = '偏低'
        else:
            diversity_status = '偏高'
            
        self.results['alpha_diversity']['status'] = diversity_status
    
    def calculate_bf_ratio(self):
        """计算拟杆菌门/厚壁菌门比值"""
        if 'Phylum' not in self.asv_table.columns:
            self.results['bf_ratio'] = {'value': None, 'status': '无法计算'}
            return
        
        # 按门水平汇总
        phylum_abundance = {}
        for idx, row in self.asv_table.iterrows():
            phylum = row['Phylum']
            if pd.notna(phylum) and phylum != 'Unclassified':
                # 清理门名称
                phylum = phylum.replace('p__', '').strip()
                if phylum not in phylum_abundance:
                    phylum_abundance[phylum] = 0
                phylum_abundance[phylum] += row[self.sample_id]
        
        # 计算B/F比值
        bacteroidetes = phylum_abundance.get('Bacteroidetes', 0) + phylum_abundance.get('Bacteroidota', 0)
        firmicutes = phylum_abundance.get('Firmicutes', 0) + phylum_abundance.get('Bacillota', 0)
        
        if firmicutes > 0:
            bf_ratio = bacteroidetes / firmicutes
            
            # 评估状态（参考范围：0.84-4.94）
            if 0.84 <= bf_ratio <= 4.94:
                status = '正常'
            elif bf_ratio < 0.84:
                status = '偏低（厚壁菌门占优）'
            else:
                status = '偏高（拟杆菌门占优）'
        else:
            bf_ratio = None
            status = '无法计算'
        
        self.results['bf_ratio'] = {
            'value': float(bf_ratio) if bf_ratio else None,
            'bacteroidetes': int(bacteroidetes),
            'firmicutes': int(firmicutes),
            'status': status
        }
    
    def analyze_composition(self):
        """分析物种组成"""
        composition = {}
        
        # 分析各分类水平
        for level in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
            if level not in self.asv_table.columns:
                continue
                
            level_abundance = {}
            total_reads = 0
            
            for idx, row in self.asv_table.iterrows():
                taxon = row[level]
                if pd.notna(taxon) and taxon != 'Unclassified':
                    # 清理分类名称
                    taxon = taxon.replace(f'{level[0].lower()}__', '').strip()
                    if taxon not in level_abundance:
                        level_abundance[taxon] = 0
                    level_abundance[taxon] += row[self.sample_id]
                    total_reads += row[self.sample_id]
            
            # 转换为相对丰度并排序
            if total_reads > 0:
                for taxon in level_abundance:
                    level_abundance[taxon] = (level_abundance[taxon] / total_reads) * 100
                
                # 只保留Top 10
                sorted_taxa = sorted(level_abundance.items(), key=lambda x: x[1], reverse=True)[:10]
                composition[level.lower()] = {
                    'taxa': [t[0] for t in sorted_taxa],
                    'abundance': [round(t[1], 3) for t in sorted_taxa]
                }
        
        self.results['composition'] = composition
    
    def calculate_core_stats(self):
        """计算核心统计信息"""
        if not self.sample_id:
            return
            
        # 基本统计
        total_reads = int(self.asv_table[self.sample_id].sum())
        total_asvs = int((self.asv_table[self.sample_id] > 0).sum())
        
        # ASV丰度分布
        asv_abundance = self.asv_table[self.sample_id][self.asv_table[self.sample_id] > 0]
        
        self.results['basic_stats'] = {
            'total_reads': total_reads,
            'total_asvs': total_asvs,
            'singleton_asvs': int((asv_abundance == 1).sum()),
            'mean_asv_abundance': float(asv_abundance.mean()) if len(asv_abundance) > 0 else 0,
            'max_asv_abundance': int(asv_abundance.max()) if len(asv_abundance) > 0 else 0
        }
    
    def save_results(self, output_dir):
        """保存结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON
        with open(output_path / 'basic_analysis.json', 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        # 保存多样性表格
        if 'alpha_diversity' in self.results:
            alpha_df = pd.DataFrame([self.results['alpha_diversity']])
            alpha_df.to_csv(output_path / 'alpha_diversity.tsv', sep='\t', index=False)
        
        print(f"基础分析完成，结果保存至: {output_path}")
        
        # 打印摘要
        if 'alpha_diversity' in self.results:
            print(f"  Shannon指数: {self.results['alpha_diversity']['shannon']:.3f}")
            print(f"  观察到的ASVs: {self.results['alpha_diversity']['observed_asvs']}")
            print(f"  多样性状态: {self.results['alpha_diversity']['status']}")
        
        if 'bf_ratio' in self.results and self.results['bf_ratio']['value']:
            print(f"  B/F比值: {self.results['bf_ratio']['value']:.2f}")
            print(f"  B/F状态: {self.results['bf_ratio']['status']}")

def main():
    parser = argparse.ArgumentParser(description='基础分析：多样性、物种组成、B/F比值')
    parser.add_argument('--input', '-i', required=True, help='ASV表路径')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    analyzer = BasicAnalyzer(args.input)
    analyzer.calculate_alpha_diversity()
    analyzer.calculate_bf_ratio()
    analyzer.analyze_composition()
    analyzer.calculate_core_stats()
    analyzer.save_results(args.output)

if __name__ == '__main__':
    main()
