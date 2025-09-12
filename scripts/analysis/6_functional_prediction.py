#!/usr/bin/env python3
"""
功能预测分析模块
使用 PICRUSt2 预测微生物群落的功能潜力
"""

import pandas as pd
import numpy as np
import json
import argparse
import subprocess
import os
from pathlib import Path
import tempfile
import shutil

class FunctionalPredictor:
    def __init__(self):
        """初始化功能预测器"""
        self.kegg_pathways = {
            # 维生素合成相关
            'ko00730': 'Thiamine metabolism (维生素B1)',
            'ko00740': 'Riboflavin metabolism (维生素B2)', 
            'ko00750': 'Vitamin B6 metabolism',
            'ko00760': 'Nicotinate metabolism (维生素B3)',
            'ko00770': 'Pantothenate metabolism (维生素B5)',
            'ko00780': 'Biotin metabolism (维生素B7)',
            'ko00790': 'Folate biosynthesis (叶酸)',
            'ko00860': 'Vitamin B12 metabolism',
            'ko00130': 'Vitamin K biosynthesis',
            
            # 短链脂肪酸
            'ko00650': 'Butanoate metabolism (丁酸)',
            'ko00640': 'Propanoate metabolism (丙酸)',
            'ko00620': 'Pyruvate metabolism (乙酸前体)',
            
            # 氨基酸代谢
            'ko00250': 'Alanine metabolism',
            'ko00260': 'Glycine and serine metabolism',
            'ko00270': 'Cysteine and methionine metabolism',
            'ko00280': 'Valine, leucine and isoleucine',
            'ko00290': 'Valine, leucine and isoleucine biosynthesis',
            'ko00300': 'Lysine biosynthesis',
            'ko00310': 'Lysine degradation',
            'ko00330': 'Arginine and proline metabolism',
            'ko00340': 'Histidine metabolism',
            'ko00350': 'Tyrosine metabolism',
            'ko00360': 'Phenylalanine metabolism',
            'ko00380': 'Tryptophan metabolism',
            
            # 其他重要功能
            'ko00190': 'Oxidative phosphorylation',
            'ko00195': 'Photosynthesis',
            'ko00680': 'Methane metabolism',
            'ko00910': 'Nitrogen metabolism',
            'ko00920': 'Sulfur metabolism'
        }
        
        self.vitamin_pathways = {
            'B1': ['ko00730'],
            'B2': ['ko00740'],
            'B3': ['ko00760'],
            'B5': ['ko00770'],
            'B6': ['ko00750'],
            'B7': ['ko00780'],
            'B9': ['ko00790'],
            'B12': ['ko00860'],
            'K': ['ko00130']
        }
        
        self.scfa_pathways = {
            'Butyrate': ['ko00650'],
            'Propionate': ['ko00640'],
            'Acetate': ['ko00620']
        }
    
    def prepare_input_for_picrust2(self, asv_table_path, temp_dir):
        """准备 PICRUSt2 输入文件"""
        print("  准备 PICRUSt2 输入...")
        
        # 读取ASV表
        df = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        
        # 分离丰度数据和序列信息
        abundance_cols = []
        for col in df.columns:
            if col not in ['Taxon', 'Confidence', 'sequence']:
                try:
                    # 检查是否是数值列
                    pd.to_numeric(df[col])
                    abundance_cols.append(col)
                except:
                    pass
        
        if not abundance_cols:
            raise ValueError("未找到丰度数据列")
        
        # 创建特征表（BIOM格式的替代）
        feature_table = df[abundance_cols].copy()
        feature_table.index.name = '#OTU ID'
        
        # 保存特征表
        feature_table_path = Path(temp_dir) / 'feature_table.tsv'
        feature_table.to_csv(feature_table_path, sep='\t')
        
        # 如果有序列信息，创建FASTA文件
        fasta_path = None
        if 'sequence' in df.columns:
            fasta_path = Path(temp_dir) / 'sequences.fasta'
            with open(fasta_path, 'w') as f:
                for asv_id, seq in zip(df.index, df['sequence']):
                    if pd.notna(seq):
                        f.write(f">{asv_id}\n{seq}\n")
        
        return feature_table_path, fasta_path, abundance_cols
    
    def run_picrust2(self, feature_table_path, fasta_path, output_dir):
        """运行 PICRUSt2 pipeline"""
        print("  运行 PICRUSt2...")
        
        try:
            # 构建命令
            cmd = [
                'picrust2_pipeline.py',
                '-s', str(fasta_path) if fasta_path else '',
                '-i', str(feature_table_path),
                '-o', str(output_dir),
                '-p', '4',  # 线程数
                '--stratified',  # 分层输出
                '--verbose'
            ]
            
            # 如果没有序列文件，使用默认参考
            if not fasta_path:
                cmd.extend(['--skip_place_seqs'])
            
            # 运行命令
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            print("  ✓ PICRUSt2 运行成功")
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"  ✗ PICRUSt2 运行失败: {e}")
            print(f"  错误输出: {e.stderr}")
            return False
        except FileNotFoundError:
            print("  ✗ 未找到 PICRUSt2，尝试简化分析...")
            return False
    
    def parse_picrust2_results(self, output_dir):
        """解析 PICRUSt2 结果"""
        results = {}
        
        # 读取 KO 预测结果
        ko_file = Path(output_dir) / 'KO_metagenome_out' / 'pred_metagenome_unstrat.tsv'
        if ko_file.exists():
            ko_df = pd.read_csv(ko_file, sep='\t', index_col=0)
            results['ko_abundances'] = ko_df.to_dict()
        
        # 读取通路预测结果
        pathway_file = Path(output_dir) / 'pathways_out' / 'path_abun_unstrat.tsv'
        if pathway_file.exists():
            pathway_df = pd.read_csv(pathway_file, sep='\t', index_col=0)
            results['pathway_abundances'] = pathway_df.to_dict()
        
        # 读取EC预测结果
        ec_file = Path(output_dir) / 'EC_metagenome_out' / 'pred_metagenome_unstrat.tsv'
        if ec_file.exists():
            ec_df = pd.read_csv(ec_file, sep='\t', index_col=0)
            results['ec_abundances'] = ec_df.to_dict()
        
        return results
    
    def analyze_vitamin_synthesis(self, ko_abundances):
        """分析维生素合成能力"""
        vitamin_scores = {}
        
        for vitamin, pathways in self.vitamin_pathways.items():
            pathway_scores = []
            for pathway in pathways:
                if pathway in ko_abundances:
                    # 计算相对丰度
                    abundance = ko_abundances[pathway]
                    if isinstance(abundance, dict):
                        abundance = sum(abundance.values()) / len(abundance)
                    pathway_scores.append(abundance)
            
            if pathway_scores:
                vitamin_scores[vitamin] = {
                    'synthesis_potential': np.mean(pathway_scores),
                    'status': self._get_synthesis_status(np.mean(pathway_scores))
                }
        
        return vitamin_scores
    
    def analyze_scfa_production(self, ko_abundances):
        """分析短链脂肪酸产生能力"""
        scfa_scores = {}
        
        for scfa, pathways in self.scfa_pathways.items():
            pathway_scores = []
            for pathway in pathways:
                if pathway in ko_abundances:
                    abundance = ko_abundances[pathway]
                    if isinstance(abundance, dict):
                        abundance = sum(abundance.values()) / len(abundance)
                    pathway_scores.append(abundance)
            
            if pathway_scores:
                scfa_scores[scfa] = {
                    'production_potential': np.mean(pathway_scores),
                    'status': self._get_production_status(np.mean(pathway_scores))
                }
        
        return scfa_scores
    
    def _get_synthesis_status(self, score):
        """评估合成能力状态"""
        if score > 1000:
            return '高'
        elif score > 100:
            return '中等'
        else:
            return '低'
    
    def _get_production_status(self, score):
        """评估产生能力状态"""
        if score > 5000:
            return '强'
        elif score > 1000:
            return '中等'
        else:
            return '弱'
    
    def simplified_functional_prediction(self, asv_table_path):
        """简化的功能预测（不依赖PICRUSt2）"""
        print("  执行简化功能预测...")
        
        # 读取ASV表
        df = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        
        # 基于已知菌属的功能特征进行预测
        functional_bacteria = {
            'Bifidobacterium': {
                'vitamins': ['B1', 'B2', 'B9', 'K'],
                'scfa': ['Acetate'],
                'functions': ['益生菌', '免疫调节']
            },
            'Lactobacillus': {
                'vitamins': ['B2', 'B6', 'B12'],
                'scfa': ['Lactate'],
                'functions': ['益生菌', '乳酸产生']
            },
            'Faecalibacterium': {
                'vitamins': ['B2'],
                'scfa': ['Butyrate'],
                'functions': ['丁酸产生', '抗炎']
            },
            'Bacteroides': {
                'vitamins': ['B7', 'K'],
                'scfa': ['Propionate', 'Acetate'],
                'functions': ['多糖降解', '短链脂肪酸产生']
            },
            'Akkermansia': {
                'vitamins': [],
                'scfa': ['Propionate'],
                'functions': ['黏蛋白降解', '代谢调节']
            },
            'Roseburia': {
                'vitamins': [],
                'scfa': ['Butyrate'],
                'functions': ['丁酸产生']
            },
            'Prevotella': {
                'vitamins': ['B1'],
                'scfa': ['Propionate', 'Acetate'],
                'functions': ['纤维降解']
            }
        }
        
        # 提取菌属丰度
        genus_abundance = {}
        if 'Taxon' in df.columns:
            for idx, row in df.iterrows():
                if pd.notna(row['Taxon']):
                    # 提取属名
                    taxon_parts = row['Taxon'].split(';')
                    for part in taxon_parts:
                        if part.startswith('g__'):
                            genus = part.replace('g__', '')
                            if genus and genus != 'unclassified':
                                # 获取丰度
                                abundance_cols = [col for col in df.columns 
                                                if col not in ['Taxon', 'Confidence']]
                                if abundance_cols:
                                    abundance = row[abundance_cols[0]]
                                    if genus not in genus_abundance:
                                        genus_abundance[genus] = 0
                                    genus_abundance[genus] += abundance
        
        # 基于菌属预测功能
        vitamin_potential = {f'B{i}': 0 for i in range(1, 13)}
        vitamin_potential['K'] = 0
        scfa_potential = {'Butyrate': 0, 'Propionate': 0, 'Acetate': 0}
        key_functions = []
        
        for genus, abundance in genus_abundance.items():
            if genus in functional_bacteria:
                features = functional_bacteria[genus]
                
                # 维生素合成
                for vitamin in features['vitamins']:
                    vitamin_potential[vitamin] += abundance
                
                # 短链脂肪酸
                for scfa in features['scfa']:
                    if scfa in scfa_potential:
                        scfa_potential[scfa] += abundance
                
                # 关键功能
                if abundance > 100:  # 丰度阈值
                    key_functions.extend(features['functions'])
        
        # 标准化评分
        total_abundance = sum(genus_abundance.values()) if genus_abundance else 1
        
        vitamin_scores = {}
        for vitamin, score in vitamin_potential.items():
            normalized_score = (score / total_abundance) * 100
            vitamin_scores[vitamin] = {
                'synthesis_potential': normalized_score,
                'status': '高' if normalized_score > 5 else '中等' if normalized_score > 1 else '低'
            }
        
        scfa_scores = {}
        for scfa, score in scfa_potential.items():
            normalized_score = (score / total_abundance) * 100
            scfa_scores[scfa] = {
                'production_potential': normalized_score,
                'status': '强' if normalized_score > 10 else '中等' if normalized_score > 5 else '弱'
            }
        
        return {
            'vitamin_synthesis': vitamin_scores,
            'scfa_production': scfa_scores,
            'key_functions': list(set(key_functions)),
            'method': 'simplified'
        }
    
    def predict(self, asv_table_path, output_dir):
        """执行功能预测"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # 创建临时目录
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                # 准备输入
                feature_table, fasta, samples = self.prepare_input_for_picrust2(
                    asv_table_path, temp_dir
                )
                
                # 创建PICRUSt2输出目录
                picrust2_output = Path(temp_dir) / 'picrust2_out'
                
                # 运行PICRUSt2
                if self.run_picrust2(feature_table, fasta, picrust2_output):
                    # 解析结果
                    picrust2_results = self.parse_picrust2_results(picrust2_output)
                    
                    if 'ko_abundances' in picrust2_results:
                        # 分析维生素合成
                        results['vitamin_synthesis'] = self.analyze_vitamin_synthesis(
                            picrust2_results['ko_abundances']
                        )
                        
                        # 分析短链脂肪酸
                        results['scfa_production'] = self.analyze_scfa_production(
                            picrust2_results['ko_abundances']
                        )
                        
                        results['method'] = 'picrust2'
                        results['raw_predictions'] = picrust2_results
                    else:
                        # PICRUSt2运行但无结果，使用简化方法
                        results = self.simplified_functional_prediction(asv_table_path)
                else:
                    # PICRUSt2失败，使用简化方法
                    results = self.simplified_functional_prediction(asv_table_path)
                    
            except Exception as e:
                print(f"  功能预测出错: {e}")
                # 使用简化方法作为后备
                results = self.simplified_functional_prediction(asv_table_path)
        
        # 添加总结
        results['summary'] = self._generate_summary(results)
        
        # 保存结果
        output_file = output_path / 'functional_prediction.json'
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        return results
    
    def _generate_summary(self, results):
        """生成功能预测总结"""
        summary = {
            'method_used': results.get('method', 'unknown'),
            'key_findings': []
        }
        
        # 维生素合成能力总结
        if 'vitamin_synthesis' in results:
            high_vitamins = [v for v, d in results['vitamin_synthesis'].items() 
                           if d['status'] == '高']
            if high_vitamins:
                summary['key_findings'].append(
                    f"高合成能力维生素: {', '.join(high_vitamins)}"
                )
        
        # 短链脂肪酸产生能力总结
        if 'scfa_production' in results:
            strong_scfa = [s for s, d in results['scfa_production'].items() 
                         if d['status'] == '强']
            if strong_scfa:
                summary['key_findings'].append(
                    f"强产生能力短链脂肪酸: {', '.join(strong_scfa)}"
                )
        
        # 关键功能
        if 'key_functions' in results and results['key_functions']:
            summary['key_findings'].append(
                f"关键功能: {', '.join(results['key_functions'][:5])}"
            )
        
        return summary

def main():
    parser = argparse.ArgumentParser(description='微生物功能预测分析')
    parser.add_argument('--input', '-i', required=True, help='ASV表文件路径')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    # 执行预测
    predictor = FunctionalPredictor()
    results = predictor.predict(args.input, args.output)
    
    # 打印结果摘要
    print(f"功能预测完成，结果保存至: {args.output}")
    print(f"  使用方法: {results.get('method', 'unknown')}")
    
    if 'summary' in results:
        for finding in results['summary'].get('key_findings', []):
            print(f"  - {finding}")

if __name__ == '__main__':
    main()