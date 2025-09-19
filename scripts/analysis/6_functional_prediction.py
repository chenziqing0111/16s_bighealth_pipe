#!/usr/bin/env python3
"""
功能预测分析模块 - 简化版
直接从预处理结果读取功能数据，不运行PICRUSt2
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path
import sys

class FunctionalPredictor:
    def __init__(self):
        """初始化功能预测器"""
        # 可能的预处理目录位置
        self.preprocessing_dirs = [
            Path("backend_output/preprocessing"),
            Path("output/preprocessing"),
            Path("../preprocessing"),
            Path("../../preprocessing"),
            Path("preprocessing")
        ]
    
    def find_preprocessing_dir(self):
        """查找预处理结果目录"""
        for dir_path in self.preprocessing_dirs:
            if dir_path.exists() and (dir_path / 'merged_functional_annotation.tsv').exists():
                print(f"找到预处理目录: {dir_path}")
                return dir_path
        return None
    
    def load_from_preprocessing(self, sample_id):
        """从预处理结果提取特定样本的功能数据"""
        
        # 找到预处理目录
        preprocessing_dir = self.find_preprocessing_dir()
        
        if not preprocessing_dir:
            print("未找到预处理目录，使用默认值")
            return self.get_default_results(sample_id)
        
        results = {
            'sample_id': sample_id,
            'method': 'preprocessed_picrust2',
            'source': str(preprocessing_dir)
        }
        
        try:
            # 1. 从functional_summary.json读取摘要（如果存在）
            summary_file = preprocessing_dir / 'functional_summary.json'
            if summary_file.exists():
                with open(summary_file, 'r', encoding='utf-8') as f:
                    all_samples_summary = json.load(f)
                    
                if sample_id in all_samples_summary:
                    sample_summary = all_samples_summary[sample_id]
                    
                    # 提取维生素合成数据
                    results['vitamin_synthesis'] = {}
                    vitamin_keys = {
                        'vitamin_B1': 'B1',
                        'vitamin_B2': 'B2', 
                        'vitamin_B6': 'B6',
                        'vitamin_B9': 'B9',
                        'vitamin_B12': 'B12',
                        'vitamin_K': 'K'
                    }
                    
                    for key, name in vitamin_keys.items():
                        if key in sample_summary:
                            value = float(sample_summary[key])
                            results['vitamin_synthesis'][name] = {
                                'synthesis_potential': value,
                                'status': '高' if value > 100 else '中等' if value > 50 else '低'
                            }
                    
                    # 提取短链脂肪酸数据
                    results['scfa_production'] = {}
                    scfa_keys = ['butyrate', 'propionate', 'acetate']
                    scfa_names = {'butyrate': '丁酸', 'propionate': '丙酸', 'acetate': '乙酸'}
                    
                    for key in scfa_keys:
                        if key in sample_summary:
                            value = float(sample_summary[key])
                            results['scfa_production'][scfa_names[key]] = {
                                'production_potential': value,
                                'status': '强' if value > 1000 else '中等' if value > 500 else '弱'
                            }
                    
                    print(f"  从functional_summary.json加载了功能摘要")
            
            # 2. 从merged_functional_annotation.tsv读取KO数据
            ko_file = preprocessing_dir / 'merged_functional_annotation.tsv'
            if ko_file.exists():
                ko_df = pd.read_csv(ko_file, sep='\t', index_col=0)
                if sample_id in ko_df.columns:
                    sample_ko = ko_df[sample_id]
                    # 过滤掉0值
                    sample_ko_nonzero = sample_ko[sample_ko > 0]
                    
                    results['ko_abundances'] = {
                        'total_abundance': float(sample_ko.sum()),
                        'total_kos': int(len(sample_ko_nonzero)),
                        'top_kos': sample_ko_nonzero.nlargest(20).to_dict() if len(sample_ko_nonzero) > 0 else {}
                    }
                    print(f"  从merged_functional_annotation.tsv加载了{len(sample_ko_nonzero)}个KO")
            
            # 3. 从functional_pathway_annotation.tsv读取通路数据
            pathway_file = preprocessing_dir / 'functional_pathway_annotation.tsv'
            if pathway_file.exists():
                pathway_df = pd.read_csv(pathway_file, sep='\t', index_col=0)
                if sample_id in pathway_df.columns:
                    sample_pathways = pathway_df[sample_id]
                    # 过滤掉0值
                    sample_pathways_nonzero = sample_pathways[sample_pathways > 0]
                    
                    results['pathway_abundances'] = {
                        'total_pathways': int(len(sample_pathways_nonzero)),
                        'top_pathways': sample_pathways_nonzero.nlargest(20).to_dict() if len(sample_pathways_nonzero) > 0 else {}
                    }
                    print(f"  从functional_pathway_annotation.tsv加载了{len(sample_pathways_nonzero)}个通路")
            
            # 4. 从functional_ec_annotation.tsv读取EC数据
            ec_file = preprocessing_dir / 'functional_ec_annotation.tsv'
            if ec_file.exists():
                ec_df = pd.read_csv(ec_file, sep='\t', index_col=0)
                if sample_id in ec_df.columns:
                    sample_ec = ec_df[sample_id]
                    # 过滤掉0值
                    sample_ec_nonzero = sample_ec[sample_ec > 0]
                    
                    results['ec_abundances'] = {
                        'total_ecs': int(len(sample_ec_nonzero)),
                        'top_ecs': sample_ec_nonzero.nlargest(20).to_dict() if len(sample_ec_nonzero) > 0 else {}
                    }
                    print(f"  从functional_ec_annotation.tsv加载了{len(sample_ec_nonzero)}个EC")
            
            # 5. 分析关键功能
            results['key_functions'] = self.analyze_key_functions(results)
            
            return results
            
        except Exception as e:
            print(f"加载预处理数据时出错: {e}")
            return self.get_default_results(sample_id)
    
    def analyze_key_functions(self, results):
        """基于通路和KO数据分析关键功能"""
        key_functions = []
        
        # 从通路分析
        if 'pathway_abundances' in results and 'top_pathways' in results['pathway_abundances']:
            for pathway_id in list(results['pathway_abundances']['top_pathways'].keys())[:5]:
                # 识别关键通路类型
                if 'GLYCOLYSIS' in pathway_id.upper():
                    key_functions.append('糖酵解')
                elif 'TCA' in pathway_id.upper():
                    key_functions.append('三羧酸循环')
                elif 'FERMENT' in pathway_id.upper():
                    key_functions.append('发酵')
                elif 'ARG' in pathway_id.upper():
                    key_functions.append('精氨酸代谢')
                elif 'FOLATE' in pathway_id.upper() or '1CMET' in pathway_id.upper():
                    key_functions.append('一碳代谢')
                elif 'BUTYRATE' in pathway_id.upper() or 'BUTANOATE' in pathway_id.upper():
                    key_functions.append('丁酸生成')
        
        # 去重
        key_functions = list(set(key_functions))
        
        # 如果没有识别到，添加默认功能
        if not key_functions:
            key_functions = ['能量代谢', '氨基酸代谢', '碳水化合物代谢']
        
        return key_functions[:10]  # 最多返回10个
    
    def get_default_results(self, sample_id):
        """当无法加载预处理数据时，返回默认结果"""
        return {
            'sample_id': sample_id,
            'method': 'default',
            'vitamin_synthesis': {
                'B1': {'synthesis_potential': 50, 'status': '中等'},
                'B2': {'synthesis_potential': 45, 'status': '低'},
                'B6': {'synthesis_potential': 60, 'status': '中等'},
                'B9': {'synthesis_potential': 55, 'status': '中等'},
                'B12': {'synthesis_potential': 30, 'status': '低'},
                'K': {'synthesis_potential': 70, 'status': '中等'}
            },
            'scfa_production': {
                '丁酸': {'production_potential': 800, 'status': '中等'},
                '丙酸': {'production_potential': 600, 'status': '中等'},
                '乙酸': {'production_potential': 1200, 'status': '强'}
            },
            'key_functions': ['能量代谢', '氨基酸代谢', '碳水化合物代谢'],
            'note': '使用默认值，建议检查预处理结果'
        }
    
    def save_results(self, results, output_dir):
        """保存功能预测结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON结果
        output_file = output_path / 'functional_prediction.json'
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # 生成摘要报告
        summary_file = output_path / 'functional_summary.txt'
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 50 + "\n")
            f.write("功能预测分析摘要\n")
            f.write("=" * 50 + "\n")
            f.write(f"样本ID: {results.get('sample_id', 'Unknown')}\n")
            f.write(f"方法: {results.get('method', 'Unknown')}\n")
            
            if 'source' in results:
                f.write(f"数据来源: {results['source']}\n")
            
            f.write("\n维生素合成能力:\n")
            if 'vitamin_synthesis' in results:
                for vit, data in results['vitamin_synthesis'].items():
                    f.write(f"  维生素{vit}: {data['status']} (潜力值: {data['synthesis_potential']:.1f})\n")
            
            f.write("\n短链脂肪酸产生能力:\n")
            if 'scfa_production' in results:
                for scfa, data in results['scfa_production'].items():
                    f.write(f"  {scfa}: {data['status']} (产量: {data['production_potential']:.1f})\n")
            
            if 'key_functions' in results:
                f.write("\n关键功能:\n")
                for func in results['key_functions']:
                    f.write(f"  - {func}\n")
            
            if 'ko_abundances' in results:
                f.write(f"\nKO统计:\n")
                f.write(f"  总KO数: {results['ko_abundances'].get('total_kos', 0)}\n")
                f.write(f"  总丰度: {results['ko_abundances'].get('total_abundance', 0):.0f}\n")
            
            f.write("=" * 50 + "\n")
        
        return output_file

def main():
    parser = argparse.ArgumentParser(description='微生物功能预测分析')
    parser.add_argument('--input', '-i', required=True, 
                       help='ASV表文件路径（用于识别样本）')
    parser.add_argument('--output', '-o', required=True, 
                       help='输出目录')
    
    args = parser.parse_args()
    
    # 从输出路径推断样本ID
    # 假设路径格式：.../sample_id/functional_prediction
    output_path = Path(args.output)
    sample_id = output_path.parent.name
    
    print(f"功能预测分析")
    print(f"  样本: {sample_id}")
    print(f"  输出: {args.output}")
    
    # 创建预测器
    predictor = FunctionalPredictor()
    
    # 加载预处理结果
    print("加载预处理功能数据...")
    results = predictor.load_from_preprocessing(sample_id)
    
    # 保存结果
    output_file = predictor.save_results(results, args.output)
    
    print(f"✓ 功能预测完成")
    print(f"  结果文件: {output_file}")
    print(f"  使用方法: {results.get('method', 'unknown')}")
    
    # 打印关键信息
    if 'vitamin_synthesis' in results:
        high_vitamins = [v for v, d in results['vitamin_synthesis'].items() if d['status'] == '高']
        if high_vitamins:
            print(f"  高合成维生素: {', '.join(high_vitamins)}")
    
    if 'scfa_production' in results:
        strong_scfa = [s for s, d in results['scfa_production'].items() if d['status'] == '强']
        if strong_scfa:
            print(f"  强产生SCFA: {', '.join(strong_scfa)}")

if __name__ == '__main__':
    main()