#!/usr/bin/env python3
"""
菌群健康评估：有益菌、有害菌、条件致病菌评分
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path

class BacteriaEvaluator:
    def __init__(self, asv_table_path, ranges_path):
        """初始化评估器"""
        self.asv_table = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        self.normal_ranges = self._load_ranges(ranges_path)
        self.sample_id = self._get_sample_id()
        self.results = {}
        
    def _load_ranges(self, path):
        """加载正常值范围"""
        if path and Path(path).exists():
            with open(path, 'r', encoding='utf-8') as f:
                return json.load(f)
        else:
            # 默认范围
            return self._get_default_ranges()
    
    def _get_default_ranges(self):
        """默认正常值范围"""
        return {
            "beneficial": {
                "Bifidobacterium": [0.183, 14.6],
                "Lactobacillus": [0.008, 0.03],
                "Faecalibacterium": [0.636, 10.97],
                "Akkermansia": [0.01, 3.89],
                "Prevotella": [0.014, 66.13],
                "Roseburia": [0.1, 5.0],
                "Coprococcus": [0.009, 0.57],
                "Butyricimonas": [0.014, 0.86],
                "Odoribacter": [0.013, 0.61],
                "Alistipes": [0.01, 2.89]
            },
            "harmful": {
                "Escherichia": [0, 0.5],
                "Shigella": [0, 0.14],
                "Salmonella": [0, 0.01],
                "Clostridium_difficile": [0, 0.001],
                "Staphylococcus": [0, 0.1],
                "Klebsiella": [0, 0.11],
                "Enterococcus": [0, 0.5],
                "Fusobacterium": [0, 0.1],
                "Campylobacter": [0, 0.01],
                "Helicobacter": [0, 0.01]
            },
            "conditional": {
                "Veillonella": [0, 0.22],
                "Streptococcus": [0, 2.0],
                "Bacteroides": [0, 30.0],
                "Eggerthella": [0, 0.1],
                "Peptostreptococcus": [0, 0.1],
                "Haemophilus": [0, 0.5]
            }
        }
    
    def _get_sample_id(self):
        """获取样本ID"""
        tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        sample_cols = [col for col in self.asv_table.columns if col not in tax_cols]
        return sample_cols[0] if sample_cols else None
    
    def _get_bacteria_abundance(self, bacteria_name, level='Genus'):
        """获取特定细菌的相对丰度"""
        if level not in self.asv_table.columns:
            return 0
            
        total_abundance = 0
        total_reads = self.asv_table[self.sample_id].sum()
        
        for idx, row in self.asv_table.iterrows():
            taxon = str(row[level]) if pd.notna(row[level]) else ''
            # 清理分类名称
            taxon = taxon.replace(f'{level[0].lower()}__', '').strip()
            
            # 模糊匹配（处理不同的命名方式）
            if bacteria_name.lower() in taxon.lower() or taxon.lower() in bacteria_name.lower():
                total_abundance += row[self.sample_id]
        
        # 返回相对丰度（百分比）
        return (total_abundance / total_reads * 100) if total_reads > 0 else 0
    
    def evaluate_beneficial_bacteria(self):
        """评估有益菌"""
        beneficial_results = {}
        total_score = 0
        max_score = 0
        
        for bacteria, normal_range in self.normal_ranges['beneficial'].items():
            abundance = self._get_bacteria_abundance(bacteria)
            min_val, max_val = normal_range
            
            # 评估状态
            if min_val <= abundance <= max_val:
                status = '正常'
                score = 10
            elif abundance < min_val:
                status = '偏低'
                score = 5 * (abundance / min_val) if min_val > 0 else 0
            else:
                status = '偏高'
                score = 8  # 有益菌偏高通常不是大问题
            
            beneficial_results[bacteria] = {
                'abundance': round(abundance, 4),
                'normal_range': normal_range,
                'status': status,
                'score': round(score, 2)
            }
            
            total_score += score
            max_score += 10
        
        # 计算总体评分
        overall_score = (total_score / max_score * 100) if max_score > 0 else 0
        
        self.results['beneficial_bacteria'] = {
            'bacteria': beneficial_results,
            'overall_score': round(overall_score, 1),
            'evaluation': self._get_evaluation(overall_score)
        }
    
    def evaluate_harmful_bacteria(self):
        """评估有害菌"""
        harmful_results = {}
        total_penalty = 0
        
        for bacteria, normal_range in self.normal_ranges['harmful'].items():
            abundance = self._get_bacteria_abundance(bacteria)
            min_val, max_val = normal_range
            
            # 评估状态
            if abundance <= max_val:
                status = '正常'
                penalty = 0
            else:
                status = '超标'
                # 超标越多，扣分越多
                penalty = min(20, (abundance - max_val) / max_val * 10) if max_val > 0 else 20
            
            harmful_results[bacteria] = {
                'abundance': round(abundance, 4),
                'threshold': max_val,
                'status': status,
                'penalty': round(penalty, 2)
            }
            
            total_penalty += penalty
        
        # 计算危害评分（100分制，扣分制）
        harm_score = max(0, 100 - total_penalty)
        
        self.results['harmful_bacteria'] = {
            'bacteria': harmful_results,
            'harm_score': round(harm_score, 1),
            'risk_level': self._get_risk_level(harm_score)
        }
    
    def evaluate_conditional_bacteria(self):
        """评估条件致病菌"""
        conditional_results = {}
        warning_count = 0
        
        for bacteria, normal_range in self.normal_ranges['conditional'].items():
            abundance = self._get_bacteria_abundance(bacteria)
            min_val, max_val = normal_range
            
            # 评估状态
            if abundance <= max_val:
                status = '正常'
            else:
                status = '需关注'
                warning_count += 1
            
            conditional_results[bacteria] = {
                'abundance': round(abundance, 4),
                'threshold': max_val,
                'status': status
            }
        
        self.results['conditional_bacteria'] = {
            'bacteria': conditional_results,
            'warning_count': warning_count,
            'attention_needed': warning_count > 3
        }
    
    def calculate_overall_health_score(self):
        """计算整体菌群健康评分"""
        beneficial_score = self.results.get('beneficial_bacteria', {}).get('overall_score', 50)
        harm_score = self.results.get('harmful_bacteria', {}).get('harm_score', 50)
        conditional_penalty = self.results.get('conditional_bacteria', {}).get('warning_count', 0) * 2
        
        # 综合评分（权重：有益菌40%，有害菌40%，条件致病菌20%）
        overall_score = (beneficial_score * 0.4 + 
                        harm_score * 0.4 + 
                        max(0, 100 - conditional_penalty * 5) * 0.2)
        
        self.results['overall_health'] = {
            'score': round(overall_score, 1),
            'grade': self._get_health_grade(overall_score),
            'components': {
                'beneficial_contribution': round(beneficial_score * 0.4, 1),
                'harmful_contribution': round(harm_score * 0.4, 1),
                'conditional_contribution': round(max(0, 100 - conditional_penalty * 5) * 0.2, 1)
            }
        }
    
    def _get_evaluation(self, score):
        """获取评价"""
        if score >= 80:
            return '优秀'
        elif score >= 60:
            return '良好'
        elif score >= 40:
            return '一般'
        else:
            return '较差'
    
    def _get_risk_level(self, score):
        """获取风险等级"""
        if score >= 80:
            return '低风险'
        elif score >= 60:
            return '中低风险'
        elif score >= 40:
            return '中高风险'
        else:
            return '高风险'
    
    def _get_health_grade(self, score):
        """获取健康等级"""
        if score >= 90:
            return 'A+ (优秀)'
        elif score >= 80:
            return 'A (良好)'
        elif score >= 70:
            return 'B (正常)'
        elif score >= 60:
            return 'C (亚健康)'
        else:
            return 'D (需改善)'
    
    def generate_recommendations(self):
        """生成改善建议"""
        recommendations = []
        
        # 基于有益菌状态
        if self.results['beneficial_bacteria']['overall_score'] < 60:
            recommendations.append({
                'category': '益生菌补充',
                'suggestion': '建议补充含双歧杆菌和乳酸杆菌的益生菌制剂',
                'foods': ['酸奶', '发酵食品', '泡菜', '味噌']
            })
        
        # 基于有害菌状态
        if self.results['harmful_bacteria']['harm_score'] < 70:
            recommendations.append({
                'category': '抑制有害菌',
                'suggestion': '增加膳食纤维摄入，减少高脂高糖食物',
                'foods': ['全谷物', '蔬菜', '水果', '豆类']
            })
        
        # 基于条件致病菌
        if self.results['conditional_bacteria']['attention_needed']:
            recommendations.append({
                'category': '免疫调节',
                'suggestion': '增强免疫力，保持肠道平衡',
                'foods': ['大蒜', '生姜', '绿茶', '蘑菇类']
            })
        
        self.results['recommendations'] = recommendations
    
    def save_results(self, output_dir):
        """保存结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON
        with open(output_path / 'bacteria_evaluation.json', 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        # 生成简要报告
        summary = []
        summary.append(f"整体健康评分: {self.results['overall_health']['score']:.1f} - {self.results['overall_health']['grade']}")
        summary.append(f"有益菌状态: {self.results['beneficial_bacteria']['evaluation']}")
        summary.append(f"有害菌风险: {self.results['harmful_bacteria']['risk_level']}")
        
        with open(output_path / 'bacteria_summary.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(summary))
        
        print(f"菌群评估完成，结果保存至: {output_path}")
        for line in summary:
            print(f"  {line}")

def main():
    parser = argparse.ArgumentParser(description='菌群健康评估')
    parser.add_argument('--input', '-i', required=True, help='ASV表路径')
    parser.add_argument('--ranges', '-r', help='正常值范围JSON文件')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    evaluator = BacteriaEvaluator(args.input, args.ranges)
    evaluator.evaluate_beneficial_bacteria()
    evaluator.evaluate_harmful_bacteria()
    evaluator.evaluate_conditional_bacteria()
    evaluator.calculate_overall_health_score()
    evaluator.generate_recommendations()
    evaluator.save_results(args.output)

if __name__ == '__main__':
    main()
