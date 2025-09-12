#!/usr/bin/env python3
"""
生物年龄预测：基于肠道菌群组成
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path

class AgePredictor:
    def __init__(self, asv_table_path, markers_path=None):
        """初始化预测器"""
        self.asv_table = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        self.sample_id = self._get_sample_id()
        self.age_markers = self._load_markers(markers_path)
        self.results = {}
        
    def _get_sample_id(self):
        """获取样本ID"""
        tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        sample_cols = [col for col in self.asv_table.columns if col not in tax_cols]
        return sample_cols[0] if sample_cols else None
    
    def _load_markers(self, path):
        """加载年龄相关标记菌"""
        if path and Path(path).exists():
            with open(path, 'r', encoding='utf-8') as f:
                return json.load(f)
        else:
            return self._get_default_markers()
    
    def _get_default_markers(self):
        """默认年龄标记菌（基于文献）"""
        return {
            'youth_associated': {
                # 年轻相关菌（负相关）
                'Bifidobacterium': {'weight': -2.5, 'optimal_range': [5, 15]},
                'Lactobacillus': {'weight': -1.8, 'optimal_range': [0.5, 3]},
                'Prevotella': {'weight': -1.5, 'optimal_range': [10, 40]},
                'Faecalibacterium': {'weight': -2.0, 'optimal_range': [5, 15]},
                'Akkermansia': {'weight': -1.2, 'optimal_range': [1, 5]}
            },
            'aging_associated': {
                # 衰老相关菌（正相关）
                'Escherichia': {'weight': 2.0, 'threshold': 1.0},
                'Enterococcus': {'weight': 1.8, 'threshold': 0.5},
                'Streptococcus': {'weight': 1.5, 'threshold': 2.0},
                'Clostridium': {'weight': 1.3, 'threshold': 3.0},
                'Staphylococcus': {'weight': 1.6, 'threshold': 0.1},
                'Klebsiella': {'weight': 1.4, 'threshold': 0.5}
            },
            'baseline_age': 40  # 基准年龄
        }
    
    def _get_bacteria_abundance(self, bacteria_name, level='Genus'):
        """获取细菌相对丰度"""
        if level not in self.asv_table.columns:
            return 0
            
        total_abundance = 0
        total_reads = self.asv_table[self.sample_id].sum()
        
        for idx, row in self.asv_table.iterrows():
            taxon = str(row[level]) if pd.notna(row[level]) else ''
            taxon = taxon.replace(f'{level[0].lower()}__', '').strip()
            
            if bacteria_name.lower() in taxon.lower() or taxon.lower() in bacteria_name.lower():
                total_abundance += row[self.sample_id]
        
        return (total_abundance / total_reads * 100) if total_reads > 0 else 0
    
    def predict_biological_age(self):
        """预测生物年龄"""
        baseline_age = self.age_markers['baseline_age']
        age_adjustment = 0
        marker_details = {}
        
        # 计算年轻相关菌的影响
        youth_score = 0
        for bacteria, params in self.age_markers['youth_associated'].items():
            abundance = self._get_bacteria_abundance(bacteria)
            optimal_min, optimal_max = params['optimal_range']
            
            # 在最佳范围内得分最高
            if optimal_min <= abundance <= optimal_max:
                score = 1.0
                status = '最佳'
            elif abundance < optimal_min:
                score = abundance / optimal_min if optimal_min > 0 else 0
                status = '偏低'
            else:
                score = max(0, 1 - (abundance - optimal_max) / optimal_max)
                status = '偏高'
            
            youth_score += score * abs(params['weight'])
            age_adjustment += score * params['weight']
            
            marker_details[bacteria] = {
                'abundance': round(abundance, 3),
                'optimal_range': params['optimal_range'],
                'status': status,
                'contribution': round(score * params['weight'], 2)
            }
        
        # 计算衰老相关菌的影响
        aging_score = 0
        for bacteria, params in self.age_markers['aging_associated'].items():
            abundance = self._get_bacteria_abundance(bacteria)
            threshold = params['threshold']
            
            # 超过阈值增加年龄
            if abundance > threshold:
                score = min(2, abundance / threshold)  # 最多2倍影响
                status = '偏高'
            else:
                score = 0
                status = '正常'
            
            aging_score += score * params['weight']
            age_adjustment += score * params['weight']
            
            marker_details[bacteria] = {
                'abundance': round(abundance, 3),
                'threshold': threshold,
                'status': status,
                'contribution': round(score * params['weight'], 2)
            }
        
        # 计算最终生物年龄
        biological_age = baseline_age + age_adjustment
        biological_age = max(20, min(90, biological_age))  # 限制在20-90岁范围
        
        self.results['age_prediction'] = {
            'biological_age': round(biological_age, 1),
            'baseline_age': baseline_age,
            'age_adjustment': round(age_adjustment, 1),
            'youth_score': round(youth_score, 2),
            'aging_score': round(aging_score, 2),
            'marker_details': marker_details
        }
    
    def analyze_age_status(self, chronological_age=None):
        """分析年龄状态"""
        bio_age = self.results['age_prediction']['biological_age']
        
        if chronological_age:
            age_diff = bio_age - chronological_age
            
            if age_diff < -5:
                status = '年轻态'
                description = f'肠道年龄比实际年龄年轻{abs(age_diff):.1f}岁'
            elif age_diff > 5:
                status = '衰老态'
                description = f'肠道年龄比实际年龄老{age_diff:.1f}岁'
            else:
                status = '同龄态'
                description = '肠道年龄与实际年龄相符'
            
            self.results['age_status'] = {
                'chronological_age': chronological_age,
                'biological_age': bio_age,
                'age_difference': round(age_diff, 1),
                'status': status,
                'description': description
            }
        else:
            # 没有实际年龄时，基于生物年龄给出评估
            if bio_age < 35:
                status = '年轻型'
            elif bio_age < 50:
                status = '中年型'
            elif bio_age < 65:
                status = '成熟型'
            else:
                status = '老年型'
            
            self.results['age_status'] = {
                'biological_age': bio_age,
                'status': status,
                'description': f'肠道微生物组呈现{status}特征'
            }
    
    def generate_rejuvenation_advice(self):
        """生成抗衰老建议"""
        advice = {
            'dietary': [],
            'lifestyle': [],
            'supplements': []
        }
        
        # 基于标记菌状态生成建议
        marker_details = self.results['age_prediction']['marker_details']
        
        # 检查年轻相关菌
        low_youth_bacteria = []
        for bacteria in self.age_markers['youth_associated']:
            if bacteria in marker_details and marker_details[bacteria]['status'] == '偏低':
                low_youth_bacteria.append(bacteria)
        
        if low_youth_bacteria:
            if 'Bifidobacterium' in low_youth_bacteria:
                advice['dietary'].append('增加酸奶、发酵乳制品摄入')
                advice['supplements'].append('补充双歧杆菌益生菌')
            
            if 'Lactobacillus' in low_youth_bacteria:
                advice['dietary'].append('增加泡菜、酸菜等发酵蔬菜')
                advice['supplements'].append('补充乳酸菌益生菌')
            
            if 'Prevotella' in low_youth_bacteria:
                advice['dietary'].append('增加全谷物、豆类摄入')
            
            if 'Faecalibacterium' in low_youth_bacteria:
                advice['dietary'].append('增加膳食纤维、抗性淀粉')
                advice['supplements'].append('补充益生元（菊粉、低聚果糖）')
        
        # 检查衰老相关菌
        high_aging_bacteria = []
        for bacteria in self.age_markers['aging_associated']:
            if bacteria in marker_details and marker_details[bacteria]['status'] == '偏高':
                high_aging_bacteria.append(bacteria)
        
        if high_aging_bacteria:
            advice['dietary'].append('减少高脂高糖食物')
            advice['dietary'].append('增加多酚类食物（绿茶、蓝莓、黑巧克力）')
            advice['lifestyle'].append('规律运动，每周至少150分钟中等强度运动')
            advice['lifestyle'].append('保证充足睡眠（7-8小时）')
        
        # 通用抗衰老建议
        advice['lifestyle'].append('管理压力，练习冥想或瑜伽')
        advice['lifestyle'].append('保持社交活动')
        advice['dietary'].append('地中海饮食模式')
        
        self.results['rejuvenation_advice'] = advice
    
    def calculate_aging_rate(self):
        """计算衰老速度评分"""
        youth_score = self.results['age_prediction']['youth_score']
        aging_score = self.results['age_prediction']['aging_score']
        
        # 衰老速度 = 衰老菌得分 / (年轻菌得分 + 1)
        aging_rate = aging_score / (youth_score + 1)
        
        if aging_rate < 0.5:
            rate_level = '缓慢'
            description = '肠道微生物衰老速度较慢，保持良好'
        elif aging_rate < 1.0:
            rate_level = '正常'
            description = '肠道微生物衰老速度正常'
        elif aging_rate < 2.0:
            rate_level = '较快'
            description = '肠道微生物衰老速度较快，需要关注'
        else:
            rate_level = '快速'
            description = '肠道微生物衰老速度过快，建议干预'
        
        self.results['aging_rate'] = {
            'rate_score': round(aging_rate, 2),
            'rate_level': rate_level,
            'description': description
        }
    
    def save_results(self, output_dir):
        """保存结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON
        with open(output_path / 'age_prediction.json', 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        # 生成报告
        report = []
        report.append(f"生物年龄预测: {self.results['age_prediction']['biological_age']:.1f}岁")
        
        if 'age_status' in self.results:
            report.append(f"年龄状态: {self.results['age_status']['status']}")
            report.append(f"说明: {self.results['age_status']['description']}")
        
        if 'aging_rate' in self.results:
            report.append(f"\n衰老速度: {self.results['aging_rate']['rate_level']}")
            report.append(f"评分: {self.results['aging_rate']['rate_score']:.2f}")
        
        # 关键标记菌
        report.append("\n关键年龄标记菌:")
        marker_details = self.results['age_prediction']['marker_details']
        
        # 找出影响最大的菌
        sorted_markers = sorted(marker_details.items(), 
                              key=lambda x: abs(x[1]['contribution']), 
                              reverse=True)[:5]
        
        for bacteria, details in sorted_markers:
            contribution = details['contribution']
            if contribution < 0:
                effect = '减龄'
            else:
                effect = '增龄'
            report.append(f"  - {bacteria}: {details['abundance']:.2f}% ({effect} {abs(contribution):.1f}岁)")
        
        with open(output_path / 'age_report.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"年龄预测完成，结果保存至: {output_path}")
        for line in report[:5]:
            print(f"  {line}")

def main():
    parser = argparse.ArgumentParser(description='生物年龄预测')
    parser.add_argument('--input', '-i', required=True, help='ASV表路径')
    parser.add_argument('--markers', '-m', help='年龄标记菌JSON文件')
    parser.add_argument('--age', '-a', type=int, help='实际年龄（可选）')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    predictor = AgePredictor(args.input, args.markers)
    predictor.predict_biological_age()
    predictor.analyze_age_status(args.age)
    predictor.calculate_aging_rate()
    predictor.generate_rejuvenation_advice()
    predictor.save_results(args.output)

if __name__ == '__main__':
    main()
