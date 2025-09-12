#!/usr/bin/env python3
"""
疾病风险评估：基于菌群-疾病关联数据库
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path

class DiseaseRiskAssessor:
    def __init__(self, asv_table_path, database_path=None):
        """初始化评估器"""
        self.asv_table = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        self.sample_id = self._get_sample_id()
        self.disease_db = self._load_database(database_path)
        self.results = {}
        
    def _get_sample_id(self):
        """获取样本ID"""
        tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        sample_cols = [col for col in self.asv_table.columns if col not in tax_cols]
        return sample_cols[0] if sample_cols else None
    
    def _load_database(self, path):
        """加载疾病-菌群关联数据库"""
        if path and Path(path).exists():
            # 加载外部数据库（如gutMDisorder）
            return pd.read_csv(path)
        else:
            # 使用内置简化数据库
            return self._get_builtin_database()
    
    def _get_builtin_database(self):
        """内置的疾病-菌群关联数据"""
        data = {
            # IBD（炎症性肠病）
            'IBD': {
                'beneficial': ['Faecalibacterium', 'Akkermansia', 'Coprococcus', 'Roseburia'],
                'harmful': ['Escherichia', 'Streptococcus', 'Eggerthella', 'Klebsiella'],
                'weight': 1.5
            },
            # IBS（肠易激综合征）
            'IBS': {
                'beneficial': ['Bifidobacterium', 'Lactobacillus', 'Faecalibacterium'],
                'harmful': ['Veillonella', 'Streptococcus', 'Ruminococcus_gnavus'],
                'weight': 1.2
            },
            # 糖尿病
            'Diabetes': {
                'beneficial': ['Prevotella', 'Roseburia', 'Akkermansia'],
                'harmful': ['Peptostreptococcus', 'Clostridium', 'Desulfovibrio'],
                'weight': 1.3
            },
            # 肥胖/代谢综合征
            'Obesity': {
                'beneficial': ['Akkermansia', 'Bacteroides', 'Prevotella'],
                'harmful': ['Firmicutes', 'Staphylococcus', 'Enterobacteriaceae'],
                'weight': 1.2
            },
            # 结直肠癌
            'CRC': {
                'beneficial': ['Bifidobacterium', 'Lactobacillus', 'Faecalibacterium'],
                'harmful': ['Fusobacterium', 'Porphyromonas', 'Bacteroides_fragilis', 'Peptostreptococcus'],
                'weight': 1.8
            },
            # 肝病
            'Liver': {
                'beneficial': ['Ruminococcus', 'Akkermansia', 'Faecalibacterium'],
                'harmful': ['Klebsiella', 'Escherichia', 'Veillonella'],
                'weight': 1.4
            },
            # 心血管疾病
            'CVD': {
                'beneficial': ['Prevotella', 'Roseburia', 'Faecalibacterium'],
                'harmful': ['Klebsiella', 'Streptococcus', 'Enterobacter'],
                'weight': 1.3
            },
            # 高血压
            'Hypertension': {
                'beneficial': ['Butyricimonas', 'Akkermansia', 'Faecalibacterium'],
                'harmful': ['Desulfovibrio', 'Klebsiella', 'Streptococcus'],
                'weight': 1.1
            },
            # 抑郁症
            'Depression': {
                'beneficial': ['Bifidobacterium', 'Lactobacillus', 'Faecalibacterium'],
                'harmful': ['Eggerthella', 'Veillonella', 'Alistipes'],
                'weight': 1.2
            },
            # 阿尔茨海默病
            'Alzheimer': {
                'beneficial': ['Bifidobacterium', 'Odoribacter', 'Faecalibacterium'],
                'harmful': ['Escherichia', 'Staphylococcus', 'Bacteroides'],
                'weight': 1.5
            },
            # 便秘
            'Constipation': {
                'beneficial': ['Lactobacillus', 'Bifidobacterium', 'Prevotella'],
                'harmful': ['Clostridium', 'Veillonella', 'Methanobrevibacter'],
                'weight': 1.0
            },
            # 息肉
            'Polyps': {
                'beneficial': ['Bifidobacterium', 'Lactobacillus'],
                'harmful': ['Fusobacterium', 'Peptostreptococcus', 'Dorea'],
                'weight': 1.3
            },
            # 痛风
            'Gout': {
                'beneficial': ['Faecalibacterium', 'Bifidobacterium', 'Roseburia'],
                'harmful': ['Bacteroides', 'Prevotella', 'Escherichia'],
                'weight': 1.1
            },
            # 湿疹
            'Eczema': {
                'beneficial': ['Bifidobacterium', 'Lactobacillus', 'Ruminococcus'],
                'harmful': ['Bacteroides', 'Veillonella', 'Escherichia'],
                'weight': 1.0
            }
        }
        return data
    
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
    
    def assess_disease_risk(self):
        """评估所有疾病风险"""
        disease_risks = {}
        
        for disease, markers in self.disease_db.items():
            # 计算有益菌和有害菌得分
            beneficial_score = 0
            harmful_score = 0
            
            # 有益菌评分
            for bacteria in markers['beneficial']:
                abundance = self._get_bacteria_abundance(bacteria)
                # 丰度越高，保护作用越强
                beneficial_score += min(10, abundance)
            
            # 有害菌评分
            for bacteria in markers['harmful']:
                abundance = self._get_bacteria_abundance(bacteria)
                # 丰度越高，风险越大
                harmful_score += min(10, abundance * 2)  # 有害菌权重更高
            
            # 计算风险分数（0-100）
            # 风险 = 有害菌得分 - 有益菌得分，并标准化
            raw_risk = harmful_score - beneficial_score
            normalized_risk = max(0, min(100, 50 + raw_risk * 2))
            
            # 应用疾病权重
            weighted_risk = normalized_risk * markers.get('weight', 1.0)
            final_risk = min(100, weighted_risk)
            
            # 确定风险等级
            if final_risk < 30:
                risk_level = '低风险'
                color = 'green'
            elif final_risk < 60:
                risk_level = '中风险'
                color = 'yellow'
            else:
                risk_level = '高风险'
                color = 'red'
            
            disease_risks[disease] = {
                'risk_score': round(final_risk, 1),
                'risk_level': risk_level,
                'color': color,
                'beneficial_score': round(beneficial_score, 2),
                'harmful_score': round(harmful_score, 2),
                'key_findings': self._get_key_findings(disease, markers)
            }
        
        self.results['disease_risks'] = disease_risks
        
        # 计算整体健康风险
        avg_risk = np.mean([d['risk_score'] for d in disease_risks.values()])
        high_risk_count = sum(1 for d in disease_risks.values() if d['risk_level'] == '高风险')
        
        self.results['overall_risk'] = {
            'average_risk_score': round(avg_risk, 1),
            'high_risk_diseases': high_risk_count,
            'health_status': self._get_health_status(avg_risk)
        }
    
    def _get_key_findings(self, disease, markers):
        """获取关键发现"""
        findings = []
        
        # 检查关键有害菌
        for bacteria in markers['harmful'][:3]:  # 前3个最重要的
            abundance = self._get_bacteria_abundance(bacteria)
            if abundance > 1.0:  # 超过1%认为较高
                findings.append(f"{bacteria}偏高 ({abundance:.2f}%)")
        
        # 检查关键有益菌
        for bacteria in markers['beneficial'][:3]:
            abundance = self._get_bacteria_abundance(bacteria)
            if abundance < 0.1:  # 低于0.1%认为较低
                findings.append(f"{bacteria}偏低 ({abundance:.2f}%)")
        
        return findings[:3]  # 最多返回3个发现
    
    def _get_health_status(self, avg_risk):
        """获取健康状态"""
        if avg_risk < 30:
            return '健康'
        elif avg_risk < 50:
            return '亚健康'
        elif avg_risk < 70:
            return '需关注'
        else:
            return '需干预'
    
    def generate_prevention_advice(self):
        """生成预防建议"""
        advice = {}
        
        for disease, risk_data in self.results['disease_risks'].items():
            if risk_data['risk_level'] in ['中风险', '高风险']:
                disease_advice = self._get_disease_specific_advice(disease)
                if disease_advice:
                    advice[disease] = disease_advice
        
        self.results['prevention_advice'] = advice
    
    def _get_disease_specific_advice(self, disease):
        """获取特定疾病的预防建议"""
        advice_db = {
            'IBD': {
                'diet': '增加膳食纤维，避免辛辣刺激食物',
                'lifestyle': '规律作息，适度运动，管理压力',
                'supplements': '益生菌、益生元、Omega-3脂肪酸'
            },
            'Diabetes': {
                'diet': '低糖低脂饮食，增加全谷物摄入',
                'lifestyle': '规律运动，控制体重',
                'supplements': '膳食纤维、益生菌'
            },
            'CRC': {
                'diet': '增加蔬果摄入，减少红肉和加工肉类',
                'lifestyle': '定期筛查，戒烟限酒',
                'supplements': '膳食纤维、益生菌、维生素D'
            },
            'Depression': {
                'diet': '地中海饮食，富含Omega-3的食物',
                'lifestyle': '规律运动，充足睡眠，社交活动',
                'supplements': '益生菌、B族维生素、Omega-3'
            },
            'CVD': {
                'diet': '低盐低脂饮食，增加蔬果和全谷物',
                'lifestyle': '规律运动，控制血压和血脂',
                'supplements': '膳食纤维、Omega-3、辅酶Q10'
            }
        }
        
        return advice_db.get(disease, {
            'diet': '均衡饮食，增加膳食纤维',
            'lifestyle': '规律作息，适度运动',
            'supplements': '益生菌、益生元'
        })
    
    def save_results(self, output_dir):
        """保存结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON
        with open(output_path / 'disease_risk_assessment.json', 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        # 生成风险报告
        high_risk_diseases = [d for d, r in self.results['disease_risks'].items() 
                             if r['risk_level'] == '高风险']
        
        report = []
        report.append(f"整体健康状态: {self.results['overall_risk']['health_status']}")
        report.append(f"平均风险分: {self.results['overall_risk']['average_risk_score']:.1f}")
        
        if high_risk_diseases:
            report.append(f"\n高风险疾病 ({len(high_risk_diseases)}):")
            for disease in high_risk_diseases:
                risk_data = self.results['disease_risks'][disease]
                report.append(f"  - {disease}: {risk_data['risk_score']:.1f}分")
                if risk_data['key_findings']:
                    report.append(f"    关键发现: {', '.join(risk_data['key_findings'])}")
        
        with open(output_path / 'risk_report.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(report))
        
        print(f"疾病风险评估完成，结果保存至: {output_path}")
        for line in report[:5]:  # 打印前5行
            print(f"  {line}")

def main():
    parser = argparse.ArgumentParser(description='疾病风险评估')
    parser.add_argument('--input', '-i', required=True, help='ASV表路径')
    parser.add_argument('--database', '-d', help='疾病-菌群关联数据库')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    assessor = DiseaseRiskAssessor(args.input, args.database)
    assessor.assess_disease_risk()
    assessor.generate_prevention_advice()
    assessor.save_results(args.output)

if __name__ == '__main__':
    main()
