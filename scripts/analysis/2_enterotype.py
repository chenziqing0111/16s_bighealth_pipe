#!/usr/bin/env python3
"""
肠型分析：基于属水平组成进行聚类
"""

import pandas as pd
import numpy as np
import json
import argparse
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

class EnterotypeAnalyzer:
    def __init__(self, asv_table_path):
        """初始化分析器"""
        self.asv_table = pd.read_csv(asv_table_path, sep='\t', index_col=0)
        self.sample_id = self._get_sample_id()
        self.results = {}
        
    def _get_sample_id(self):
        """获取样本ID"""
        tax_cols = ['Taxon', 'Confidence', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        sample_cols = [col for col in self.asv_table.columns if col not in tax_cols]
        return sample_cols[0] if sample_cols else None
    
    def prepare_genus_data(self):
        """准备属水平数据"""
        if 'Genus' not in self.asv_table.columns:
            self.results['error'] = '缺少属水平分类信息'
            return None
            
        # 汇总到属水平
        genus_abundance = {}
        total_reads = 0
        
        for idx, row in self.asv_table.iterrows():
            genus = row['Genus']
            if pd.notna(genus) and genus != 'Unclassified':
                genus = genus.replace('g__', '').strip()
                if genus:
                    if genus not in genus_abundance:
                        genus_abundance[genus] = 0
                    genus_abundance[genus] += row[self.sample_id]
                    total_reads += row[self.sample_id]
        
        # 转换为相对丰度
        if total_reads > 0:
            for genus in genus_abundance:
                genus_abundance[genus] = (genus_abundance[genus] / total_reads) * 100
        
        return genus_abundance
    
    def determine_enterotype(self):
        """确定肠型"""
        genus_data = self.prepare_genus_data()
        
        if not genus_data:
            return
        
        # 关键属的丰度
        bacteroides = genus_data.get('Bacteroides', 0)
        prevotella = genus_data.get('Prevotella', 0) + genus_data.get('Prevotella_9', 0)
        ruminococcus = genus_data.get('Ruminococcus', 0) + genus_data.get('Ruminococcus_1', 0) + genus_data.get('Ruminococcus_2', 0)
        
        # 简单判定规则
        dominant_genera = {
            'Bacteroides': bacteroides,
            'Prevotella': prevotella,
            'Ruminococcus': ruminococcus
        }
        
        # 找出优势属
        dominant = max(dominant_genera, key=dominant_genera.get)
        
        # 确定肠型
        if dominant == 'Bacteroides':
            enterotype = '拟杆菌型'
            description = '以拟杆菌属为主导，常见于高蛋白、高脂肪饮食人群'
        elif dominant == 'Prevotella':
            enterotype = '普氏菌型'
            description = '以普氏菌属为主导，常见于高纤维、植物性饮食人群'
        else:
            enterotype = '瘤胃球菌型'
            description = '以瘤胃球菌属为主导，混合型饮食模式'
        
        self.results['enterotype'] = {
            'type': enterotype,
            'description': description,
            'dominant_genus': dominant,
            'key_genera_abundance': {
                'Bacteroides': round(bacteroides, 2),
                'Prevotella': round(prevotella, 2),
                'Ruminococcus': round(ruminococcus, 2)
            },
            'confidence': self._calculate_confidence(dominant_genera)
        }
        
        # 保存所有属的丰度（用于详细分析）
        sorted_genera = sorted(genus_data.items(), key=lambda x: x[1], reverse=True)[:20]
        self.results['genus_profile'] = {
            'genera': [g[0] for g in sorted_genera],
            'abundance': [round(g[1], 3) for g in sorted_genera]
        }
    
    def _calculate_confidence(self, dominant_genera):
        """计算肠型判定的置信度"""
        values = list(dominant_genera.values())
        max_val = max(values)
        total = sum(values)
        
        if total == 0:
            return 'low'
        
        # 如果最高丰度占比超过50%，高置信度
        if max_val / total > 0.5:
            return 'high'
        elif max_val / total > 0.33:
            return 'medium'
        else:
            return 'low'
    
    def analyze_enterotype_features(self):
        """分析肠型相关特征"""
        if 'enterotype' not in self.results:
            return
            
        enterotype = self.results['enterotype']['type']
        
        # 不同肠型的特征
        features = {
            '拟杆菌型': {
                'diet_recommendation': '增加膳食纤维摄入，减少动物性脂肪',
                'health_implications': '可能与肥胖风险增加相关，建议控制热量摄入',
                'beneficial_foods': ['全谷物', '蔬菜', '水果', '豆类'],
                'avoid_foods': ['红肉', '加工肉类', '高脂乳制品']
            },
            '普氏菌型': {
                'diet_recommendation': '保持高纤维饮食，适量增加优质蛋白',
                'health_implications': '通常与较好的血糖控制相关',
                'beneficial_foods': ['糙米', '燕麦', '绿叶蔬菜', '坚果'],
                'avoid_foods': ['精制糖', '加工食品']
            },
            '瘤胃球菌型': {
                'diet_recommendation': '均衡饮食，注意食物多样性',
                'health_implications': '中间型，需要根据其他指标综合评估',
                'beneficial_foods': ['发酵食品', '益生元食物', '鱼类'],
                'avoid_foods': ['过度加工食品', '含糖饮料']
            }
        }
        
        if enterotype in features:
            self.results['enterotype_features'] = features[enterotype]
    
    def save_results(self, output_dir):
        """保存结果"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON
        with open(output_path / 'enterotype_analysis.json', 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        print(f"肠型分析完成，结果保存至: {output_path}")
        
        if 'enterotype' in self.results:
            print(f"  肠型: {self.results['enterotype']['type']}")
            print(f"  优势属: {self.results['enterotype']['dominant_genus']}")
            print(f"  置信度: {self.results['enterotype']['confidence']}")

def main():
    parser = argparse.ArgumentParser(description='肠型分析')
    parser.add_argument('--input', '-i', required=True, help='ASV表路径')
    parser.add_argument('--output', '-o', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    analyzer = EnterotypeAnalyzer(args.input)
    analyzer.determine_enterotype()
    analyzer.analyze_enterotype_features()
    analyzer.save_results(args.output)

if __name__ == '__main__':
    main()
