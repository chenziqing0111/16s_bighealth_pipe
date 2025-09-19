#!/usr/bin/env python3
"""
中文注释模块
为分析结果添加中文注释，包括细菌、代谢通路、酶等
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import logging

class ChineseAnnotator:
    def __init__(self, database_dir='database'):
        """
        初始化中文注释器
        
        Args:
            database_dir: 注释数据库目录
        """
        self.database_dir = Path(database_dir)
        
        # 设置日志（必须在使用前初始化）
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
        self.logger = logging.getLogger(__name__)
        
        # 加载注释数据
        self.annotations = self._load_annotations()
    
    def _load_annotations(self):
        """加载所有注释数据"""
        annotations = {}
        
        # 加载细菌注释
        bacteria_file = self.database_dir / 'core_bacteria_annotations.json'
        if bacteria_file.exists():
            with open(bacteria_file, 'r', encoding='utf-8') as f:
                annotations['bacteria'] = json.load(f)
            self.logger.info(f"加载了 {len(annotations['bacteria'])} 个细菌注释")
        else:
            self.logger.warning(f"未找到细菌注释文件: {bacteria_file}")
            annotations['bacteria'] = {}
        
        # 加载通路注释
        pathway_file = self.database_dir / 'core_pathway_translations.json'
        if pathway_file.exists():
            with open(pathway_file, 'r', encoding='utf-8') as f:
                annotations['pathways'] = json.load(f)
            self.logger.info(f"加载了 {len(annotations['pathways'])} 个通路注释")
        else:
            self.logger.warning(f"未找到通路注释文件: {pathway_file}")
            annotations['pathways'] = {}
        
        # 加载EC注释
        ec_file = self.database_dir / 'core_ec_translations.json'
        if ec_file.exists():
            with open(ec_file, 'r', encoding='utf-8') as f:
                annotations['ecs'] = json.load(f)
            self.logger.info(f"加载了 {len(annotations['ecs'])} 个EC注释")
        else:
            self.logger.warning(f"未找到EC注释文件: {ec_file}")
            annotations['ecs'] = {}
        
        return annotations
    
    def annotate_bacteria(self, bacteria_data):
        """
        为细菌数据添加中文注释
        
        Args:
            bacteria_data: 细菌评估结果（来自3_bacteria_eval.py）
        
        Returns:
            添加了中文注释的数据
        """
        annotated_data = bacteria_data.copy() if isinstance(bacteria_data, dict) else bacteria_data
        
        # 处理有益菌
        if 'beneficial_bacteria' in annotated_data and 'bacteria' in annotated_data['beneficial_bacteria']:
            for bacteria_name in annotated_data['beneficial_bacteria']['bacteria']:
                if bacteria_name in self.annotations['bacteria']:
                    annotation = self.annotations['bacteria'][bacteria_name]
                    annotated_data['beneficial_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': annotation.get('cn_name', bacteria_name),
                        'description': annotation.get('description', ''),
                        'functions': annotation.get('functions', []),
                        'health_impact': annotation.get('health_impact', ''),
                        'food_sources': annotation.get('food_sources', [])
                    })
                else:
                    # 未找到注释的细菌，提供默认值
                    annotated_data['beneficial_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': bacteria_name,
                        'description': '益生菌',
                        'functions': [],
                        'health_impact': '',
                        'food_sources': []
                    })
        
        # 处理有害菌
        if 'harmful_bacteria' in annotated_data and 'bacteria' in annotated_data['harmful_bacteria']:
            for bacteria_name in annotated_data['harmful_bacteria']['bacteria']:
                if bacteria_name in self.annotations['bacteria']:
                    annotation = self.annotations['bacteria'][bacteria_name]
                    annotated_data['harmful_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': annotation.get('cn_name', bacteria_name),
                        'description': annotation.get('description', ''),
                        'functions': annotation.get('functions', []),
                        'health_impact': annotation.get('health_impact', ''),
                        'risk_factors': annotation.get('risk_factors', [])
                    })
                else:
                    annotated_data['harmful_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': bacteria_name,
                        'description': '潜在致病菌',
                        'functions': [],
                        'health_impact': '',
                        'risk_factors': []
                    })
        
        # 处理条件致病菌
        if 'conditional_bacteria' in annotated_data and 'bacteria' in annotated_data['conditional_bacteria']:
            for bacteria_name in annotated_data['conditional_bacteria']['bacteria']:
                if bacteria_name in self.annotations['bacteria']:
                    annotation = self.annotations['bacteria'][bacteria_name]
                    annotated_data['conditional_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': annotation.get('cn_name', bacteria_name),
                        'description': annotation.get('description', ''),
                        'category': annotation.get('category', 'conditional')
                    })
                else:
                    annotated_data['conditional_bacteria']['bacteria'][bacteria_name].update({
                        'cn_name': bacteria_name,
                        'description': '条件致病菌',
                        'category': 'conditional'
                    })
        
        return annotated_data
    
    def annotate_pathways(self, pathway_data):
        """
        为代谢通路数据添加中文注释
        
        Args:
            pathway_data: 通路数据（来自功能预测）
        
        Returns:
            添加了中文注释的数据
        """
        annotated_pathways = {}
        
        if isinstance(pathway_data, pd.DataFrame):
            # 处理DataFrame格式
            for pathway_id in pathway_data.index:
                if pathway_id in self.annotations['pathways']:
                    annotation = self.annotations['pathways'][pathway_id]
                    annotated_pathways[pathway_id] = {
                        'id': pathway_id,
                        'cn_name': annotation.get('cn_name', pathway_id),
                        'category': annotation.get('category', ''),
                        'description': annotation.get('description', ''),
                        'importance': annotation.get('importance', 'medium'),
                        'related_nutrients': annotation.get('related_nutrients', []),
                        'health_impact': annotation.get('health_impact', ''),
                        'abundance': float(pathway_data.loc[pathway_id].iloc[0]) if len(pathway_data.columns) > 0 else 0
                    }
                else:
                    # 尝试自动翻译
                    cn_name = self._auto_translate_pathway(pathway_id)
                    annotated_pathways[pathway_id] = {
                        'id': pathway_id,
                        'cn_name': cn_name,
                        'category': '其他',
                        'description': '',
                        'importance': 'low',
                        'abundance': float(pathway_data.loc[pathway_id].iloc[0]) if len(pathway_data.columns) > 0 else 0
                    }
        
        elif isinstance(pathway_data, dict):
            # 处理字典格式
            for pathway_id, abundance in pathway_data.items():
                if pathway_id in self.annotations['pathways']:
                    annotation = self.annotations['pathways'][pathway_id]
                    annotated_pathways[pathway_id] = {
                        'id': pathway_id,
                        'cn_name': annotation.get('cn_name', pathway_id),
                        'category': annotation.get('category', ''),
                        'description': annotation.get('description', ''),
                        'importance': annotation.get('importance', 'medium'),
                        'related_nutrients': annotation.get('related_nutrients', []),
                        'health_impact': annotation.get('health_impact', ''),
                        'abundance': float(abundance) if abundance else 0
                    }
                else:
                    cn_name = self._auto_translate_pathway(pathway_id)
                    annotated_pathways[pathway_id] = {
                        'id': pathway_id,
                        'cn_name': cn_name,
                        'category': '其他',
                        'description': '',
                        'importance': 'low',
                        'abundance': float(abundance) if abundance else 0
                    }
        
        return annotated_pathways
    
    def annotate_ecs(self, ec_data):
        """
        为EC酶数据添加中文注释
        
        Args:
            ec_data: EC数据
        
        Returns:
            添加了中文注释的数据
        """
        annotated_ecs = {}
        
        if isinstance(ec_data, pd.DataFrame):
            # 处理DataFrame格式
            for ec_id in ec_data.index:
                if ec_id in self.annotations['ecs']:
                    annotation = self.annotations['ecs'][ec_id]
                    annotated_ecs[ec_id] = {
                        'id': ec_id,
                        'cn_name': annotation.get('cn_name', ec_id),
                        'category': annotation.get('category', ''),
                        'function': annotation.get('function', ''),
                        'importance': annotation.get('importance', 'medium'),
                        'health_relevance': annotation.get('health_relevance', ''),
                        'abundance': float(ec_data.loc[ec_id].iloc[0]) if len(ec_data.columns) > 0 else 0
                    }
                else:
                    # 基于EC分类自动生成描述
                    category = self._get_ec_category(ec_id)
                    annotated_ecs[ec_id] = {
                        'id': ec_id,
                        'cn_name': ec_id,
                        'category': category,
                        'function': '',
                        'importance': 'low',
                        'abundance': float(ec_data.loc[ec_id].iloc[0]) if len(ec_data.columns) > 0 else 0
                    }
        
        elif isinstance(ec_data, dict):
            # 处理字典格式
            for ec_id, abundance in ec_data.items():
                if ec_id in self.annotations['ecs']:
                    annotation = self.annotations['ecs'][ec_id]
                    annotated_ecs[ec_id] = {
                        'id': ec_id,
                        'cn_name': annotation.get('cn_name', ec_id),
                        'category': annotation.get('category', ''),
                        'function': annotation.get('function', ''),
                        'importance': annotation.get('importance', 'medium'),
                        'health_relevance': annotation.get('health_relevance', ''),
                        'abundance': float(abundance) if abundance else 0
                    }
                else:
                    category = self._get_ec_category(ec_id)
                    annotated_ecs[ec_id] = {
                        'id': ec_id,
                        'cn_name': ec_id,
                        'category': category,
                        'function': '',
                        'importance': 'low',
                        'abundance': float(abundance) if abundance else 0
                    }
        
        return annotated_ecs
    
    def _auto_translate_pathway(self, pathway_id):
        """自动翻译通路名称"""
        translations = {
            'PWY': '通路',
            'SYN': '合成',
            'BIOSYNTHESIS': '生物合成',
            'DEGRADATION': '降解',
            'FERMENTATION': '发酵',
            'METABOLISM': '代谢',
            'CAT': '分解代谢',
            'ANAERO': '厌氧',
            'GLYCOLYSIS': '糖酵解',
            'TCA': '三羧酸循环',
            'OXIDO': '氧化',
            'REDUCTION': '还原'
        }
        
        result = pathway_id
        for eng, chn in translations.items():
            if eng in pathway_id.upper():
                result = result.replace(eng, chn)
        
        return result
    
    def _get_ec_category(self, ec_id):
        """根据EC编号获取酶类别"""
        if not ec_id.startswith('EC:'):
            return '未分类'
        
        try:
            first_digit = ec_id.split('.')[0].replace('EC:', '')
            categories = {
                '1': '氧化还原酶',
                '2': '转移酶',
                '3': '水解酶',
                '4': '裂解酶',
                '5': '异构酶',
                '6': '连接酶'
            }
            return categories.get(first_digit, '未分类')
        except:
            return '未分类'
    
    def annotate_nutrition_metabolism(self, functional_data):
        """
        为营养代谢数据添加中文注释
        
        Args:
            functional_data: 功能预测数据
        
        Returns:
            注释后的营养代谢数据
        """
        nutrition_annotations = {
            'vitamins': {
                'B1': {'cn_name': '维生素B1（硫胺素）', 'function': '能量代谢，神经功能'},
                'B2': {'cn_name': '维生素B2（核黄素）', 'function': '能量代谢，抗氧化'},
                'B3': {'cn_name': '维生素B3（烟酸）', 'function': '能量代谢，DNA修复'},
                'B5': {'cn_name': '维生素B5（泛酸）', 'function': '脂肪代谢，激素合成'},
                'B6': {'cn_name': '维生素B6（吡哆醇）', 'function': '蛋白质代谢，神经递质'},
                'B7': {'cn_name': '维生素B7（生物素）', 'function': '脂肪酸合成，糖异生'},
                'B9': {'cn_name': '维生素B9（叶酸）', 'function': 'DNA合成，红细胞生成'},
                'B12': {'cn_name': '维生素B12（钴胺素）', 'function': '红细胞生成，神经功能'},
                'K': {'cn_name': '维生素K', 'function': '凝血，骨代谢'}
            },
            'scfa': {
                'acetate': {'cn_name': '乙酸', 'function': '能量供应，调节血脂'},
                'propionate': {'cn_name': '丙酸', 'function': '调节血糖，抑制胆固醇'},
                'butyrate': {'cn_name': '丁酸', 'function': '肠道健康，抗炎'}
            },
            'amino_acids': {
                'arginine': {'cn_name': '精氨酸', 'function': '免疫调节，伤口愈合'},
                'glutamine': {'cn_name': '谷氨酰胺', 'function': '肠道修复，免疫支持'},
                'tryptophan': {'cn_name': '色氨酸', 'function': '5-羟色胺前体，睡眠调节'}
            },
            'pigments': {
                'carotenoid': {'cn_name': '类胡萝卜素', 'function': '抗氧化，维生素A前体'},
                'heme': {'cn_name': '血红素', 'function': '氧气运输，电子传递'},
                'melanin': {'cn_name': '黑色素', 'function': '抗氧化，UV防护'}
            }
        }
        
        annotated_data = {}
        
        # 处理维生素合成
        if 'vitamin_synthesis' in functional_data:
            annotated_data['vitamin_synthesis'] = {}
            for vitamin, data in functional_data['vitamin_synthesis'].items():
                if vitamin in nutrition_annotations['vitamins']:
                    annotated_data['vitamin_synthesis'][vitamin] = {
                        **data,
                        **nutrition_annotations['vitamins'][vitamin]
                    }
                else:
                    annotated_data['vitamin_synthesis'][vitamin] = data
        
        # 处理短链脂肪酸
        if 'scfa_production' in functional_data:
            annotated_data['scfa_production'] = {}
            for scfa, data in functional_data['scfa_production'].items():
                scfa_key = scfa.lower().replace('丁酸', 'butyrate').replace('丙酸', 'propionate').replace('乙酸', 'acetate')
                if scfa_key in nutrition_annotations['scfa']:
                    annotated_data['scfa_production'][scfa] = {
                        **data,
                        **nutrition_annotations['scfa'][scfa_key]
                    }
                else:
                    annotated_data['scfa_production'][scfa] = data
        
        return annotated_data
    
    def process_sample_analysis(self, sample_dir):
        """
        处理单个样本的所有分析结果，添加中文注释
        
        Args:
            sample_dir: 样本分析结果目录
        
        Returns:
            所有注释后的数据
        """
        sample_path = Path(sample_dir)
        annotated_results = {}
        
        # 1. 处理细菌评估结果
        bacteria_file = sample_path / 'bacteria_scores' / 'bacteria_evaluation.json'
        if bacteria_file.exists():
            with open(bacteria_file, 'r', encoding='utf-8') as f:
                bacteria_data = json.load(f)
            annotated_results['bacteria'] = self.annotate_bacteria(bacteria_data)
            self.logger.info("✓ 细菌评估数据已注释")
        
        # 2. 处理功能预测结果
        functional_file = sample_path / 'functional' / 'functional_prediction.json'
        if functional_file.exists():
            with open(functional_file, 'r', encoding='utf-8') as f:
                functional_data = json.load(f)
            
            # 注释营养代谢
            if 'vitamin_synthesis' in functional_data or 'scfa_production' in functional_data:
                annotated_results['nutrition'] = self.annotate_nutrition_metabolism(functional_data)
                self.logger.info("✓ 营养代谢数据已注释")
        
        # 3. 处理通路数据
        pathway_file = sample_path.parent.parent / 'preprocessing' / 'functional_pathway_annotation.tsv'
        if pathway_file.exists():
            pathway_df = pd.read_csv(pathway_file, sep='\t', index_col=0)
            # 只取该样本的数据
            sample_name = sample_path.name
            if sample_name in pathway_df.columns:
                sample_pathways = pathway_df[[sample_name]]
                # 取TOP 20通路
                top_pathways = sample_pathways.nlargest(20, sample_name)
                annotated_results['pathways'] = self.annotate_pathways(top_pathways)
                self.logger.info(f"✓ TOP 20 代谢通路已注释")
        
        # 4. 处理EC数据
        ec_file = sample_path.parent.parent / 'preprocessing' / 'functional_ec_annotation.tsv'
        if ec_file.exists():
            ec_df = pd.read_csv(ec_file, sep='\t', index_col=0)
            sample_name = sample_path.name
            if sample_name in ec_df.columns:
                sample_ecs = ec_df[[sample_name]]
                # 取TOP 20酶
                top_ecs = sample_ecs.nlargest(20, sample_name)
                annotated_results['ecs'] = self.annotate_ecs(top_ecs)
                self.logger.info(f"✓ TOP 20 酶已注释")
        
        # 5. 保存综合注释结果
        output_file = sample_path / 'cn_annotations.json'
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(annotated_results, f, ensure_ascii=False, indent=2)
        
        self.logger.info(f"所有注释已保存至: {output_file}")
        
        return annotated_results

def main():
    parser = argparse.ArgumentParser(description='为分析结果添加中文注释')
    parser.add_argument('--sample-dir', '-s', required=True,
                      help='样本分析结果目录')
    parser.add_argument('--database', '-d', default='database',
                      help='注释数据库目录（默认: database）')
    parser.add_argument('--output', '-o', help='输出文件路径（可选）')
    
    args = parser.parse_args()
    
    # 创建注释器
    annotator = ChineseAnnotator(database_dir=args.database)
    
    # 处理样本
    results = annotator.process_sample_analysis(args.sample_dir)
    
    # 如果指定了输出路径，额外保存一份
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
        print(f"注释结果已保存至: {output_path}")
    
    print("中文注释完成！")

if __name__ == '__main__':
    main()