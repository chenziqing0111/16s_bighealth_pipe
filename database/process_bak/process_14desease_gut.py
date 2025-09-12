#!/usr/bin/env python3
"""
处理gutMD.csv文件，提取14种重点疾病的菌群关联数据
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from collections import defaultdict

def process_gutmd_for_14_diseases(input_file='gutMD.csv', output_dir='database'):
    """处理gutMD数据库，提取14种疾病的关联"""
    
    # 14种疾病的映射关系
    disease_mapping = {
        # 消化系统疾病
        'IBD': {
            'name_cn': '炎症性肠病',
            'keywords': [
                'Inflammatory Bowel Disease', 'IBD',
                'Crohn Disease', "Crohn's Disease", 'Crohn',
                'Colitis,Ulcerative', 'Ulcerative Colitis', 
                'Active ulcerative colitis', 'Pediatric ulcerative colitis',
                'Active Crohn', 'Inactive Crohn', 'Pediatric Crohn'
            ]
        },
        'IBS': {
            'name_cn': '肠易激综合征',
            'keywords': [
                'Irritable Bowel Syndrome', 'IBS',
                'Diarrhea-Predominant Irritable', 
                'Constipated-Irritable Bowel',
                'mixed type'
            ]
        },
        'Constipation': {
            'name_cn': '便秘',
            'keywords': ['Constipation', 'Chronic Constipation']
        },
        'CRC': {
            'name_cn': '结直肠癌/息肉',
            'keywords': [
                'Colorectal Neoplasms', 'CRC', 'Colorectal Cancer',
                'Adenoma', 'Polyp', 'adenoma', 'polyp',
                'tubular adenoma', 'serrated', 'hyperplastic'
            ]
        },
        'Celiac': {
            'name_cn': '乳糜泻',
            'keywords': ['Celiac Disease', 'Gluten', 'Diet,Gluten-Free']
        },
        
        # 代谢性疾病
        'Diabetes': {
            'name_cn': '糖尿病',
            'keywords': [
                'Diabetes Mellitus', 'Type 2 Diabetes', 'T2D',
                'Type 1 Diabetes', 'T1D', 'Prediabetic',
                'Gestational Diabetes', 'Beta cell autoimmunity'
            ]
        },
        'Obesity': {
            'name_cn': '肥胖/代谢综合征',
            'keywords': [
                'Obesity', 'Obese', 'Overweight', 
                'Metabolic Syndrome', 'MetS', 'BMI>33',
                'Weight Loss', 'Gastric Bypass'
            ]
        },
        'NAFLD': {
            'name_cn': '非酒精性脂肪肝',
            'keywords': [
                'Non-alcoholic Fatty Liver', 'NAFLD', 
                'Fatty Liver', 'Liver Disease',
                'Hepatitis', 'Cirrhosis', 'hepatocellular'
            ]
        },
        'CVD': {
            'name_cn': '心血管疾病',
            'keywords': [
                'Cardiovascular', 'CVD', 'Coronary Disease',
                'Heart Disease', 'Myocardial', 'Atherosclerosis',
                'Heart Failure', 'Angina'
            ]
        },
        
        # 神经精神疾病
        'Depression': {
            'name_cn': '抑郁症',
            'keywords': [
                'Depression', 'Depressive', 'Major Depression',
                'MDD', 'Depressed'
            ]
        },
        'Autism': {
            'name_cn': '自闭症',
            'keywords': [
                'Autism', 'ASD', 'Autism Spectrum Disorder',
                'Autistic'
            ]
        },
        'Parkinson': {
            'name_cn': '帕金森病',
            'keywords': [
                'Parkinson', "Parkinson's Disease", 'PD',
                'Parkinsonian'
            ]
        },
        
        # 免疫相关疾病
        'Eczema': {
            'name_cn': '湿疹/特应性皮炎',
            'keywords': [
                'Eczema', 'Atopic Dermatitis', 'Dermatitis',
                'Atopic', 'Atopy'
            ]
        },
        'RA': {
            'name_cn': '类风湿性关节炎',
            'keywords': [
                'Rheumatoid Arthritis', 'RA', 'Arthritis,Rheumatoid',
                'Arthritis, Rheumatoid'
            ]
        }
    }
    
    # 读取数据
    print(f"读取文件: {input_file}")
    df = pd.read_csv(input_file)
    df.columns = df.columns.str.strip()
    
    print(f"数据列: {df.columns.tolist()}")
    print(f"总记录数: {len(df)}")
    
    # 存储结果
    disease_associations = {}
    enhanced_database = {}
    
    # 处理每种疾病
    for disease_key, disease_info in disease_mapping.items():
        print(f"\n处理 {disease_key} ({disease_info['name_cn']})...")
        
        associations = {
            'beneficial': defaultdict(lambda: {'count': 0, 'pmids': set(), 'level': None}),
            'harmful': defaultdict(lambda: {'count': 0, 'pmids': set(), 'level': None}),
            'evidence_count': 0,
            'pmids': set(),
            'matched_conditions': set()
        }
        
        # 搜索所有关键词
        for keyword in disease_info['keywords']:
            # 查找疾病vs健康的对比
            disease_records = df[
                (df['Condition1'].str.contains(keyword, case=False, na=False, regex=False)) &
                (df['Condition2'].str.contains('Health|Control|Normal', case=False, na=False, regex=True))
            ]
            
            # 记录匹配到的条件
            if len(disease_records) > 0:
                associations['matched_conditions'].update(disease_records['Condition1'].unique())
            
            # 处理每条记录
            for _, row in disease_records.iterrows():
                bacteria_info = row.get('Gut Microbiota (ID)', row.get('Gut Microbiota', ''))
                
                if pd.isna(bacteria_info) or bacteria_info == '':
                    continue
                
                # 提取菌名
                if '(' in str(bacteria_info):
                    bacteria_name = bacteria_info.split('(')[0].strip()
                else:
                    bacteria_name = str(bacteria_info).strip()
                
                if not bacteria_name or bacteria_name == 'nan':
                    continue
                
                # 清理菌名
                bacteria_name = bacteria_name.replace('lactobacillus', 'Lactobacillus')
                bacteria_name = bacteria_name.replace('bifidobacterium', 'Bifidobacterium')
                
                level = row.get('Classification', 'unknown')
                alteration = str(row.get('Alteration', '')).lower()
                pmid = str(row.get('PMID', ''))
                
                # 分类到有益或有害
                if alteration in ['decrease', 'absent', 'reduced', 'lower']:
                    associations['beneficial'][bacteria_name]['count'] += 1
                    associations['beneficial'][bacteria_name]['pmids'].add(pmid)
                    associations['beneficial'][bacteria_name]['level'] = level
                    
                elif alteration in ['increase', 'present', 'higher', 'elevated']:
                    associations['harmful'][bacteria_name]['count'] += 1
                    associations['harmful'][bacteria_name]['pmids'].add(pmid)
                    associations['harmful'][bacteria_name]['level'] = level
                
                associations['evidence_count'] += 1
                associations['pmids'].add(pmid)
        
        # 排序并选择最相关的菌
        beneficial_sorted = sorted(
            associations['beneficial'].items(),
            key=lambda x: x[1]['count'],
            reverse=True
        )[:15]
        
        harmful_sorted = sorted(
            associations['harmful'].items(),
            key=lambda x: x[1]['count'],
            reverse=True
        )[:15]
        
        # 生成简化版本
        beneficial_list = [name for name, _ in beneficial_sorted]
        harmful_list = [name for name, _ in harmful_sorted]
        
        # 根据证据数量设置权重
        if associations['evidence_count'] > 100:
            weight = 2.0
        elif associations['evidence_count'] > 50:
            weight = 1.5
        elif associations['evidence_count'] > 20:
            weight = 1.2
        else:
            weight = 1.0
        
        # 保存结果
        enhanced_database[disease_key] = {
            'name_cn': disease_info['name_cn'],
            'beneficial': beneficial_list,
            'harmful': harmful_list,
            'weight': weight,
            'evidence_count': associations['evidence_count'],
            'pmid_count': len(associations['pmids']),
            'matched_conditions': list(associations['matched_conditions'])[:5]
        }
        
        # 详细统计
        disease_associations[disease_key] = {
            'name_cn': disease_info['name_cn'],
            'beneficial_details': {
                name: {
                    'count': data['count'],
                    'level': data['level'],
                    'pmids': list(data['pmids'])[:5]
                }
                for name, data in beneficial_sorted
            },
            'harmful_details': {
                name: {
                    'count': data['count'],
                    'level': data['level'],
                    'pmids': list(data['pmids'])[:5]
                }
                for name, data in harmful_sorted
            },
            'statistics': {
                'total_evidence': associations['evidence_count'],
                'unique_pmids': len(associations['pmids']),
                'beneficial_bacteria_count': len(beneficial_list),
                'harmful_bacteria_count': len(harmful_list)
            }
        }
        
        # 打印统计
        print(f"  - 证据数: {associations['evidence_count']}")
        print(f"  - 有益菌: {len(beneficial_list)}")
        print(f"  - 有害菌: {len(harmful_list)}")
        if associations['matched_conditions']:
            print(f"  - 匹配条件示例: {list(associations['matched_conditions'])[:3]}")
    
    # 保存文件
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # 1. 完整版（JSON）
    with open(output_path / 'disease_associations_full.json', 'w', encoding='utf-8') as f:
        json.dump(disease_associations, f, indent=2, ensure_ascii=False)
    
    # 2. 简化版（JSON）
    with open(output_path / 'disease_associations.json', 'w', encoding='utf-8') as f:
        json.dump(enhanced_database, f, indent=2, ensure_ascii=False)
    
    # 3. Python代码版本
    with open(output_path / 'disease_db_code.py', 'w', encoding='utf-8') as f:
        f.write("#!/usr/bin/env python3\n")
        f.write('"""从gutMD提取的14种疾病菌群关联数据库"""\n\n')
        f.write("DISEASE_DATABASE = {\n")
        for disease_key, data in enhanced_database.items():
            f.write(f"    '{disease_key}': {{\n")
            f.write(f"        'name_cn': '{data['name_cn']}',\n")
            f.write(f"        'beneficial': {data['beneficial']},\n")
            f.write(f"        'harmful': {data['harmful']},\n")
            f.write(f"        'weight': {data['weight']},\n")
            f.write(f"        'evidence_count': {data['evidence_count']}\n")
            f.write("    },\n")
        f.write("}\n")
    
    # 4. 生成统计报告
    with open(output_path / 'statistics_report.txt', 'w', encoding='utf-8') as f:
        f.write("14种疾病菌群关联数据统计报告\n")
        f.write("="*50 + "\n\n")
        
        total_evidence = 0
        for disease_key, data in enhanced_database.items():
            f.write(f"{disease_key} ({data['name_cn']})\n")
            f.write(f"  证据数: {data['evidence_count']}\n")
            f.write(f"  文献数: {data['pmid_count']}\n")
            f.write(f"  有益菌: {', '.join(data['beneficial'][:5])}...\n")
            f.write(f"  有害菌: {', '.join(data['harmful'][:5])}...\n")
            f.write("-"*30 + "\n")
            total_evidence += data['evidence_count']
        
        f.write(f"\n总计: {total_evidence} 条证据\n")
    
    print("\n" + "="*50)
    print("处理完成！文件已保存到:", output_path)
    print("- disease_associations_full.json (完整版)")
    print("- disease_associations.json (简化版)")
    print("- disease_db_code.py (Python代码)")
    print("- statistics_report.txt (统计报告)")
    
    # 打印总结
    print("\n=== 14种疾病统计摘要 ===")
    for disease_key, data in enhanced_database.items():
        status = "✓" if data['evidence_count'] > 0 else "✗"
        print(f"{status} {disease_key:12} ({data['name_cn']:10}): {data['evidence_count']:4}条证据")
    
    return enhanced_database

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='处理gutMD数据库提取14种疾病关联')
    parser.add_argument('--input', '-i', default='gutMD.csv', help='输入文件路径')
    parser.add_argument('--output', '-o', default='database', help='输出目录')
    
    args = parser.parse_args()
    
    # 执行处理
    result = process_gutmd_for_14_diseases(args.input, args.output)