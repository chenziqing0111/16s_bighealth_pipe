#!/usr/bin/env python3
"""
肠道微生物检测报告生成器
集成中文注释，生成增强版HTML报告
"""

import json
import pandas as pd
from pathlib import Path
import argparse
from datetime import datetime
import sys

class ReportGenerator:
    def __init__(self, template_dir='scripts/report'):
        """
        初始化报告生成器
        
        Args:
            template_dir: 模板文件目录
        """
        self.template_dir = Path(template_dir)
        self.report_date = datetime.now().strftime('%Y-%m-%d')
        
    def load_sample_data(self, sample_dir):
        """加载样本的所有分析数据"""
        sample_path = Path(sample_dir)
        data = {
            'sample_id': sample_path.name,
            'report_date': self.report_date
        }
        
        # 加载各个分析模块的结果
        files_to_load = {
            'basic': 'diversity/basic_analysis.json',
            'diversity': 'diversity/basic_analysis.json',
            'enterotype': 'enterotype/enterotype_analysis.json',
            'bacteria': 'bacteria_scores/bacteria_evaluation.json',
            'disease': 'disease_risk/disease_risk_assessment.json',
            'age': 'age_prediction/age_prediction.json',
            'functional': 'functional_prediction/functional_prediction.json',
            'cn_annotations': 'cn_annotations.json'
        }
        
        for key, filepath in files_to_load.items():
            file_path = sample_path / filepath
            if file_path.exists():
                with open(file_path, 'r', encoding='utf-8') as f:
                    data[key] = json.load(f)
                    print(f"  ✓ 加载 {key}: {filepath}")
            else:
                if key not in ['cn_annotations']:
                    print(f"  ⚠ 未找到 {filepath}")
                data[key] = {}
        
        # 从 basic 中提取 composition 数据
        if 'basic' in data and 'composition' in data['basic']:
            data['composition'] = data['basic']['composition']
        
        return data
    
    def load_template(self):
        """加载HTML模板"""
        template_file = self.template_dir / 'report_template.html'
        if not template_file.exists():
            raise FileNotFoundError(f"模板文件不存在: {template_file}")
        
        with open(template_file, 'r', encoding='utf-8') as f:
            template = f.read()
        
        return template
    
    def load_styles(self):
        """加载CSS样式"""
        styles_file = self.template_dir / 'report_styles.css'
        if styles_file.exists():
            with open(styles_file, 'r', encoding='utf-8') as f:
                return f.read()
        return ""
    
    def load_scripts(self):
        """加载JavaScript脚本"""
        scripts_file = self.template_dir / 'report_scripts.js'
        if scripts_file.exists():
            with open(scripts_file, 'r', encoding='utf-8') as f:
                return f.read()
        return ""
    
    def process_data(self, data):
        """处理数据，准备模板变量"""
        # 提取关键数据
        template_vars = {
            'sample_id': data['sample_id'],
            'report_date': data['report_date'],
            
            # 综合评分
            'overall_score': data.get('bacteria', {}).get('overall_health', {}).get('score', 0),
            'health_grade': data.get('bacteria', {}).get('overall_health', {}).get('grade', '未评估'),
            
            # 各项评分
            'beneficial_score': data.get('bacteria', {}).get('beneficial_bacteria', {}).get('overall_score', 0),
            'harmful_score': data.get('bacteria', {}).get('harmful_bacteria', {}).get('harm_score', 0),
            'diversity_score': data.get('basic', {}).get('alpha_diversity', {}).get('shannon', 0),
            
            # 将完整数据转换为JSON传递给JavaScript
            'report_data_json': json.dumps(data, ensure_ascii=False)
        }
        
        return template_vars
    
    def generate_html(self, sample_data, embed_resources=True):
        """
        生成HTML报告
        """
        # 加载模板
        template = self.load_template()
        
        # 处理数据
        template_vars = self.process_data(sample_data)
        
        # 如果需要嵌入资源
        if embed_resources:
            styles = self.load_styles()
            scripts = self.load_scripts()
            
            # 替换外部引用为嵌入式
            template = template.replace(
                '<link rel="stylesheet" href="report_styles.css">',
                f'<style>\n{styles}\n</style>'
            )
            template = template.replace(
                '<script src="report_scripts.js"></script>',
                f'<script>\n{scripts}\n</script>'
            )
        
        # 替换模板变量
        for key, value in template_vars.items():
            if isinstance(value, (int, float)):
                template = template.replace(f'{{{{{key}}}}}', str(round(value, 2)))
            else:
                template = template.replace(f'{{{{{key}}}}}', str(value))
        
        return template
    
    def save_report(self, html_content, output_path):
        """保存HTML报告"""
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✓ 报告已生成: {output}")
        return str(output)
    
    def generate_report(self, sample_dir, output_path, embed_resources=True):
        """
        完整的报告生成流程
        """
        print(f"\n开始生成报告...")
        print(f"  样本目录: {sample_dir}")
        
        # 加载数据
        sample_data = self.load_sample_data(sample_dir)
        
        # 检查疾病风险数据
        if 'disease' in sample_data and sample_data['disease']:
            print(f"  ✓ 疾病风险数据已加载")
            if 'detailed_analysis' in sample_data['disease']:
                disease_count = len(sample_data['disease']['detailed_analysis'])
                print(f"    包含 {disease_count} 种疾病的详细分析")
        else:
            print(f"  ⚠ 疾病风险数据为空")
        
        # 生成HTML
        html_content = self.generate_html(sample_data, embed_resources)
        
        # 保存报告
        report_path = self.save_report(html_content, output_path)
        
        # 生成摘要
        self.generate_summary(sample_data, output_path)
        
        return report_path
    
    def generate_summary(self, sample_data, report_path):
        """生成报告摘要JSON"""
        summary = {
            'sample_id': sample_data['sample_id'],
            'report_date': sample_data['report_date'],
            'report_file': Path(report_path).name,
            'overall_score': sample_data.get('bacteria', {}).get('overall_health', {}).get('score', 0),
            'health_grade': sample_data.get('bacteria', {}).get('overall_health', {}).get('grade', '未评估'),
            'diversity': {
                'shannon': sample_data.get('basic', {}).get('alpha_diversity', {}).get('shannon', 0),
                'observed_asvs': sample_data.get('basic', {}).get('alpha_diversity', {}).get('observed_asvs', 0)
            },
            'bacteria_scores': {
                'beneficial': sample_data.get('bacteria', {}).get('beneficial_bacteria', {}).get('overall_score', 0),
                'harmful': sample_data.get('bacteria', {}).get('harmful_bacteria', {}).get('harm_score', 0)
            },
            'high_risk_diseases': [],
            'has_cn_annotations': bool(sample_data.get('cn_annotations'))
        }
        
        # 提取高风险疾病
        if 'disease' in sample_data and 'risk_assessment' in sample_data['disease']:
            disease_risks = sample_data['disease']['risk_assessment']
            high_risk = [info.get('disease_name', name) for name, info in disease_risks.items() 
                        if info.get('risk_level') == '高风险']
            summary['high_risk_diseases'] = high_risk
        
        # 保存摘要
        summary_path = Path(report_path).with_suffix('.summary.json')
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        
        print(f"✓ 摘要已生成: {summary_path}")
        return summary

def main():
    parser = argparse.ArgumentParser(description='生成肠道微生物检测报告')
    parser.add_argument('--sample-dir', '-s', required=True,
                       help='样本分析结果目录')
    parser.add_argument('--output', '-o', required=True,
                       help='输出HTML文件路径')
    parser.add_argument('--template-dir', '-t', default='scripts/report',
                       help='模板文件目录')
    parser.add_argument('--no-embed', action='store_true',
                       help='不嵌入CSS/JS，使用外部文件引用')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='显示详细信息')
    
    args = parser.parse_args()
    
    try:
        # 创建报告生成器
        generator = ReportGenerator(template_dir=args.template_dir)
        
        # 生成报告
        embed = not args.no_embed
        report_path = generator.generate_report(
            args.sample_dir,
            args.output,
            embed_resources=embed
        )
        
        if args.verbose:
            print(f"\n报告生成成功！")
            print(f"样本: {Path(args.sample_dir).name}")
            print(f"报告: {report_path}")
            print(f"类型: {'单文件' if embed else '多文件'}")
        
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()