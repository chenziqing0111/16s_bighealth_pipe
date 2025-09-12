#!/usr/bin/env python3
"""
肠道微生物检测报告生成器
主脚本：加载分析结果、调用AI解读、生成HTML报告
"""

import json
import argparse
from pathlib import Path
from datetime import datetime
from jinja2 import Template
import pandas as pd
from ai_interpreter import AIInterpreter

class ReportGenerator:
    def __init__(self, api_key=None):
        """初始化报告生成器"""
        self.api_key = api_key
        self.ai_interpreter = AIInterpreter(api_key) if api_key else None
        self.script_dir = Path(__file__).parent
        
    def load_all_results(self, results_dir):
        """加载所有分析模块的结果"""
        results_path = Path(results_dir)
        results = {}
        
        # 加载各模块的JSON结果
        modules = {
            'diversity': 'diversity/basic_analysis.json',
            'enterotype': 'enterotype/enterotype_analysis.json', 
            'bacteria': 'bacteria_scores/bacteria_evaluation.json',
            'disease_risk': 'disease_risk/disease_risk_assessment.json',
            'age': 'age_prediction/age_prediction.json'
        }
        
        for key, path in modules.items():
            file_path = results_path / path
            if file_path.exists():
                with open(file_path, 'r', encoding='utf-8') as f:
                    results[key] = json.load(f)
            else:
                print(f"警告: 未找到 {key} 结果文件: {file_path}")
                results[key] = {}
        
        # 加载样本基本信息（从diversity中提取）
        results['sample_info'] = self._extract_sample_info(results)
        
        return results
    
    def _extract_sample_info(self, results):
        """提取样本基本信息"""
        info = {
            'total_reads': results.get('diversity', {}).get('basic_stats', {}).get('total_reads', 0),
            'total_asvs': results.get('diversity', {}).get('basic_stats', {}).get('total_asvs', 0),
            'analysis_date': datetime.now().strftime('%Y-%m-%d')
        }
        return info
    
    def check_ai_requirements(self, results):
        """判断需要AI解读的内容"""
        ai_needed = []
        
        # 1. 总是需要综合评价
        ai_needed.append('overall_assessment')
        
        # 2. 检查是否有异常指标需要解读
        if self._has_abnormal_bacteria(results):
            ai_needed.append('abnormal_explanation')
        
        # 3. 检查是否有高风险疾病需要解读
        if self._has_high_risk_diseases(results):
            ai_needed.append('disease_interpretation')
        
        # 4. 总是生成个性化建议
        ai_needed.append('personalized_advice')
        
        return ai_needed
    
    def _has_abnormal_bacteria(self, results):
        """检查是否有异常细菌"""
        bacteria = results.get('bacteria', {})
        
        # 检查有益菌是否偏低
        beneficial = bacteria.get('beneficial_bacteria', {}).get('bacteria', {})
        for name, data in beneficial.items():
            if data.get('status') in ['偏低', '偏高']:
                return True
        
        # 检查有害菌是否超标
        harmful = bacteria.get('harmful_bacteria', {}).get('bacteria', {})
        for name, data in harmful.items():
            if data.get('status') == '超标':
                return True
                
        return False
    
    def _has_high_risk_diseases(self, results):
        """检查是否有高风险疾病"""
        diseases = results.get('disease_risk', {}).get('disease_risks', {})
        for disease, data in diseases.items():
            if data.get('risk_score', 0) > 70:
                return True
        return False
    
    def prepare_ai_context(self, results):
        """准备AI解读所需的上下文数据"""
        context = {
            # 多样性数据
            'diversity_score': results.get('diversity', {}).get('alpha_diversity', {}).get('shannon', 0),
            'diversity_status': results.get('diversity', {}).get('alpha_diversity', {}).get('status', '未知'),
            'observed_asvs': results.get('diversity', {}).get('alpha_diversity', {}).get('observed_asvs', 0),
            
            # B/F比值
            'bf_ratio': results.get('diversity', {}).get('bf_ratio', {}).get('value', 0),
            'bf_status': results.get('diversity', {}).get('bf_ratio', {}).get('status', '未知'),
            
            # 肠型
            'enterotype': results.get('enterotype', {}).get('enterotype', {}).get('type', '未知'),
            
            # 细菌评分
            'beneficial_score': results.get('bacteria', {}).get('beneficial_bacteria', {}).get('overall_score', 0),
            'harm_score': results.get('bacteria', {}).get('harmful_bacteria', {}).get('harm_score', 0),
            
            # 异常细菌
            'abnormal_bacteria': self._get_abnormal_bacteria(results),
            
            # 高风险疾病
            'high_risk_diseases': self._get_high_risk_diseases(results),
            
            # 年龄相关
            'biological_age': results.get('age', {}).get('age_prediction', {}).get('biological_age', 0),
            'age_status': results.get('age', {}).get('age_status', {}).get('status', '未知')
        }
        
        return context
    
    def _get_abnormal_bacteria(self, results):
        """获取异常细菌列表"""
        abnormal = []
        
        bacteria = results.get('bacteria', {})
        
        # 有益菌异常
        beneficial = bacteria.get('beneficial_bacteria', {}).get('bacteria', {})
        for name, data in beneficial.items():
            if data.get('status') in ['偏低', '偏高']:
                abnormal.append({
                    'name': name,
                    'type': '有益菌',
                    'status': data.get('status'),
                    'value': data.get('abundance')
                })
        
        # 有害菌异常
        harmful = bacteria.get('harmful_bacteria', {}).get('bacteria', {})
        for name, data in harmful.items():
            if data.get('status') == '超标':
                abnormal.append({
                    'name': name,
                    'type': '有害菌',
                    'status': '超标',
                    'value': data.get('abundance')
                })
        
        return abnormal[:5]  # 最多返回5个
    
    def _get_high_risk_diseases(self, results):
        """获取高风险疾病列表"""
        high_risk = []
        
        diseases = results.get('disease_risk', {}).get('disease_risks', {})
        for disease, data in diseases.items():
            if data.get('risk_score', 0) > 70:
                high_risk.append({
                    'name': disease,
                    'score': data.get('risk_score'),
                    'level': data.get('risk_level')
                })
        
        return sorted(high_risk, key=lambda x: x['score'], reverse=True)[:3]  # 最多返回3个
    
    def generate_report(self, sample_id, results_dir, output_path, metadata=None):
        """生成完整报告"""
        print(f"开始生成报告: {sample_id}")
        
        # 1. 加载分析结果
        print("  加载分析结果...")
        results = self.load_all_results(results_dir)
        
        # 2. 判断需要AI解读的内容
        ai_texts = {}
        if self.ai_interpreter:
            print("  检查AI解读需求...")
            ai_needed = self.check_ai_requirements(results)
            
            if ai_needed:
                print(f"  需要AI解读: {', '.join(ai_needed)}")
                context = self.prepare_ai_context(results)
                ai_texts = self.ai_interpreter.generate_interpretations(context, ai_needed)
        else:
            print("  未配置AI解读，使用默认文本")
            ai_texts = self._get_default_texts()
        
        # 3. 准备模板数据
        template_data = {
            'sample_id': sample_id,
            'metadata': metadata or {},
            'results': results,
            'ai_texts': ai_texts,
            'report_date': datetime.now().strftime('%Y年%m月%d日'),
            'charts_data': self._prepare_charts_data(results)
        }
        
        # 4. 加载并渲染模板
        print("  渲染HTML模板...")
        template_path = self.script_dir / 'report_template.html'
        with open(template_path, 'r', encoding='utf-8') as f:
            template = Template(f.read())
        
        html = template.render(**template_data)
        
        # 5. 内联CSS和JS（便于单文件分发）
        html = self._inline_resources(html)
        
        # 6. 保存报告
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"  报告已生成: {output}")
        return str(output)
    
    def _get_default_texts(self):
        """获取默认文本（无AI时使用）"""
        return {
            'overall_assessment': """
                根据检测结果，您的肠道微生物群落整体处于相对平衡状态。
                建议继续保持良好的饮食习惯，适当增加膳食纤维的摄入，
                有助于维持肠道微生物的多样性和稳定性。
            """,
            'personalized_advice': """
                1. 饮食建议：增加全谷物、蔬菜和水果的摄入
                2. 生活方式：保持规律作息，适度运动
                3. 补充建议：可适当补充益生菌和益生元
            """,
            'abnormal_explanation': """
                检测发现部分细菌指标异常，这可能与近期饮食结构、
                生活压力或用药史有关。建议调整饮食，必要时咨询医生。
            """,
            'disease_interpretation': """
                基于当前的菌群状态，某些疾病风险指标偏高。
                这仅作为健康参考，不能作为临床诊断依据。
                建议定期体检，关注相关健康指标。
            """
        }
    
    def _prepare_charts_data(self, results):
        """准备图表数据"""
        charts = {}
        
        # 多样性数据
        diversity = results.get('diversity', {}).get('alpha_diversity', {})
        charts['diversity'] = {
            'shannon': diversity.get('shannon', 0),
            'simpson': diversity.get('simpson', 0),
            'chao1': diversity.get('chao1', 0),
            'observed_asvs': diversity.get('observed_asvs', 0)
        }
        
        # B/F比值
        bf = results.get('diversity', {}).get('bf_ratio', {})
        charts['bf_ratio'] = {
            'value': bf.get('value', 0),
            'bacteroidetes': bf.get('bacteroidetes', 0),
            'firmicutes': bf.get('firmicutes', 0)
        }
        
        # 物种组成（Top 10）
        composition = results.get('diversity', {}).get('composition', {})
        if 'genus' in composition:
            charts['genus_composition'] = {
                'labels': composition['genus'].get('taxa', [])[:10],
                'values': composition['genus'].get('abundance', [])[:10]
            }
        
        # 疾病风险
        disease_risks = results.get('disease_risk', {}).get('disease_risks', {})
        charts['disease_risks'] = [
            {'name': name, 'score': data.get('risk_score', 0)}
            for name, data in disease_risks.items()
        ]
        
        return charts
    
    def _inline_resources(self, html):
        """内联CSS和JS资源"""
        # 读取CSS
        css_path = self.script_dir / 'report_styles.css'
        if css_path.exists():
            with open(css_path, 'r', encoding='utf-8') as f:
                css_content = f.read()
            html = html.replace(
                '<link rel="stylesheet" href="report_styles.css">',
                f'<style>{css_content}</style>'
            )
        
        # 读取JS
        js_path = self.script_dir / 'report_scripts.js'
        if js_path.exists():
            with open(js_path, 'r', encoding='utf-8') as f:
                js_content = f.read()
            html = html.replace(
                '<script src="report_scripts.js"></script>',
                f'<script>{js_content}</script>'
            )
        
        return html

def main():
    parser = argparse.ArgumentParser(description='生成肠道微生物检测报告')
    parser.add_argument('--sample_id', '-s', required=True, help='样本ID')
    parser.add_argument('--results', '-r', required=True, help='分析结果目录')
    parser.add_argument('--output', '-o', required=True, help='输出HTML文件路径')
    parser.add_argument('--metadata', '-m', help='样本元数据文件')
    parser.add_argument('--api_key', help='DeepSeek API密钥')
    
    args = parser.parse_args()
    
    # 加载元数据
    metadata = None
    if args.metadata and Path(args.metadata).exists():
        metadata = pd.read_csv(args.metadata, sep='\t')
        # 提取当前样本的元数据
        sample_meta = metadata[metadata.iloc[:, 0] == args.sample_id]
        if not sample_meta.empty:
            metadata = sample_meta.iloc[0].to_dict()
    
    # 生成报告
    generator = ReportGenerator(api_key=args.api_key)
    generator.generate_report(
        sample_id=args.sample_id,
        results_dir=args.results,
        output_path=args.output,
        metadata=metadata
    )

if __name__ == '__main__':
    main()
