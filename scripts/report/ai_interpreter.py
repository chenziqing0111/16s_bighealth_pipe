#!/usr/bin/env python3
"""
AI解读模块：使用DeepSeek API生成个性化文本
"""

import json
import time
from typing import Dict, List, Any
import requests

class AIInterpreter:
    def __init__(self, api_key: str, model: str = "deepseek-chat"):
        """初始化AI解读器"""
        self.api_key = api_key
        self.model = model
        self.api_url = "https://api.deepseek.com/v1/chat/completions"
        
        # 内置提示词模板
        self.prompts = {
            'overall_assessment': {
                'system': "你是一位专业的肠道微生物检测报告解读专家，擅长用通俗易懂的语言解释复杂的检测结果。",
                'template': """
                基于以下肠道微生物检测数据，生成一段200字左右的综合健康评价：
                
                检测数据：
                - Shannon多样性指数：{diversity_score}（状态：{diversity_status}）
                - 观察到的ASVs数量：{observed_asvs}
                - B/F比值：{bf_ratio}（状态：{bf_status}）
                - 肠型：{enterotype}
                - 有益菌评分：{beneficial_score}/100
                - 有害菌评分：{harm_score}/100
                - 生物年龄：{biological_age}岁（{age_status}）
                
                要求：
                1. 首先评价整体肠道健康状况（优秀/良好/一般/需改善）
                2. 指出主要的优势和问题
                3. 简要说明这些指标的健康意义
                4. 语言专业但通俗易懂
                5. 不要使用过于绝对的表述
                """
            },
            
            'abnormal_explanation': {
                'system': "你是一位肠道微生物专家，能够准确解释细菌异常的原因和影响。",
                'template': """
                以下细菌检测结果出现异常，请解释可能的原因和健康影响（150字以内）：
                
                异常细菌：
                {abnormal_bacteria}
                
                要求：
                1. 解释这些细菌异常的可能原因（饮食、生活方式、用药等）
                2. 说明对健康的潜在影响
                3. 避免引起过度担忧
                4. 提供改善方向
                """
            },
            
            'disease_interpretation': {
                'system': "你是一位预防医学专家，擅长解读疾病风险并提供预防建议。",
                'template': """
                基于肠道菌群检测，以下疾病风险偏高，请解释其关联性（200字以内）：
                
                高风险疾病：
                {high_risk_diseases}
                
                相关菌群状态：
                - 有益菌评分：{beneficial_score}/100
                - 有害菌评分：{harm_score}/100
                - 主要异常菌：{abnormal_bacteria}
                
                要求：
                1. 解释菌群与这些疾病风险的关联
                2. 强调这只是风险评估，不是诊断
                3. 提供预防建议的方向
                4. 语气客观，避免制造恐慌
                """
            },
            
            'personalized_advice': {
                'system': "你是一位营养和健康管理专家，擅长根据肠道菌群状态提供个性化建议。",
                'template': """
                根据以下肠道菌群检测结果，生成个性化的健康改善建议（300字以内）：
                
                主要问题：
                - 异常细菌：{abnormal_bacteria}
                - 高风险疾病：{high_risk_diseases}
                - 肠型：{enterotype}
                - 多样性状态：{diversity_status}
                
                要求：
                1. 饮食建议（具体到食物类别，3-5条）
                2. 生活方式建议（作息、运动等，2-3条）
                3. 补充剂建议（益生菌、益生元等，如需要）
                4. 建议要具体、可执行
                5. 按优先级排列
                6. 考虑肠型特点给出针对性建议
                """
            }
        }
    
    def generate_interpretations(self, context: Dict[str, Any], needed: List[str]) -> Dict[str, str]:
        """批量生成所需的AI解读文本"""
        results = {}
        
        for key in needed:
            if key in self.prompts:
                print(f"    生成{key}...")
                text = self._generate_single_interpretation(key, context)
                results[key] = text
                
                # 避免API限流
                time.sleep(1)
        
        return results
    
    def _generate_single_interpretation(self, interpretation_type: str, context: Dict[str, Any]) -> str:
        """生成单个解读文本"""
        prompt_config = self.prompts[interpretation_type]
        
        # 格式化用户提示词
        user_prompt = prompt_config['template'].format(**self._format_context(context))
        
        # 构建请求
        messages = [
            {"role": "system", "content": prompt_config['system']},
            {"role": "user", "content": user_prompt}
        ]
        
        try:
            # 调用API
            response = self._call_api(messages)
            return self._clean_response(response)
        except Exception as e:
            print(f"      AI生成失败: {e}")
            return self._get_fallback_text(interpretation_type)
    
    def _format_context(self, context: Dict[str, Any]) -> Dict[str, str]:
        """格式化上下文数据为字符串"""
        formatted = {}
        
        for key, value in context.items():
            if isinstance(value, list):
                # 格式化列表（如异常细菌列表）
                if key == 'abnormal_bacteria':
                    formatted[key] = '\n'.join([
                        f"- {item['name']}（{item['type']}）：{item['status']}，检测值{item['value']}%"
                        for item in value
                    ])
                elif key == 'high_risk_diseases':
                    formatted[key] = '\n'.join([
                        f"- {item['name']}：风险评分{item['score']}分（{item['level']}）"
                        for item in value
                    ])
                else:
                    formatted[key] = ', '.join(str(v) for v in value)
            else:
                formatted[key] = str(value)
        
        return formatted
    
    def _call_api(self, messages: List[Dict[str, str]]) -> str:
        """调用DeepSeek API"""
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        
        data = {
            "model": self.model,
            "messages": messages,
            "temperature": 0.7,
            "max_tokens": 500
        }
        
        response = requests.post(
            self.api_url,
            headers=headers,
            json=data,
            timeout=30
        )
        
        if response.status_code == 200:
            result = response.json()
            return result['choices'][0]['message']['content']
        else:
            raise Exception(f"API请求失败: {response.status_code} - {response.text}")
    
    def _clean_response(self, text: str) -> str:
        """清理API响应文本"""
        # 去除多余的空白
        text = text.strip()
        
        # 去除可能的编号开头
        lines = text.split('\n')
        cleaned_lines = []
        for line in lines:
            line = line.strip()
            # 去除"1. "这样的编号
            if line and line[0].isdigit() and line[1:3] == '. ':
                line = line[3:]
            cleaned_lines.append(line)
        
        return '\n'.join(cleaned_lines)
    
    def _get_fallback_text(self, interpretation_type: str) -> str:
        """获取备用文本（API失败时使用）"""
        fallbacks = {
            'overall_assessment': """
                您的肠道微生物检测显示整体健康状况良好。菌群多样性处于正常范围，
                有益菌和有害菌的平衡基本正常。建议继续保持健康的生活方式，
                定期进行肠道健康检测，及时了解肠道微生态的变化。
            """,
            
            'abnormal_explanation': """
                检测发现部分细菌指标出现异常。这可能与近期的饮食习惯、
                生活压力或用药史有关。建议调整饮食结构，增加膳食纤维摄入，
                必要时可在医生指导下补充益生菌。
            """,
            
            'disease_interpretation': """
                基于肠道菌群组成，某些疾病的风险指标偏高。请注意，
                这仅是基于菌群特征的风险评估，不能作为疾病诊断的依据。
                建议结合其他检查结果，必要时咨询专业医生。
            """,
            
            'personalized_advice': """
                饮食建议：增加全谷物、新鲜蔬果的摄入，减少高脂高糖食物。
                生活建议：保持规律作息，每周至少150分钟中等强度运动。
                补充建议：可适当补充益生菌和益生元，促进肠道健康。
            """
        }
        
        return fallbacks.get(interpretation_type, "暂无相关解读。")

# 测试函数
def test_interpreter():
    """测试AI解读器（需要提供API密钥）"""
    # 模拟上下文数据
    test_context = {
        'diversity_score': 3.2,
        'diversity_status': '正常',
        'observed_asvs': 450,
        'bf_ratio': 2.5,
        'bf_status': '正常',
        'enterotype': '拟杆菌型',
        'beneficial_score': 65,
        'harm_score': 30,
        'biological_age': 45,
        'age_status': '年轻态',
        'abnormal_bacteria': [
            {'name': 'Bifidobacterium', 'type': '有益菌', 'status': '偏低', 'value': 0.5}
        ],
        'high_risk_diseases': [
            {'name': 'IBD', 'score': 75, 'level': '高风险'}
        ]
    }
    
    # 需要生成的解读
    needed = ['overall_assessment', 'personalized_advice']
    
    # 初始化解读器（需要真实的API密钥）
    interpreter = AIInterpreter(api_key="your_api_key_here")
    
    # 生成解读
    results = interpreter.generate_interpretations(test_context, needed)
    
    for key, text in results.items():
        print(f"\n{key}:")
        print(text)

if __name__ == '__main__':
    # 运行测试
    # test_interpreter()
    pass
