// ============ 导航功能 ============
function toggleNav() {
    const sidebar = document.getElementById('navSidebar');
    const mainContent = document.getElementById('mainContent');
    sidebar.classList.toggle('active');
    mainContent.classList.toggle('nav-active');
}

function scrollToSection(sectionId) {
    const section = document.getElementById(sectionId);
    if (section) {
        section.scrollIntoView({ behavior: 'smooth', block: 'start' });
        
        // 更新活动导航项
        document.querySelectorAll('.nav-item').forEach(item => {
            item.classList.remove('active');
        });
        event.target.classList.add('active');
    }
}

// ============ 数据渲染 ============
document.addEventListener('DOMContentLoaded', function() {
    // 渲染综合评分
    renderOverallScore();
    
    // 渲染多样性数据
    renderDiversity();
    
    // 渲染细菌分析（使用中文注释）
    renderBacteriaAnalysis();
    
    // 渲染营养代谢（使用中文注释）
    renderNutritionMetabolism();
    
    // 渲染疾病风险
    renderDiseaseRisk();
    
    // 渲染免疫力评估
    renderImmunity();
    
    // 渲染饮食建议
    renderRecommendations();
    
    // 初始化图表
    initCharts();
});

// 渲染综合评分
function renderOverallScore() {
    const score = reportData.bacteria?.overall_health?.score || 0;
    const grade = reportData.bacteria?.overall_health?.grade || '未评估';
    
    // 设置评分圆圈颜色
    const scoreCircle = document.getElementById('overallScoreCircle');
    if (scoreCircle) {
        if (score >= 80) {
            scoreCircle.className = 'score-circle score-good';
        } else if (score >= 60) {
            scoreCircle.className = 'score-circle score-warning';
        } else {
            scoreCircle.className = 'score-circle score-danger';
        }
    }
    
    // 设置健康等级徽章
    const gradeElement = document.getElementById('overallGrade');
    if (gradeElement) {
        gradeElement.textContent = grade;
        if (grade.includes('优秀') || grade.includes('良好')) {
            gradeElement.className = 'status-badge status-normal';
        } else if (grade.includes('亚健康')) {
            gradeElement.className = 'status-badge status-warning';
        } else {
            gradeElement.className = 'status-badge status-danger';
        }
    }
}

// 渲染多样性数据
function renderDiversity() {
    const diversity = reportData.diversity?.alpha_diversity || {};
    
    // 创建多样性指标表格
    const metricsHtml = `
        <table class="table table-striped">
            <tbody>
                <tr>
                    <td>Shannon指数</td>
                    <td>${(diversity.shannon || 0).toFixed(3)}</td>
                    <td><span class="status-badge status-normal">正常</span></td>
                </tr>
                <tr>
                    <td>Simpson指数</td>
                    <td>${(diversity.simpson || 0).toFixed(3)}</td>
                    <td><span class="status-badge status-normal">正常</span></td>
                </tr>
                <tr>
                    <td>Chao1估计</td>
                    <td>${(diversity.chao1 || 0).toFixed(1)}</td>
                    <td><span class="status-badge status-normal">正常</span></td>
                </tr>
                <tr>
                    <td>观察到的ASVs</td>
                    <td>${diversity.observed_asvs || 0}</td>
                    <td><span class="status-badge ${diversity.status === '正常' ? 'status-normal' : 'status-warning'}">${diversity.status || '正常'}</span></td>
                </tr>
            </tbody>
        </table>
    `;
    
    const metricsDiv = document.getElementById('diversityMetrics');
    if (metricsDiv) {
        metricsDiv.innerHTML = metricsHtml;
    }
}

// 渲染细菌分析（使用中文注释）
function renderBacteriaAnalysis() {
    // 优先使用中文注释数据
    const bacteria = reportData.cn_annotations?.bacteria || reportData.bacteria || {};
    
    // 渲染有益菌卡片
    const beneficial = bacteria.beneficial_bacteria?.bacteria || {};
    const beneficialCards = document.getElementById('beneficialBacteriaCards');
    if (beneficialCards) {
        let cardsHtml = '';
        Object.entries(beneficial).slice(0, 6).forEach(([name, info]) => {
            const cnName = info.cn_name || name;
            const description = info.description || '益生菌';
            const abundance = info.abundance || 0;
            const status = info.status || '未知';
            const statusClass = status === '正常' ? 'status-normal' : 'status-warning';
            
            cardsHtml += `
                <div class="col-md-6">
                    <div class="bacteria-card">
                        <h5 class="bacteria-name">${cnName}</h5>
                        <p class="bacteria-latin">${name}</p>
                        <p class="bacteria-description">${description}</p>
                        <div class="d-flex justify-content-between mt-2">
                            <span>丰度: ${abundance.toFixed(3)}%</span>
                            <span class="status-badge ${statusClass}">${status}</span>
                        </div>
                    </div>
                </div>
            `;
        });
        beneficialCards.innerHTML = cardsHtml;
    }
    
    // 渲染有害菌表格
    const harmful = bacteria.harmful_bacteria?.bacteria || {};
    const harmfulTable = document.querySelector('#harmfulBacteriaTable tbody');
    if (harmfulTable) {
        let tableHtml = '';
        Object.entries(harmful).forEach(([name, info]) => {
            const cnName = info.cn_name || name;
            const abundance = info.abundance || 0;
            const threshold = info.threshold || 0.5;
            const status = info.status || '未知';
            const statusClass = status === '正常' ? 'status-normal' : 'status-danger';
            
            tableHtml += `
                <tr>
                    <td>
                        ${cnName}
                        <br><small class="text-muted">${name}</small>
                    </td>
                    <td>${abundance.toFixed(4)}%</td>
                    <td>&lt; ${threshold}%</td>
                    <td><span class="status-badge ${statusClass}">${status}</span></td>
                </tr>
            `;
        });
        harmfulTable.innerHTML = tableHtml;
    }
}

// 渲染营养代谢（使用中文注释）
function renderNutritionMetabolism() {
    const nutrition = reportData.cn_annotations?.nutrition || reportData.functional || {};
    const pathways = reportData.cn_annotations?.pathways || {};
    
    // 渲染维生素合成
    const vitamins = nutrition.vitamin_synthesis || {};
    const vitaminDiv = document.getElementById('vitaminSynthesis');
    if (vitaminDiv) {
        let vitaminHtml = '';
        Object.entries(vitamins).forEach(([vitamin, info]) => {
            const cnName = info.cn_name || `维生素${vitamin}`;
            const potential = info.synthesis_potential || 0;
            const status = info.status || '未知';
            const color = status === '高' ? '#4CAF50' : status === '中等' ? '#FF9800' : '#F44336';
            
            vitaminHtml += `
                <div class="col-md-4">
                    <div class="metric-card">
                        <h6>${cnName}</h6>
                        <div class="progress">
                            <div class="progress-bar" style="width: ${Math.min(100, potential * 10)}%; background: ${color}">
                                ${potential.toFixed(1)}
                            </div>
                        </div>
                        <small>${status}合成能力</small>
                    </div>
                </div>
            `;
        });
        vitaminDiv.innerHTML = vitaminHtml;
    }
    
    // 渲染短链脂肪酸
    const scfa = nutrition.scfa_production || {};
    const scfaDiv = document.getElementById('scfaProduction');
    if (scfaDiv) {
        let scfaHtml = '';
        Object.entries(scfa).forEach(([name, info]) => {
            const cnName = info.cn_name || name;
            const potential = info.production_potential || 0;
            const status = info.status || '未知';
            
            scfaHtml += `
                <div class="col-md-4">
                    <div class="metric-card">
                        <div class="metric-icon"><i class="fas fa-flask"></i></div>
                        <h6>${cnName}</h6>
                        <div class="metric-value">${potential.toFixed(1)}</div>
                        <div class="metric-label">${status}产生能力</div>
                    </div>
                </div>
            `;
        });
        scfaDiv.innerHTML = scfaHtml;
    }
    
    // 渲染代谢通路表格
    const pathwaysTable = document.querySelector('#pathwaysTable tbody');
    if (pathwaysTable && pathways) {
        let tableHtml = '';
        const sortedPathways = Object.entries(pathways)
            .sort((a, b) => (b[1].abundance || 0) - (a[1].abundance || 0))
            .slice(0, 5);
        
        sortedPathways.forEach(([id, info]) => {
            const cnName = info.cn_name || id;
            const category = info.category || '其他';
            const abundance = info.abundance || 0;
            const importance = info.importance || 'medium';
            const importanceClass = importance === 'high' ? 'status-danger' : 
                                   importance === 'medium' ? 'status-warning' : 'status-info';
            
            tableHtml += `
                <tr>
                    <td>
                        ${cnName}
                        <br><small class="text-muted">${id}</small>
                    </td>
                    <td>${category}</td>
                    <td>${abundance.toFixed(0)}</td>
                    <td><span class="status-badge ${importanceClass}">${importance}</span></td>
                </tr>
            `;
        });
        pathwaysTable.innerHTML = tableHtml;
    }
}

// 渲染疾病风险
function renderDiseaseRisk() {
    const disease = reportData.disease?.risk_assessment || {};
    const diseaseCards = document.getElementById('diseaseRiskCards');
    
    if (diseaseCards) {
        let cardsHtml = '';
        Object.entries(disease).forEach(([name, info]) => {
            const riskScore = info.risk_score || 0;
            const riskLevel = info.risk_level || '未知';
            
            let color, icon;
            if (riskLevel === '低风险') {
                color = '#4CAF50';
                icon = 'fa-check-circle';
            } else if (riskLevel === '中风险') {
                color = '#FF9800';
                icon = 'fa-exclamation-triangle';
            } else {
                color = '#F44336';
                icon = 'fa-times-circle';
            }
            
            cardsHtml += `
                <div class="col-md-4">
                    <div class="metric-card" style="border-color: ${color}">
                        <i class="fas ${icon}" style="color: ${color}; font-size: 2rem;"></i>
                        <h6 class="mt-2">${name}</h6>
                        <div style="color: ${color}; font-size: 1.5rem; font-weight: bold;">
                            ${riskScore.toFixed(1)}%
                        </div>
                        <small>${riskLevel}</small>
                    </div>
                </div>
            `;
        });
        diseaseCards.innerHTML = cardsHtml;
    }
}

// 渲染免疫力评估
function renderImmunity() {
    const beneficialScore = reportData.bacteria?.beneficial_bacteria?.overall_score || 0;
    const immunityScore = Math.min(100, beneficialScore * 1.2);
    
    let immunityLevel, color;
    if (immunityScore >= 80) {
        immunityLevel = '优秀';
        color = '#4CAF50';
    } else if (immunityScore >= 60) {
        immunityLevel = '良好';
        color = '#66BB6A';
    } else if (immunityScore >= 40) {
        immunityLevel = '一般';
        color = '#FF9800';
    } else {
        immunityLevel = '较弱';
        color = '#F44336';
    }
    
    const scoreElement = document.getElementById('immunityScore');
    const levelElement = document.getElementById('immunityLevel');
    
    if (scoreElement) scoreElement.textContent = immunityScore.toFixed(0) + '分';
    if (levelElement) {
        levelElement.textContent = `免疫力水平: ${immunityLevel}`;
        levelElement.style.color = color;
    }
    
    // 渲染影响因素
    const factorsDiv = document.getElementById('immunityFactors');
    if (factorsDiv) {
        factorsDiv.innerHTML = `
            <h5>影响因素分析</h5>
            <ul>
                <li>有益菌水平: ${beneficialScore.toFixed(1)}分</li>
                <li>菌群多样性: ${reportData.diversity?.alpha_diversity?.status || '正常'}</li>
                <li>代谢功能: 正常</li>
            </ul>
        `;
    }
}

// 渲染饮食建议
function renderRecommendations() {
    const recommendations = reportData.bacteria?.recommendations || [];
    
    // 基于分析结果生成个性化建议
    const dietRecommendations = [
        {
            category: '增加摄入',
            items: ['酸奶', '发酵食品', '全谷物', '豆类', '新鲜蔬菜'],
            reason: '促进有益菌生长',
            color: '#4CAF50'
        },
        {
            category: '适量摄入',
            items: ['瘦肉', '鸡蛋', '坚果', '水果'],
            reason: '平衡营养',
            color: '#FF9800'
        },
        {
            category: '减少摄入',
            items: ['高脂食物', '精制糖', '加工食品', '酒精'],
            reason: '减少有害菌滋生',
            color: '#F44336'
        }
    ];
    
    const dietDiv = document.getElementById('dietRecommendations');
    if (dietDiv) {
        let html = '';
        dietRecommendations.forEach(rec => {
            const itemsList = rec.items.map(item => `<li>${item}</li>`).join('');
            html += `
                <div class="col-md-4">
                    <div class="recommendation-card" style="border-left-color: ${rec.color}">
                        <h5 style="color: ${rec.color}">${rec.category}</h5>
                        <ul style="list-style: none; padding-left: 0;">
                            ${itemsList}
                        </ul>
                        <small class="text-muted">${rec.reason}</small>
                    </div>
                </div>
            `;
        });
        dietDiv.innerHTML = html;
    }
    
    // 生活方式建议
    const lifestyleDiv = document.getElementById('lifestyleAdvice');
    if (lifestyleDiv) {
        lifestyleDiv.innerHTML = `
            <ul>
                <li>保持规律作息，每天睡眠7-8小时</li>
                <li>适量运动，每周至少150分钟中等强度运动</li>
                <li>管理压力，可尝试冥想或瑜伽</li>
                <li>定期复查，监测肠道健康状况</li>
                <li>避免滥用抗生素，保护肠道菌群</li>
            </ul>
        `;
    }
}

// 初始化图表
function initCharts() {
    // 多样性雷达图
    const diversityCtx = document.getElementById('diversityChart');
    if (diversityCtx && typeof Chart !== 'undefined') {
        const diversity = reportData.diversity?.alpha_diversity || {};
        
        new Chart(diversityCtx.getContext('2d'), {
            type: 'radar',
            data: {
                labels: ['Shannon', 'Simpson', 'Chao1', 'ASVs', 'Evenness'],
                datasets: [{
                    label: '多样性指标',
                    data: [
                        (diversity.shannon || 0) / 10 * 100,
                        (diversity.simpson || 0) * 100,
                        (diversity.chao1 || 0) / 2000 * 100,
                        (diversity.observed_asvs || 0) / 2000 * 100,
                        (diversity.evenness || 0.5) * 100
                    ],
                    backgroundColor: 'rgba(102, 187, 106, 0.2)',
                    borderColor: 'rgba(102, 187, 106, 1)',
                    pointBackgroundColor: 'rgba(102, 187, 106, 1)',
                    pointBorderColor: '#fff',
                    pointHoverBackgroundColor: '#fff',
                    pointHoverBorderColor: 'rgba(102, 187, 106, 1)'
                }]
            },
            options: {
                scales: {
                    r: {
                        beginAtZero: true,
                        max: 100,
                        ticks: {
                            stepSize: 20
                        }
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    }
                }
            }
        });
    }
}

// 自动调整导航栏高亮
window.addEventListener('scroll', function() {
    const sections = document.querySelectorAll('.report-section');
    const navItems = document.querySelectorAll('.nav-item');
    
    let current = '';
    sections.forEach(section => {
        const rect = section.getBoundingClientRect();
        if (rect.top <= 100 && rect.bottom >= 100) {
            current = section.getAttribute('id');
        }
    });
    
    navItems.forEach(item => {
        item.classList.remove('active');
        if (item.getAttribute('onclick')?.includes(current)) {
            item.classList.add('active');
        }
    });
});