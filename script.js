// 声明一个变量来存储工具数据，初始为空数组
let toolsData = [];

// 当页面加载完成后初始化应用
document.addEventListener('DOMContentLoaded', function() {
    console.log("DOM加载完成，开始初始化应用...");
    
    // 检查toolsData是否已从data.js加载
    console.log("toolsData长度:", window.toolsData ? window.toolsData.length : "未定义");
    
    // 如果window.toolsData存在，则使用它
    if (window.toolsData) {
        toolsData = window.toolsData;
        console.log("已从data.js加载数据:", toolsData.length, "条记录");
        initApp();
    } else {
        console.error("未找到toolsData，请确保data.js已正确加载");
        // 显示错误信息
        document.getElementById('toolsContainer').innerHTML = `
            <div class="error-message">
                <i class="fas fa-exclamation-triangle"></i>
                <h3>数据加载错误</h3>
                <p>无法加载工具数据。请确保data.js文件存在且格式正确。</p>
                <p>检查浏览器控制台获取更多信息。</p>
            </div>
        `;
    }
});

// 初始化应用
function initApp() {
    console.log("初始化应用...");
    
    // 获取DOM元素
    const toolsContainer = document.getElementById('toolsContainer');
    const searchInput = document.getElementById('searchInput');
    const categoryFilter = document.getElementById('categoryFilter');
    const serverFilter = document.getElementById('serverFilter');
    const resetFiltersBtn = document.getElementById('resetFilters');
    const noResults = document.getElementById('noResults');
    
    // 统计元素
    const totalToolsEl = document.getElementById('totalTools');
    const visibleToolsEl = document.getElementById('visibleTools');
    const serverCountEl = document.getElementById('serverCount');
    const categoryCountEl = document.getElementById('categoryCount');
    
    // 初始化变量
    let filteredTools = [...toolsData];
    let categories = [];
    let servers = [];
    
    // 提取唯一的类别和服务器
    function extractUniqueValues() {
        console.log("提取唯一值...");
        
        // 提取类别
        const categorySet = new Set();
        toolsData.forEach(tool => {
            if (tool.category && tool.category.trim() !== '') {
                categorySet.add(tool.category);
            }
        });
        categories = Array.from(categorySet).sort();
        console.log("找到类别:", categories.length, "个");
        
        // 提取服务器
        const serverSet = new Set();
        toolsData.forEach(tool => {
            if (tool['Server Name'] && tool['Server Name'].trim() !== '') {
                serverSet.add(tool['Server Name']);
            }
        });
        servers = Array.from(serverSet).sort();
        console.log("找到服务器:", servers.length, "个");
    }
    
    // 填充筛选器下拉菜单
    function populateFilters() {
        console.log("填充筛选器...");
        
        // 填充类别筛选器
        categoryFilter.innerHTML = '<option value="">All Categories</option>';
        categories.forEach(category => {
            const option = document.createElement('option');
            option.value = category;
            option.textContent = category;
            categoryFilter.appendChild(option);
        });
        
        // 填充服务器筛选器
        serverFilter.innerHTML = '<option value="">All Servers</option>';
        servers.forEach(server => {
            const option = document.createElement('option');
            option.value = server;
            option.textContent = server;
            serverFilter.appendChild(option);
        });
    }
    
    // 渲染工具卡片
    function renderTools() {
        console.log("渲染工具卡片...", filteredTools.length, "个工具");
        
        toolsContainer.innerHTML = '';
        
        if (filteredTools.length === 0) {
            noResults.style.display = 'block';
            return;
        }
        
        noResults.style.display = 'none';
        
        filteredTools.forEach(tool => {
            const toolCard = createToolCard(tool);
            toolsContainer.appendChild(toolCard);
        });
    }
    
    // 创建工具卡片
    function createToolCard(tool) {
        const card = document.createElement('div');
        card.className = 'tool-card';
        
        // 确保所有字段都有值
        const toolId = tool.IDX || 'N/A';
        const toolName = tool['Tool Name'] || 'Unnamed Tool';
        const description = tool.Description || 'No description available';
        const category = tool.category || 'Uncategorized';
        const serverName = tool['Server Name'] || 'Unknown Server';
        
        // 创建卡片内容
        card.innerHTML = `
            <div class="card-header">
                <span class="tool-id">ID: ${escapeHtml(toolId)}</span>
                <h3 class="tool-name">${escapeHtml(toolName)}</h3>
                <p class="tool-description">${escapeHtml(description)}</p>
            </div>
            <div class="card-body">
                <div class="tool-meta">
                    <span class="category-badge">
                        <i class="fas fa-tag"></i>
                        ${escapeHtml(category)}
                    </span>
                    <span class="server-badge">
                        <i class="fas fa-server"></i>
                        ${escapeHtml(serverName)}
                    </span>
                </div>
                <div class="tool-details">
                    <p><strong>Tool Name:</strong> ${escapeHtml(toolName)}</p>
                    <p><strong>Category:</strong> ${escapeHtml(category)}</p>
                    <p><strong>Server:</strong> ${escapeHtml(serverName)}</p>
                </div>
            </div>
        `;
        
        return card;
    }
    
    // 更新统计信息
    function updateStats() {
        totalToolsEl.textContent = toolsData.length;
        visibleToolsEl.textContent = filteredTools.length;
        serverCountEl.textContent = servers.length;
        categoryCountEl.textContent = categories.length;
        
        console.log("统计信息更新:");
        console.log("- 总工具数:", toolsData.length);
        console.log("- 可见工具数:", filteredTools.length);
        console.log("- 服务器数:", servers.length);
        console.log("- 类别数:", categories.length);
    }
    
    // 过滤工具
    function filterTools() {
        const searchTerm = searchInput.value.toLowerCase();
        const selectedCategory = categoryFilter.value;
        const selectedServer = serverFilter.value;
        
        console.log("过滤工具:", {
            searchTerm,
            selectedCategory,
            selectedServer
        });
        
        filteredTools = toolsData.filter(tool => {
            // 检查搜索词
            const toolName = (tool['Tool Name'] || '').toLowerCase();
            const toolDesc = (tool.Description || '').toLowerCase();
            const matchesSearch = !searchTerm || 
                toolName.includes(searchTerm) ||
                toolDesc.includes(searchTerm);
            
            // 检查类别筛选
            const toolCategory = tool.category || '';
            const matchesCategory = !selectedCategory || 
                toolCategory === selectedCategory;
            
            // 检查服务器筛选
            const toolServer = tool['Server Name'] || '';
            const matchesServer = !selectedServer || 
                toolServer === selectedServer;
            
            return matchesSearch && matchesCategory && matchesServer;
        });
        
        console.log("过滤后剩余:", filteredTools.length, "个工具");
        renderTools();
        updateStats();
    }
    
    // 设置事件监听器
    function setupEventListeners() {
        console.log("设置事件监听器...");
        
        // 搜索输入事件
        searchInput.addEventListener('input', filterTools);
        
        // 筛选器变更事件
        categoryFilter.addEventListener('change', filterTools);
        serverFilter.addEventListener('change', filterTools);
        
        // 重置筛选器按钮
        resetFiltersBtn.addEventListener('click', function() {
            searchInput.value = '';
            categoryFilter.value = '';
            serverFilter.value = '';
            filterTools();
        });
    }
    
    // 辅助函数：转义HTML以防止XSS攻击
    function escapeHtml(text) {
        if (text === null || text === undefined) {
            return '';
        }
        const div = document.createElement('div');
        div.textContent = text.toString();
        return div.innerHTML;
    }
    
    // 执行初始化步骤
    console.log("开始应用初始化...");
    extractUniqueValues();
    populateFilters();
    renderTools();
    updateStats();
    setupEventListeners();
    console.log("应用初始化完成！");
}
