// 工具数据已经从data.js加载为toolsData变量

// 初始化应用
document.addEventListener('DOMContentLoaded', function() {
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
    
    // 初始化应用
    function initApp() {
        // 提取所有唯一的类别和服务器
        extractUniqueValues();
        
        // 填充筛选器下拉菜单
        populateFilters();
        
        // 渲染工具卡片
        renderTools();
        
        // 更新统计信息
        updateStats();
        
        // 添加事件监听器
        setupEventListeners();
    }
    
    // 提取唯一的类别和服务器
    function extractUniqueValues() {
        // 提取类别
        const categorySet = new Set();
        toolsData.forEach(tool => {
            if (tool.category) {
                categorySet.add(tool.category);
            }
        });
        categories = Array.from(categorySet).sort();
        
        // 提取服务器
        const serverSet = new Set();
        toolsData.forEach(tool => {
            if (tool['Server Name']) {
                serverSet.add(tool['Server Name']);
            }
        });
        servers = Array.from(serverSet).sort();
    }
    
    // 填充筛选器下拉菜单
    function populateFilters() {
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
        
        // 创建卡片内容
        card.innerHTML = `
            <div class="card-header">
                <span class="tool-id">ID: ${tool.IDX}</span>
                <h3 class="tool-name">${escapeHtml(tool['Tool Name'])}</h3>
                <p class="tool-description">${escapeHtml(tool.Description)}</p>
            </div>
            <div class="card-body">
                <div class="tool-meta">
                    <span class="category-badge">
                        <i class="fas fa-tag"></i>
                        ${escapeHtml(tool.category)}
                    </span>
                    <span class="server-badge">
                        <i class="fas fa-server"></i>
                        ${escapeHtml(tool['Server Name'])}
                    </span>
                </div>
                <div class="tool-details">
                    <p><strong>Tool Name:</strong> ${escapeHtml(tool['Tool Name'])}</p>
                    <p><strong>Category:</strong> ${escapeHtml(tool.category)}</p>
                    <p><strong>Server:</strong> ${escapeHtml(tool['Server Name'])}</p>
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
    }
    
    // 过滤工具
    function filterTools() {
        const searchTerm = searchInput.value.toLowerCase();
        const selectedCategory = categoryFilter.value;
        const selectedServer = serverFilter.value;
        
        filteredTools = toolsData.filter(tool => {
            // 检查搜索词
            const matchesSearch = !searchTerm || 
                tool['Tool Name'].toLowerCase().includes(searchTerm) ||
                tool.Description.toLowerCase().includes(searchTerm);
            
            // 检查类别筛选
            const matchesCategory = !selectedCategory || 
                tool.category === selectedCategory;
            
            // 检查服务器筛选
            const matchesServer = !selectedServer || 
                tool['Server Name'] === selectedServer;
            
            return matchesSearch && matchesCategory && matchesServer;
        });
        
        renderTools();
        updateStats();
    }
    
    // 设置事件监听器
    function setupEventListeners() {
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
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }
    
    // 初始化应用
    initApp();
});
