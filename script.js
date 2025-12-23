document.addEventListener('DOMContentLoaded', function() {
    // 全局变量
    let allTools = [];
    let filteredTools = [];
    let currentPage = 1;
    const toolsPerPage = 50;
    let currentSortColumn = 'IDX';
    let currentSortDirection = 'asc';
    
    // DOM 元素
    const tableBody = document.getElementById('table-body');
    const searchInput = document.getElementById('search-input');
    const clearSearchBtn = document.getElementById('clear-search');
    const categoryFilter = document.getElementById('category-filter');
    const totalToolsSpan = document.getElementById('total-tools');
    const shownToolsSpan = document.getElementById('shown-tools');
    const pageNumSpan = document.getElementById('page-num');
    const totalPagesSpan = document.getElementById('total-pages');
    const prevPageBtn = document.getElementById('prev-page');
    const nextPageBtn = document.getElementById('next-page');
    const totalRowsSpan = document.getElementById('total-rows');
    const noResultsDiv = document.getElementById('no-results');
    const updateDateSpan = document.getElementById('update-date');
    
    // 初始化 - 从data.js加载数据
    function initialize() {
        // 设置更新日期为当前日期
        const now = new Date();
        updateDateSpan.textContent = `${now.getFullYear()}年${now.getMonth() + 1}月${now.getDate()}日`;
        
        // 从全局变量tools获取数据（来自data.js）
        if (typeof tools !== 'undefined' && tools.length > 0) {
            allTools = tools;
            totalRowsSpan.textContent = allTools.length;
            totalToolsSpan.textContent = allTools.length;
            
            // 填充分类筛选器
            populateCategoryFilter();
            
            // 初始显示所有工具
            filteredTools = [...allTools];
            
            // 排序并显示
            sortTools(currentSortColumn, currentSortDirection);
            updateDisplay();
            
            // 设置表头排序事件
            setupSorting();
        } else {
            console.error('未找到工具数据。请确保data.js已正确加载。');
            tableBody.innerHTML = '<tr><td colspan="4" style="text-align: center; color: red;">数据加载失败，请检查data.js文件</td></tr>';
        }
    }
    
    // 填充分类筛选器
    function populateCategoryFilter() {
        const categories = new Set();
        allTools.forEach(tool => {
            if (tool.category) {
                categories.add(tool.category);
            }
        });
        
        // 清空现有选项（除了"所有分类"）
        while (categoryFilter.options.length > 1) {
            categoryFilter.remove(1);
        }
        
        // 按字母顺序排序并添加选项
        Array.from(categories).sort().forEach(category => {
            const option = document.createElement('option');
            option.value = category;
            option.textContent = category;
            categoryFilter.appendChild(option);
        });
    }
    
    // 设置表头排序
    function setupSorting() {
        const headers = document.querySelectorAll('#tools-table th');
        headers.forEach(header => {
            header.addEventListener('click', function() {
                const column = this.getAttribute('data-sort');
                if (column === currentSortColumn) {
                    // 切换排序方向
                    currentSortDirection = currentSortDirection === 'asc' ? 'desc' : 'asc';
                } else {
                    // 新列，默认升序
                    currentSortColumn = column;
                    currentSortDirection = 'asc';
                }
                
                // 更新排序图标
                headers.forEach(h => {
                    const icon = h.querySelector('i');
                    if (h.getAttribute('data-sort') === currentSortColumn) {
                        icon.className = currentSortDirection === 'asc' ? 'fas fa-sort-up' : 'fas fa-sort-down';
                    } else {
                        icon.className = 'fas fa-sort';
                    }
                });
                
                // 排序并显示
                sortTools(currentSortColumn, currentSortDirection);
                updateDisplay();
            });
        });
    }
    
    // 排序工具
    function sortTools(column, direction) {
        filteredTools.sort((a, b) => {
            let aValue = a[column];
            let bValue = b[column];
            
            // 处理可能的空值
            if (aValue === undefined || aValue === null) aValue = '';
            if (bValue === undefined || bValue === null) bValue = '';
            
            // 如果是ID列，按数字排序
            if (column === 'IDX') {
                aValue = parseInt(aValue) || 0;
                bValue = parseInt(bValue) || 0;
                return direction === 'asc' ? aValue - bValue : bValue - aValue;
            }
            
            // 其他列按字符串排序
            aValue = String(aValue).toLowerCase();
            bValue = String(bValue).toLowerCase();
            
            if (aValue < bValue) return direction === 'asc' ? -1 : 1;
            if (aValue > bValue) return direction === 'asc' ? 1 : -1;
            return 0;
        });
    }
    
    // 过滤工具（搜索和分类筛选）
    function filterTools() {
        const searchTerm = searchInput.value.trim().toLowerCase();
        const selectedCategory = categoryFilter.value;
        
        filteredTools = allTools.filter(tool => {
            // 分类筛选
            if (selectedCategory && tool.category !== selectedCategory) {
                return false;
            }
            
            // 搜索筛选
            if (searchTerm) {
                const toolName = (tool['Tool Name'] || '').toLowerCase();
                const description = (tool.Description || '').toLowerCase();
                const category = (tool.category || '').toLowerCase();
                
                return toolName.includes(searchTerm) || 
                       description.includes(searchTerm) || 
                       category.includes(searchTerm);
            }
            
            return true;
        });
        
        // 重置到第一页
        currentPage = 1;
        
        // 重新排序
        sortTools(currentSortColumn, currentSortDirection);
        updateDisplay();
    }
    
    // 更新显示
    function updateDisplay() {
        // 更新工具数量统计
        shownToolsSpan.textContent = filteredTools.length;
        
        // 计算分页
        const totalPages = Math.ceil(filteredTools.length / toolsPerPage);
        totalPagesSpan.textContent = totalPages;
        
        // 更新分页按钮状态
        prevPageBtn.disabled = currentPage <= 1;
        nextPageBtn.disabled = currentPage >= totalPages;
        
        // 如果没有结果，显示提示
        if (filteredTools.length === 0) {
            tableBody.innerHTML = '';
            noResultsDiv.style.display = 'block';
            pageNumSpan.textContent = '0';
            return;
        } else {
            noResultsDiv.style.display = 'none';
        }
        
        // 确保当前页面有效
        if (currentPage > totalPages && totalPages > 0) {
            currentPage = totalPages;
        }
        
        pageNumSpan.textContent = currentPage;
        
        // 获取当前页的工具
        const startIndex = (currentPage - 1) * toolsPerPage;
        const endIndex = Math.min(startIndex + toolsPerPage, filteredTools.length);
        const pageTools = filteredTools.slice(startIndex, endIndex);
        
        // 清空表格
        tableBody.innerHTML = '';
        
        // 填充表格
        pageTools.forEach(tool => {
            const row = document.createElement('tr');
            
            // ID 列
            const idCell = document.createElement('td');
            idCell.textContent = tool.IDX || '';
            row.appendChild(idCell);
            
            // 工具名称列
            const nameCell = document.createElement('td');
            const toolName = tool['Tool Name'] || '';
            nameCell.textContent = toolName;
            row.appendChild(nameCell);
            
            // 描述列
            const descCell = document.createElement('td');
            const description = tool.Description || '';
            descCell.textContent = description;
            row.appendChild(descCell);
            
            // 分类列
            const categoryCell = document.createElement('td');
            const category = tool.category || '';
            
            if (category) {
                const badge = document.createElement('span');
                badge.textContent = category;
                badge.className = 'category-badge ' + getCategoryClass(category);
                categoryCell.appendChild(badge);
            }
            
            row.appendChild(categoryCell);
            
            tableBody.appendChild(row);
        });
        
        // 高亮搜索关键词
        if (searchInput.value.trim()) {
            highlightSearchTerms(searchInput.value.trim());
        }
    }
    
    // 获取分类对应的CSS类
    function getCategoryClass(category) {
        const lowerCategory = category.toLowerCase();
        if (lowerCategory.includes('computational')) return 'category-computational';
        if (lowerCategory.includes('database')) return 'category-databases';
        if (lowerCategory.includes('model')) return 'category-models';
        if (lowerCategory.includes('literature')) return 'category-literature';
        if (lowerCategory.includes('wet') || lowerCategory.includes('lab')) return 'category-wetlab';
        return 'category-computational'; // 默认
    }
    
    // 高亮搜索关键词
    function highlightSearchTerms(searchTerm) {
        const searchTerms = searchTerm.toLowerCase().split(' ').filter(term => term.length > 0);
        
        if (searchTerms.length === 0) return;
        
        const rows = tableBody.querySelectorAll('tr');
        
        rows.forEach(row => {
            const cells = row.querySelectorAll('td');
            cells.forEach((cell, index) => {
                // 跳过ID列
                if (index === 0) return;
                
                const originalText = cell.textContent;
                let highlightedText = originalText;
                
                searchTerms.forEach(term => {
                    if (term.length < 2) return;
                    
                    const regex = new RegExp(`(${term})`, 'gi');
                    highlightedText = highlightedText.replace(regex, '<span class="highlight">$1</span>');
                });
                
                if (highlightedText !== originalText) {
                    cell.innerHTML = highlightedText;
                }
            });
        });
    }
    
    // 事件监听器
    searchInput.addEventListener('input', function() {
        clearSearchBtn.style.visibility = this.value ? 'visible' : 'hidden';
        filterTools();
    });
    
    clearSearchBtn.addEventListener('click', function() {
        searchInput.value = '';
        clearSearchBtn.style.visibility = 'hidden';
        filterTools();
        searchInput.focus();
    });
    
    categoryFilter.addEventListener('change', filterTools);
    
    prevPageBtn.addEventListener('click', function() {
        if (currentPage > 1) {
            currentPage--;
            updateDisplay();
        }
    });
    
    nextPageBtn.addEventListener('click', function() {
        const totalPages = Math.ceil(filteredTools.length / toolsPerPage);
        if (currentPage < totalPages) {
            currentPage++;
            updateDisplay();
        }
    });
    
    // 键盘快捷键
    document.addEventListener('keydown', function(e) {
        // Ctrl+F 或 Cmd+F 聚焦搜索框
        if ((e.ctrlKey || e.metaKey) && e.key === 'f') {
            e.preventDefault();
            searchInput.focus();
            searchInput.select();
        }
        
        // Escape 清除搜索
        if (e.key === 'Escape' && document.activeElement === searchInput) {
            searchInput.value = '';
            filterTools();
        }
    });
    
    // 初始化应用
    initialize();
});
