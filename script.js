// Wait for DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {
    // Get DOM elements
    const toolsGrid = document.getElementById('tools-grid');
    const searchInput = document.getElementById('search-input');
    const serverFilter = document.getElementById('server-filter');
    const categoryFilter = document.getElementById('category-filter');
    const resetFiltersBtn = document.getElementById('reset-filters');
    const gridViewBtn = document.getElementById('grid-view');
    const listViewBtn = document.getElementById('list-view');
    const noResults = document.getElementById('no-results');
    const activeFiltersContainer = document.getElementById('active-filters');
    const serverList = document.getElementById('server-list');
    const categoryList = document.getElementById('category-list');
    
    // Stats elements
    const totalTools = document.getElementById('total-tools');
    const totalCount = document.getElementById('total-count');
    const serverCount = document.getElementById('server-count');
    const categoryCount = document.getElementById('category-count');
    const visibleCount = document.getElementById('visible-count');
    const filteredCount = document.getElementById('filtered-count');
    const toolCountFooter = document.getElementById('tool-count-footer');
    const serverCountFooter = document.getElementById('server-count-footer');
    
    // Initialize variables
    let filteredTools = [...toolsData];
    let currentView = 'grid';
    
    // Initialize the application
    function init() {
        // Set initial stats
        updateStats();
        
        // Populate server and category filter dropdowns
        populateFilters();
        
        // Populate server and category lists in sidebar
        populateSidebarLists();
        
        // Display all tools initially
        renderTools();
        
        // Set up event listeners
        setupEventListeners();
    }
    
    // Set up all event listeners
    function setupEventListeners() {
        // Search input
        searchInput.addEventListener('input', filterTools);
        
        // Filter dropdowns
        serverFilter.addEventListener('change', filterTools);
        categoryFilter.addEventListener('change', filterTools);
        
        // Reset filters button
        resetFiltersBtn.addEventListener('click', resetFilters);
        
        // View controls
        gridViewBtn.addEventListener('click', () => switchView('grid'));
        listViewBtn.addEventListener('click', () => switchView('list'));
        
        // Server and category list items in sidebar
        serverList.addEventListener('click', handleServerListClick);
        categoryList.addEventListener('click', handleCategoryListClick);
    }
    
    // Populate server and category filter dropdowns
    function populateFilters() {
        // Get unique servers and categories
        const servers = [...new Set(toolsData.map(tool => tool['Server Name']))].sort();
        const categories = [...new Set(toolsData.map(tool => tool.category))].sort();
        
        // Clear existing options (except the first "All" option)
        while (serverFilter.options.length > 1) {
            serverFilter.remove(1);
        }
        
        while (categoryFilter.options.length > 1) {
            categoryFilter.remove(1);
        }
        
        // Add server options
        servers.forEach(server => {
            const option = document.createElement('option');
            option.value = server;
            option.textContent = server;
            serverFilter.appendChild(option);
        });
        
        // Add category options
        categories.forEach(category => {
            const option = document.createElement('option');
            option.value = category;
            option.textContent = category;
            categoryFilter.appendChild(option);
        });
    }
    
    // Populate server and category lists in sidebar
    function populateSidebarLists() {
        // Clear existing lists
        serverList.innerHTML = '';
        categoryList.innerHTML = '';
        
        // Count tools per server
        const serverCounts = {};
        toolsData.forEach(tool => {
            const server = tool['Server Name'];
            serverCounts[server] = (serverCounts[server] || 0) + 1;
        });
        
        // Count tools per category
        const categoryCounts = {};
        toolsData.forEach(tool => {
            const category = tool.category;
            categoryCounts[category] = (categoryCounts[category] || 0) + 1;
        });
        
        // Sort servers by count (descending) and take top 8
        const topServers = Object.entries(serverCounts)
            .sort((a, b) => b[1] - a[1])
            .slice(0, 8);
        
        // Sort categories by count (descending) and take top 8
        const topCategories = Object.entries(categoryCounts)
            .sort((a, b) => b[1] - a[1])
            .slice(0, 8);
        
        // Populate server list
        topServers.forEach(([server, count]) => {
            const serverItem = document.createElement('div');
            serverItem.className = 'server-item';
            serverItem.dataset.server = server;
            serverItem.innerHTML = `
                <span class="server-name">${server}</span>
                <span class="server-count">${count}</span>
            `;
            serverList.appendChild(serverItem);
        });
        
        // Populate category list
        topCategories.forEach(([category, count]) => {
            const categoryItem = document.createElement('div');
            categoryItem.className = 'category-item';
            categoryItem.dataset.category = category;
            categoryItem.innerHTML = `
                <span class="category-name">${category}</span>
                <span class="category-count">${count}</span>
            `;
            categoryList.appendChild(categoryItem);
        });
    }
    
    // Filter tools based on search input and filters
    function filterTools() {
        const searchTerm = searchInput.value.toLowerCase();
        const selectedServer = serverFilter.value;
        const selectedCategory = categoryFilter.value;
        
        // Filter tools
        filteredTools = toolsData.filter(tool => {
            // Check search term
            const matchesSearch = searchTerm === '' || 
                tool['Tool Name'].toLowerCase().includes(searchTerm) ||
                tool.Description.toLowerCase().includes(searchTerm) ||
                tool['Server Name'].toLowerCase().includes(searchTerm) ||
                tool.category.toLowerCase().includes(searchTerm);
            
            // Check server filter
            const matchesServer = selectedServer === 'all' || tool['Server Name'] === selectedServer;
            
            // Check category filter
            const matchesCategory = selectedCategory === 'all' || tool.category === selectedCategory;
            
            return matchesSearch && matchesServer && matchesCategory;
        });
        
        // Update active filters display
        updateActiveFilters(searchTerm, selectedServer, selectedCategory);
        
        // Render filtered tools
        renderTools();
        
        // Update stats
        updateStats();
    }
    
    // Update active filters display
    function updateActiveFilters(searchTerm, selectedServer, selectedCategory) {
        // Clear current active filters
        activeFiltersContainer.innerHTML = '';
        
        // Add search term filter if present
        if (searchTerm) {
            const searchFilter = document.createElement('div');
            searchFilter.className = 'filter-tag';
            searchFilter.innerHTML = `
                <span>Search: "${searchTerm}"</span>
                <button class="remove-filter" data-type="search"><i class="fas fa-times"></i></button>
            `;
            activeFiltersContainer.appendChild(searchFilter);
        }
        
        // Add server filter if not "all"
        if (selectedServer !== 'all') {
            const serverFilterTag = document.createElement('div');
            serverFilterTag.className = 'filter-tag';
            serverFilterTag.innerHTML = `
                <span>Server: ${selectedServer}</span>
                <button class="remove-filter" data-type="server"><i class="fas fa-times"></i></button>
            `;
            activeFiltersContainer.appendChild(serverFilterTag);
        }
        
        // Add category filter if not "all"
        if (selectedCategory !== 'all') {
            const categoryFilterTag = document.createElement('div');
            categoryFilterTag.className = 'filter-tag';
            categoryFilterTag.innerHTML = `
                <span>Category: ${selectedCategory}</span>
                <button class="remove-filter" data-type="category"><i class="fas fa-times"></i></button>
            `;
            activeFiltersContainer.appendChild(categoryFilterTag);
        }
        
        // If no active filters, show message
        if (activeFiltersContainer.children.length === 0) {
            const noFilters = document.createElement('p');
            noFilters.className = 'no-filters';
            noFilters.textContent = 'No active filters';
            activeFiltersContainer.appendChild(noFilters);
        } else {
            // Add event listeners to remove filter buttons
            const removeButtons = activeFiltersContainer.querySelectorAll('.remove-filter');
            removeButtons.forEach(button => {
                button.addEventListener('click', function() {
                    const filterType = this.dataset.type;
                    
                    if (filterType === 'search') {
                        searchInput.value = '';
                    } else if (filterType === 'server') {
                        serverFilter.value = 'all';
                    } else if (filterType === 'category') {
                        categoryFilter.value = 'all';
                    }
                    
                    filterTools();
                });
            });
        }
    }
    
    // Render tools to the grid
    function renderTools() {
        // Clear current tools
        toolsGrid.innerHTML = '';
        
        // Show/hide no results message
        if (filteredTools.length === 0) {
            noResults.style.display = 'block';
            toolsGrid.style.display = 'none';
            return;
        } else {
            noResults.style.display = 'none';
            toolsGrid.style.display = 'grid';
        }
        
        // Apply view class
        if (currentView === 'list') {
            toolsGrid.classList.add('list-layout');
        } else {
            toolsGrid.classList.remove('list-layout');
        }
        
        // Create tool cards
        filteredTools.forEach(tool => {
            const toolCard = document.createElement('div');
            toolCard.className = `tool-card ${currentView}-view`;
            
            // Truncate description if too long
            let description = tool.Description;
            if (description.length > 200 && currentView === 'grid') {
                description = description.substring(0, 200) + '...';
            }
            
            toolCard.innerHTML = `
                <div class="tool-main">
                    <div class="tool-idx">IDX: ${tool.IDX}</div>
                    <h3 class="tool-name">${tool['Tool Name']}</h3>
                    <p class="tool-description">${description}</p>
                </div>
                <div class="tool-meta">
                    <span class="tool-category">
                        <i class="fas fa-tag"></i> ${tool.category}
                    </span>
                    <span class="tool-server">
                        <i class="fas fa-server"></i> ${tool['Server Name']}
                    </span>
                </div>
            `;
            
            toolsGrid.appendChild(toolCard);
        });
    }
    
    // Switch between grid and list view
    function switchView(view) {
        currentView = view;
        
        // Update active button
        if (view === 'grid') {
            gridViewBtn.classList.add('active');
            listViewBtn.classList.remove('active');
        } else {
            listViewBtn.classList.add('active');
            gridViewBtn.classList.remove('active');
        }
        
        // Re-render tools with new view
        renderTools();
    }
    
    // Reset all filters
    function resetFilters() {
        searchInput.value = '';
        serverFilter.value = 'all';
        categoryFilter.value = 'all';
        
        filterTools();
    }
    
    // Handle server list item clicks in sidebar
    function handleServerListClick(e) {
        const serverItem = e.target.closest('.server-item');
        if (serverItem) {
            const server = serverItem.dataset.server;
            serverFilter.value = server;
            filterTools();
        }
    }
    
    // Handle category list item clicks in sidebar
    function handleCategoryListClick(e) {
        const categoryItem = e.target.closest('.category-item');
        if (categoryItem) {
            const category = categoryItem.dataset.category;
            categoryFilter.value = category;
            filterTools();
        }
    }
    
    // Update statistics
    function updateStats() {
        // Count unique servers and categories
        const uniqueServers = new Set(toolsData.map(tool => tool['Server Name'])).size;
        const uniqueCategories = new Set(toolsData.map(tool => tool.category)).size;
        
        // Update all stat elements
        totalTools.textContent = `${toolsData.length}+`;
        totalCount.textContent = toolsData.length;
        serverCount.textContent = uniqueServers;
        categoryCount.textContent = uniqueCategories;
        visibleCount.textContent = filteredTools.length;
        filteredCount.textContent = `(${filteredTools.length})`;
        toolCountFooter.textContent = toolsData.length;
        serverCountFooter.textContent = uniqueServers;
    }
    
    // Initialize the application
    init();
});
