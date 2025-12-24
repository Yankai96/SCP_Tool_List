// Wait for the DOM to be fully loaded
document.addEventListener('DOMContentLoaded', function() {
    // Get DOM elements
    const searchInput = document.getElementById('searchInput');
    const clearSearchBtn = document.getElementById('clearSearch');
    const categoryFilter = document.getElementById('categoryFilter');
    const serverFilter = document.getElementById('serverFilter');
    const sortSelect = document.getElementById('sortSelect');
    const tableBody = document.getElementById('tableBody');
    const totalCount = document.getElementById('totalCount');
    const filteredCount = document.getElementById('filteredCount');
    const footerCount = document.getElementById('footerCount');
    const noResults = document.getElementById('noResults');
    const tableHeaders = document.querySelectorAll('#toolsTable th');
    
    // Check if toolsData is available
    if (typeof toolsData === 'undefined') {
        console.error('toolsData is not defined. Make sure data.js is loaded correctly.');
        tableBody.innerHTML = '<tr><td colspan="5" style="text-align: center; color: red;">Error: Data not loaded. Please check data.js file.</td></tr>';
        return;
    }
    
    // Initialize variables
    let filteredTools = [...toolsData];
    let currentSort = { column: 'idx', direction: 'asc' };
    
    // Populate filter dropdowns with unique values
    function populateFilters() {
        // Get unique categories
        const categories = [...new Set(toolsData.map(tool => tool.category))].sort();
        categories.forEach(category => {
            const option = document.createElement('option');
            option.value = category;
            option.textContent = category;
            categoryFilter.appendChild(option);
        });
        
        // Get unique server names
        const servers = [...new Set(toolsData.map(tool => tool['Server Name']))].sort();
        servers.forEach(server => {
            const option = document.createElement('option');
            option.value = server;
            option.textContent = server;
            serverFilter.appendChild(option);
        });
    }
    
    // Update counts display
    function updateCounts() {
        const total = toolsData.length;
        const filtered = filteredTools.length;
        
        totalCount.textContent = `Total: ${total} tools`;
        filteredCount.textContent = `Showing: ${filtered} tools`;
        footerCount.textContent = total;
    }
    
    // Render the table with current filtered tools
    function renderTable() {
        // Clear current table content
        tableBody.innerHTML = '';
        
        if (filteredTools.length === 0) {
            noResults.style.display = 'block';
            return;
        }
        
        noResults.style.display = 'none';
        
        // Create table rows for each tool
        filteredTools.forEach(tool => {
            const row = document.createElement('tr');
            
            // Create category badge with appropriate class
            const categoryClass = getCategoryClass(tool.category);
            const categoryBadge = `<span class="category-badge ${categoryClass}">${tool.category}</span>`;
            
            // Create server badge
            const serverBadge = `<span class="server-badge">${tool['Server Name']}</span>`;
            
            row.innerHTML = `
                <td>${tool.IDX}</td>
                <td><strong>${tool['Tool Name']}</strong></td>
                <td>${tool.Description}</td>
                <td>${categoryBadge}</td>
                <td>${serverBadge}</td>
            `;
            
            tableBody.appendChild(row);
        });
        
        updateCounts();
    }
    
    // Get CSS class for category badge
    function getCategoryClass(category) {
        if (category.includes('Database')) return 'category-Databases';
        if (category.includes('Computational')) return 'category-Computational';
        if (category.includes('Model')) return 'category-Model';
        if (category.includes('Literature')) return 'category-Literature';
        if (category.includes('Wet-lab')) return 'category-Wet-lab';
        return 'category-Databases'; // default
    }
    
    // Filter tools based on search and filter criteria
    function filterTools() {
        const searchTerm = searchInput.value.toLowerCase().trim();
        const selectedCategory = categoryFilter.value;
        const selectedServer = serverFilter.value;
        
        filteredTools = toolsData.filter(tool => {
            // Search filter
            const matchesSearch = searchTerm === '' || 
                tool['Tool Name'].toLowerCase().includes(searchTerm) ||
                tool.Description.toLowerCase().includes(searchTerm) ||
                tool.category.toLowerCase().includes(searchTerm) ||
                tool['Server Name'].toLowerCase().includes(searchTerm);
            
            // Category filter
            const matchesCategory = selectedCategory === '' || tool.category === selectedCategory;
            
            // Server filter
            const matchesServer = selectedServer === '' || tool['Server Name'] === selectedServer;
            
            return matchesSearch && matchesCategory && matchesServer;
        });
        
        // Apply current sort
        sortTools();
        renderTable();
    }
    
    // Sort tools based on current sort column and direction
    function sortTools() {
        filteredTools.sort((a, b) => {
            let aValue, bValue;
            
            switch (currentSort.column) {
                case 'name':
                    aValue = a['Tool Name'].toLowerCase();
                    bValue = b['Tool Name'].toLowerCase();
                    break;
                case 'category':
                    aValue = a.category.toLowerCase();
                    bValue = b.category.toLowerCase();
                    break;
                case 'server':
                    aValue = a['Server Name'].toLowerCase();
                    bValue = b['Server Name'].toLowerCase();
                    break;
                case 'idx':
                default:
                    aValue = a.IDX;
                    bValue = b.IDX;
                    break;
            }
            
            // Compare values
            if (aValue < bValue) return currentSort.direction === 'asc' ? -1 : 1;
            if (aValue > bValue) return currentSort.direction === 'asc' ? 1 : -1;
            return 0;
        });
    }
    
    // Update sort UI indicators
    function updateSortIndicators() {
        // Remove all sort classes from headers
        tableHeaders.forEach(header => {
            header.classList.remove('sort-asc', 'sort-desc');
        });
        
        // Add sort class to appropriate header
        const headerIndex = currentSort.column === 'idx' ? 0 : 
                          currentSort.column === 'name' ? 1 :
                          currentSort.column === 'category' ? 3 :
                          currentSort.column === 'server' ? 4 : 0;
        
        if (tableHeaders[headerIndex]) {
            tableHeaders[headerIndex].classList.add(`sort-${currentSort.direction}`);
        }
    }
    
    // Event Listeners
    
    // Search input
    searchInput.addEventListener('input', function() {
        filterTools();
        clearSearchBtn.style.display = this.value ? 'block' : 'none';
    });
    
    // Clear search button
    clearSearchBtn.addEventListener('click', function() {
        searchInput.value = '';
        filterTools();
        this.style.display = 'none';
        searchInput.focus();
    });
    
    // Category filter
    categoryFilter.addEventListener('change', filterTools);
    
    // Server filter
    serverFilter.addEventListener('change', filterTools);
    
    // Sort select
    sortSelect.addEventListener('change', function() {
        const sortValue = this.value;
        
        if (sortValue === 'name') {
            currentSort.column = 'name';
        } else if (sortValue === 'category') {
            currentSort.column = 'category';
        } else if (sortValue === 'server') {
            currentSort.column = 'server';
        } else {
            currentSort.column = 'idx';
        }
        
        sortTools();
        renderTable();
        updateSortIndicators();
    });
    
    // Table header click for sorting
    tableHeaders.forEach((header, index) => {
        header.addEventListener('click', function() {
            let column;
            
            switch (index) {
                case 0: column = 'idx'; break;
                case 1: column = 'name'; break;
                case 3: column = 'category'; break;
                case 4: column = 'server'; break;
                default: return;
            }
            
            // If clicking the same column, toggle direction
            if (currentSort.column === column) {
                currentSort.direction = currentSort.direction === 'asc' ? 'desc' : 'asc';
            } else {
                // Otherwise set to new column with ascending direction
                currentSort.column = column;
                currentSort.direction = 'asc';
            }
            
            // Update sort select to match
            sortSelect.value = column;
            
            // Sort and render
            sortTools();
            renderTable();
            updateSortIndicators();
        });
    });
    
    // Initialize the page
    function init() {
        populateFilters();
        sortTools();
        renderTable();
        updateSortIndicators();
        
        // Add keyboard shortcut for search (Ctrl/Cmd + F)
        document.addEventListener('keydown', function(e) {
            if ((e.ctrlKey || e.metaKey) && e.key === 'f') {
                e.preventDefault();
                searchInput.focus();
            }
            
            // Escape key clears search
            if (e.key === 'Escape' && document.activeElement === searchInput) {
                searchInput.value = '';
                filterTools();
                clearSearchBtn.style.display = 'none';
            }
        });
    }
    
    // Start the application
    init();
});
