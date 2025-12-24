# Bioinformatics Tools Database

A web-based interface for browsing and searching a comprehensive collection of bioinformatics tools, databases, and model services.

## Features

- **Search functionality**: Search by tool name, description, category, or server name
- **Filtering**: Filter tools by category and server name
- **Sorting**: Sort by IDX, tool name, category, or server name
- **Responsive design**: Works on desktop, tablet, and mobile devices
- **Visual indicators**: Color-coded categories for easy scanning

## Deployment to GitHub Pages

1. Create a new GitHub repository
2. Upload all files to the repository:
   - `index.html`
   - `style.css`
   - `script.js`
   - `data.js`
   - `README.md`
3. Go to repository Settings > Pages
4. Under "Source", select "Deploy from a branch"
5. Select the main branch and root folder
6. Click Save
7. Your site will be published at `https://[your-username].github.io/[repository-name]`

## Data Format

The data is stored in `data.js` as a JavaScript array with the following structure:

```javascript
var toolsData = [
    {
        "IDX": 1,
        "Tool Name": "tool_name",
        "Description": "Tool description",
        "category": "Category name",
        "Server Name": "Server name"
    },
    // ... more tools
];
