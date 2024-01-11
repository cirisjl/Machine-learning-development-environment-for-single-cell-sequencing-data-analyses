import React, { useState, useMemo } from 'react';

const DatasetSelectionDialog = ({ datasets, onSelect, multiple, onClose , isVisible }) => {
    const [selectedDatasets, setSelectedDatasets] = useState([]);

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };
  
    const handleSelectClick = () => {
      onSelect(selectedDatasets);
      onClose(); // Close the dialog after selection
    };
  
    const handleDatasetClick = (dataset) => {
      if (multiple) {
        // For multiple selection, toggle datasets in the array
        setSelectedDatasets(selectedDatasets.includes(dataset)
          ? selectedDatasets.filter((d) => d !== dataset)
          : [...selectedDatasets, dataset]);
      } else {
        // For single selection, set the dataset directly
        setSelectedDatasets([dataset]);
      }
    };

    // const [data, setData] = useState({
    //   initial_data: datasets,
    // });
    // const [sortColumn, setSortColumn] = useState(null);
    // const [sortDirection, setSortDirection] = useState(1);
    // const [visibleColumns, setVisibleColumns] = useState([
    //   'Species',
    //   'Library Construction Method',
    //   'Anatomical Entity',
    //   'Disease Status (Donor)',
    //   'Cell Count Estimate',
    // ]);
  
    // const [filterValues, setFilterValues] = useState({});
    // const [searchTerm, setSearchTerm] = useState('');
    // const [currentPage, setCurrentPage] = useState(1);
    // const [itemsPerPage, setItemsPerPage] = useState(25);
    // const [totalPages, setTotalPages] = useState(1);
    // const [globalSearchTerm, setGlobalSearchTerm] = useState('');
    // const [exportedData, setExportedData] = useState([]);
    // const [error, setError] = useState(null);
  
    // const requiredColumns = [
    //   'Species',
    //   'Library Construction Method',
    //   'Anatomical Entity',
    //   'Disease Status (Donor)',
    //   'Cell Count Estimate',
    // ];
  
    // const columnColors = {
    //   Species: 'btn-info',
    //   'Library Construction Method': 'btn-info',
    //   'Anatomical Entity': 'btn-info',
    //   'Disease Status (Donor)': 'btn-info',
    //   'Cell Count Estimate': 'btn-info',
    // };
  
    // const handleSort = (columnName) => {
    //   try {
    //     const newDirection = columnName === sortColumn ? -sortDirection : 1;
  
    //     const sortedData = [...data.initial_data].sort((a, b) => {
    //       if (a[columnName] < b[columnName]) return -newDirection;
    //       if (a[columnName] > b[columnName]) return newDirection;
    //       return 0;
    //     });
  
    //     setData({
    //       initial_data: sortedData,
    //     });
  
    //     setSortColumn(columnName);
    //     setSortDirection(newDirection);
    //     setError(null);
    //   } catch (error) {
    //     setError(`Error sorting data: ${error.message}`);
    //   }
    // };
  
    // const handleColumnToggle = (columnName) => {
    //   setVisibleColumns((prevVisibleColumns) => {
    //     if (prevVisibleColumns.includes(columnName)) {
    //       return prevVisibleColumns.filter((col) => col !== columnName);
    //     } else {
    //       return [...prevVisibleColumns, columnName];
    //     }
    //   });
    // };
  
    // const handleFilterChange = (column, value) => {
    //   setFilterValues((prevFilterValues) => ({
    //     ...prevFilterValues,
    //     [column]: value,
    //   }));
    //   // Clear global search term when a column filter is changed
    //   setGlobalSearchTerm('');
    // };
  
    // const handleGlobalSearchChange = (value) => {
    //   setGlobalSearchTerm(value);
    // };
  
    // const handleResetFilters = () => {
    //   setGlobalSearchTerm('');
    //   setFilterValues({});
    // };
  
    // const applySearch = (row) => {
    //   return (
    //     (globalSearchTerm === '' ||
    //       requiredColumns.some(
    //         (column) =>
    //           row[column]
    //             .toString()
    //             .toLowerCase()
    //             .includes(globalSearchTerm.toLowerCase())
    //       )) &&
    //     requiredColumns.every((column) => {
    //       const value = row[column];
    //       const filterValue = filterValues[column];
  
    //       if (
    //         filterValue === undefined ||
    //         filterValue === null ||
    //         filterValue === ''
    //       ) {
    //         return true; // No filter for this column
    //       }
  
    //       if (typeof value === 'string' && typeof filterValue === 'string') {
    //         // Case-insensitive search for strings
    //         return value.toLowerCase().includes(filterValue.toLowerCase());
    //       } else if (
    //         typeof value === 'number' &&
    //         typeof filterValue === 'string'
    //       ) {
    //         // Convert numbers to string and then search
    //         return value.toString().includes(filterValue);
    //       }
  
    //       // Non-string and non-number types are not searchable
    //       return false;
    //     })
    //   );
    // };
  
    // const renderCell = (column, value, index) => {
    //   const uniqueKey = `cell_${index}_${column}`;

    //    // Handle different types of values
    //   let displayValue;
    //   if (typeof value === 'string' || typeof value === 'number') {
    //     displayValue = value; // Directly use the value if it's a string or number
    //   } else if (value && typeof value === 'object' && value.label) {
    //     displayValue = value.label; // Use the label property if it's an object with a label
    //   } else {
    //     displayValue = ''; // Fallback to an empty string if the value is not usable
    //   }
  
    //  // Customize display for certain columns
    //   if (['Species', 'Anatomical Entity', 'Disease Status (Donor)'].includes(column)) {
    //     if (displayValue !== '') {
    //       const values = displayValue.split(',');
    //       const count = values.length;
    //       const k = column === 'Anatomical Entity' ? 'Anatomic' : column === 'Disease Status (Donor)' ? 'Disease' : 'Species';

    //       return (
    //         <td key={uniqueKey} data-toggle="tooltip" data-placement="top" title={displayValue}>
    //           {count === 1 ? (
    //             <span>{displayValue}</span>
    //           ) : (
    //             <span className='highlight'>{`${k} (${count})`}</span>
    //           )}
    //         </td>
    //       );
    //     }
    //   } else {
    //     return (
    //       <td key={uniqueKey} data-toggle="tooltip" data-placement="top" title={displayValue}>
    //         {displayValue}
    //       </td>
    //     );
    //   }
    // };
  
    // const handleItemsPerPageChange = (value) => {
    //   setItemsPerPage(value);
    //   setCurrentPage(1);
    // };
  
    // const downloadCSV = (csvContent, fileName) => {
    //   const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    //   const link = document.createElement('a');
  
    //   if (link.download !== undefined) {
    //     const url = URL.createObjectURL(blob);
    //     link.setAttribute('href', url);
    //     link.setAttribute('download', fileName);
    //     document.body.appendChild(link);
    //     link.click();
    //     document.body.removeChild(link);
    //   } else {
    //     alert(
    //       'Your browser does not support the download functionality. Please try a different browser.'
    //     );
    //   }
    // };
  
    // const formatCount = (count) => {
    //   if (count >= 1e6) {
    //     return `${(count / 1e6).toFixed(1)}M`;
    //   } else if (count >= 1e3) {
    //     return `${(count / 1e3).toFixed(1)}k`;
    //   } else {
    //     return count.toString();
    //   }
    // };
  
    // const handleExport = () => {
    //   try {
    //     const exportedData = data.initial_data
    //       .filter(applySearch)
    //       .map((result) =>
    //         visibleColumns.map((column) => result[column]).join(',')
    //       );
  
    //     setExportedData(exportedData);
  
    //     downloadCSV(exportedData.join('\n'), 'exported_data.csv');
    //     setError(null);
    //   } catch (error) {
    //     setError(`Error exporting data: ${error.message}`);
    //   }
    // };
  
    // const paginatedData = useMemo(() => {
    //   const startIndex = (currentPage - 1) * itemsPerPage;
    //   const endIndex = startIndex + itemsPerPage;
    //   setTotalPages(Math.ceil(data.initial_data.length / itemsPerPage));
  
    //   return data.initial_data.slice(startIndex, endIndex);
    // }, [data.initial_data, currentPage, itemsPerPage]);
  
  
    return (
        <div style={dialogStyle} className="dialog-container">
            <div className="dialog-backdrop" onClick={onClose} />
            <div className="dialog">
              <div id="search-panel">
                <div class="metadata-search-wrap">
                    <span class="metadata-search search-title">
                      Search by filters 
                      <a class="action advanced-opts" data-analytics-name="search-help">
                          <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="question-circle" class="svg-inline--fa fa-question-circle " role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
                            <path fill="currentColor" d="M504 256c0 136.997-111.043 248-248 248S8 392.997 8 256C8 119.083 119.043 8 256 8s248 111.083 248 248zM262.655 90c-54.497 0-89.255 22.957-116.549 63.758-3.536 5.286-2.353 12.415 2.715 16.258l34.699 26.31c5.205 3.947 12.621 3.008 16.665-2.122 17.864-22.658 30.113-35.797 57.303-35.797 20.429 0 45.698 13.148 45.698 32.958 0 14.976-12.363 22.667-32.534 33.976C247.128 238.528 216 254.941 216 296v4c0 6.627 5.373 12 12 12h56c6.627 0 12-5.373 12-12v-1.333c0-28.462 83.186-29.647 83.186-106.667 0-58.002-60.165-102-116.531-102zM256 338c-25.365 0-46 20.635-46 46 0 25.364 20.635 46 46 46s46-20.636 46-46c0-25.365-20.635-46-46-46z"></path>
                          </svg>
                      </a>
                    </span>
                    <span class="facet  "><a>organ</a></span><span id="facet-species" class="facet  "><a>species</a></span><span id="facet-disease" class="facet  "><a>disease</a></span><span class="facet  "><a>cell type</a></span>
                    <span id="more-facets-button" class=" facet">
                      <a>
                          <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="sliders-h" class="svg-inline--fa fa-sliders-h icon-left" role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
                            <path fill="currentColor" d="M496 384H160v-16c0-8.8-7.2-16-16-16h-32c-8.8 0-16 7.2-16 16v16H16c-8.8 0-16 7.2-16 16v32c0 8.8 7.2 16 16 16h80v16c0 8.8 7.2 16 16 16h32c8.8 0 16-7.2 16-16v-16h336c8.8 0 16-7.2 16-16v-32c0-8.8-7.2-16-16-16zm0-160h-80v-16c0-8.8-7.2-16-16-16h-32c-8.8 0-16 7.2-16 16v16H16c-8.8 0-16 7.2-16 16v32c0 8.8 7.2 16 16 16h336v16c0 8.8 7.2 16 16 16h32c8.8 0 16-7.2 16-16v-16h80c8.8 0 16-7.2 16-16v-32c0-8.8-7.2-16-16-16zm0-160H288V48c0-8.8-7.2-16-16-16h-32c-8.8 0-16 7.2-16 16v16H16C7.2 64 0 71.2 0 80v32c0 8.8 7.2 16 16 16h208v16c0 8.8 7.2 16 16 16h32c8.8 0 16-7.2 16-16v-16h208c8.8 0 16-7.2 16-16V80c0-8.8-7.2-16-16-16z"></path>
                          </svg>
                          More facets 
                      </a>
                    </span>
                </div>
                <form class="study-keyword-search form-horizontal">
                    <span class="text-search search-title">
                      Search by text 
                      <a class="action advanced-opts" data-analytics-name="search-help">
                          <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="question-circle" class="svg-inline--fa fa-question-circle " role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
                            <path fill="currentColor" d="M504 256c0 136.997-111.043 248-248 248S8 392.997 8 256C8 119.083 119.043 8 256 8s248 111.083 248 248zM262.655 90c-54.497 0-89.255 22.957-116.549 63.758-3.536 5.286-2.353 12.415 2.715 16.258l34.699 26.31c5.205 3.947 12.621 3.008 16.665-2.122 17.864-22.658 30.113-35.797 57.303-35.797 20.429 0 45.698 13.148 45.698 32.958 0 14.976-12.363 22.667-32.534 33.976C247.128 238.528 216 254.941 216 296v4c0 6.627 5.373 12 12 12h56c6.627 0 12-5.373 12-12v-1.333c0-28.462 83.186-29.647 83.186-106.667 0-58.002-60.165-102-116.531-102zM256 338c-25.365 0-46 20.635-46 46 0 25.364 20.635 46 46 46s46-20.636 46-46c0-25.365-20.635-46-46-46z"></path>
                          </svg>
                      </a>
                    </span>
                    <span class="input-group">
                      <input class="form-control" size="30" type="text" placeholder="Search" name="keywordText" value="" fdprocessedid="i3m9c" />
                      <div class="input-group-append">
                          <button type="submit" aria-label="Search keywords" class="btn btn-default" fdprocessedid="9dy9m9">
                            <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="search" class="svg-inline--fa fa-search " role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
                                <path fill="currentColor" d="M505 442.7L405.3 343c-4.5-4.5-10.6-7-17-7H372c27.6-35.3 44-79.7 44-128C416 93.1 322.9 0 208 0S0 93.1 0 208s93.1 208 208 208c48.3 0 92.7-16.4 128-44v16.3c0 6.4 2.5 12.5 7 17l99.7 99.7c9.4 9.4 24.6 9.4 33.9 0l28.3-28.3c9.4-9.4 9.4-24.6.1-34zM208 336c-70.7 0-128-57.2-128-128 0-70.7 57.2-128 128-128 70.7 0 128 57.2 128 128 0 70.7-57.2 128-128 128z"></path>
                            </svg>
                          </button>
                      </div>
                    </span>
                </form>
                <span>
                    <button type="button" id="download-button" class="btn btn-primary" disabled="" aria-disabled="true" aria-label="Download">
                      <span>
                          <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="download" class="svg-inline--fa fa-download icon-left" role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
                            <path fill="currentColor" d="M216 0h80c13.3 0 24 10.7 24 24v168h87.7c17.8 0 26.7 21.5 14.1 34.1L269.7 378.3c-7.5 7.5-19.8 7.5-27.3 0L90.1 226.1c-12.6-12.6-3.7-34.1 14.1-34.1H192V24c0-13.3 10.7-24 24-24zm296 376v112c0 13.3-10.7 24-24 24H24c-13.3 0-24-10.7-24-24V376c0-13.3 10.7-24 24-24h146.7l49 49c20.1 20.1 52.5 20.1 72.6 0l49-49H488c13.3 0 24 10.7 24 24zm-124 88c0-11-9-20-20-20s-20 9-20 20 9 20 20 20 20-9 20-20zm64 0c0-11-9-20-20-20s-20 9-20 20 9 20 20 20 20-9 20-20z"></path>
                          </svg>
                          Download
                      </span>
                    </button>
                </span>
              </div>
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
