import React, { useState,useEffect, useMemo } from 'react';
import FilterComponent from './filtersComponent';
import {SERVER_URL} from '../../../constants/declarations';

const DatasetSelectionDialog = ({onSelect, multiple, onClose , isVisible }) => {
    // const [selectedDatasets, setSelectedDatasets] = useState([]);

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };

    const [filters, setFilters] = useState({});
    const [results, setResults] = useState([]);
    const [pagination, setPagination] = useState({});
    const [activeFilters, setActiveFilters] = useState({});

    useEffect(() => {
      // Function to fetch data from the API
      const fetchData = async () => {
        // const query = new URLSearchParams();
        
        // for (const [key, values] of Object.entries(activeFilters)) {
        //   values.forEach(value => query.append(key, value));
        // }
        
        try {
          const response = await fetch(`${SERVER_URL}/api/datasets/search`);
          const data = await response.json();
          setFilters(data.facets);
          setResults(data.results);
          setPagination(data.pagination);
        } catch (error) {
          console.error('Error fetching data:', error);
        }
      };      
  
      fetchData();
    }, [activeFilters]); // Refetch when activeFilters change
  
    // const handleSelectClick = () => {
    //   onSelect(selectedDatasets);
    //   onClose(); // Close the dialog after selection
    // };
  
    // const handleDatasetClick = (dataset) => {
    //   if (multiple) {
    //     // For multiple selection, toggle datasets in the array
    //     setSelectedDatasets(selectedDatasets.includes(dataset)
    //       ? selectedDatasets.filter((d) => d !== dataset)
    //       : [...selectedDatasets, dataset]);
    //   } else {
    //     // For single selection, set the dataset directly
    //     setSelectedDatasets([dataset]);
    //   }
    // };

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

    const handleFilterChange = (filterCategory, filterValue) => {
      setActiveFilters(prevFilters => {
        const newFilters = { ...prevFilters };
        if (newFilters[filterCategory]) {
          if (newFilters[filterCategory].includes(filterValue)) {
            // Filter is already active, so remove it
            newFilters[filterCategory] = newFilters[filterCategory].filter(v => v !== filterValue);
          } else {
            // Add the filter to the active category
            newFilters[filterCategory].push(filterValue);
          }
        } else {
          // This category hasn't been selected yet, so add a new array with this value
          newFilters[filterCategory] = [filterValue];
        }
    
        // If the filter category array is empty, remove the category from the filters
        if (newFilters[filterCategory].length === 0) {
          delete newFilters[filterCategory];
        }
    
        return newFilters;
      });
    };
    
  
  
    return (
        <div style={dialogStyle} className="dialog-container">
            <div className="dialog-backdrop" onClick={onClose} />
            <div className="dialog">
              <div>
                
                <div className='filters-and search-container'>
                  <div className='filters-section'>
                    {/* Filters */}
                    {Object.keys(filters).map((filterName) => (
                      <FilterComponent
                        key={filterName}
                        name={filterName}
                        options={filters[filterName]}
                        activeFilters={activeFilters}
                        onFilterChange={handleFilterChange}
                      />
                    ))}
                  </div>
                  <div className='search-section'>

                  </div>
                </div>
                
                {/* Results */}
                {/* <ResultsList results={results} /> */}
                
                {/* Pagination */}
                {/* <Pagination
                  pagination={pagination}
                  onPageChange=
                /> */}
              </div>
              
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
