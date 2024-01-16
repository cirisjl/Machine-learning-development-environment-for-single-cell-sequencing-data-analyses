import React, { useState,useEffect, useMemo } from 'react';
import FilterComponent from './filtersComponent';
import {SERVER_URL} from '../../../constants/declarations';
import ResultsTable from './tableResultsComponent';
import Pagination from './tablePaginationComponent';

const DatasetSelectionDialog = ({onSelect, multiple, onClose , isVisible }) => {

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };

    const [filters, setFilters] = useState({});
    const [results, setResults] = useState([]);
    const [pagination, setPagination] = useState({});
    const [activeFilters, setActiveFilters] = useState({});


    // Function to fetch data from the API
    const fetchData = async (currentPage, currentFilters) => {

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

    useEffect(() => {   
      fetchData(pagination.page, activeFilters);
    }, [activeFilters]); // Refetch when activeFilters change
  
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

    const onPageChange = (newPage) => {
      setPagination(prev => ({ ...prev, page: newPage }));
    };

    return (
        <div style={dialogStyle} className="dialog-container">
            <div className="dialog-backdrop" onClick={onClose} />
            <div className="dialog">
              <div>
                
                <div className='filters-and search-container'>
                  <div className='filters-container'>
                      {Object.keys(filters).map((filterName) => (
                          <FilterComponent
                              key={filterName}
                              name={filterName}
                              options={filters[filterName]}
                              activeFilters={activeFilters}
                              onFilterChange={handleFilterChange}
                              className="filter"
                          />
                      ))}
                  </div>
                  <div className='search-section'>

                  </div>
                </div>
                
                {/* Results */}
                <ResultsTable data={results} />

                
                {/* Pagination */}
                <Pagination
                  pagination={pagination}
                  onPageChange=
                />
              </div>
              
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
