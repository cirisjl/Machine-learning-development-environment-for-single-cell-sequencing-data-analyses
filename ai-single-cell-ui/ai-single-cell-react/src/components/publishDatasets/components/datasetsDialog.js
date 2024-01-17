import React, { useState,useEffect, useMemo } from 'react';
import FilterComponent from './filtersComponent';
import {SERVER_URL} from '../../../constants/declarations';
import ResultsTable from './tableResultsComponent';
import Pagination from './tablePaginationComponent';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {  faQuestionCircle } from '@fortawesome/free-solid-svg-icons';
import SearchBox from '../../Header/searchBar';

const DatasetSelectionDialog = ({onSelect, multiple, onClose , isVisible }) => {

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };

    const [filters, setFilters] = useState({});
    const [results, setResults] = useState([]);
    const [pagination, setPagination] = useState({});
    const [activeFilters, setActiveFilters] = useState({});

    const [activeFilterCategory, setActiveFilterCategory] = useState(null);

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
    }, [activeFilters, pagination.page]); // Refetch when activeFilters change

    const handleFilterCategoryChange = (category) => {
      setActiveFilterCategory(prevCategory => prevCategory === category ? null : category);
  };
  
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
                
                <div className='filters-and-search-container'>
                  <div className='metadata-search-wrap filters-container'>
                  <span class="metadata-search search-title">Search by filters <FontAwesomeIcon icon={faQuestionCircle}/></span>
                      {Object.keys(filters).map((filterName) => (
                          <FilterComponent
                              key={filterName}
                              name={filterName}
                              options={filters[filterName]}
                              activeFilters={activeFilters}
                              onFilterChange={handleFilterChange}
                              isVisible={activeFilterCategory === filterName}
                              onCategoryChange={() => handleFilterCategoryChange(filterName)}
                              className="filter"
                          />
                      ))}
                  </div>
                  <div className='study-keyword-search'>
                    <span class="text-search search-title">Search by text <FontAwesomeIcon icon={faQuestionCircle} /></span>
                    <SearchBox />
                  </div>

                </div>
                
                {/* Results */}
                <ResultsTable data={results} />

                
                {/* Pagination */}
                <Pagination
                  pagination={pagination}
                  onPageChange={onPageChange}
                />
              </div>
              
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
