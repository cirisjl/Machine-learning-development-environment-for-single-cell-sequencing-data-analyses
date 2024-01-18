import React, { useState,useEffect, useMemo } from 'react';
import FilterComponent from './filtersComponent';
import {SERVER_URL} from '../../../constants/declarations';
import ResultsTable from './tableResultsComponent';
import Pagination from './tablePaginationComponent';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {  faQuestionCircle, faSliders } from '@fortawesome/free-solid-svg-icons';
import SearchBox from '../../Header/searchBar';

const DatasetSelectionDialog = ({onSelect, multiple, onClose , isVisible }) => {

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };

    const [visibleFacets, setVisibleFacets] = useState([]); // Will hold the keys of the facets to display
    const [showMoreFacets, setShowMoreFacets] = useState(false); // Toggle state for showing more facets

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

    useEffect(() => {
      // Set initial visible facets to the first four, or fewer if there aren't four
      setVisibleFacets(Object.keys(filters).slice(0, 4));
  }, [filters]); // This will update visible facets when the filters are fetched


  const toggleMoreFacets = () => {
    setShowMoreFacets(!showMoreFacets);
    // Show all facets if showMoreFacets is true, else show only the first four
    if (!showMoreFacets) {
        setVisibleFacets(Object.keys(filters));
    } else {
        setVisibleFacets(Object.keys(filters).slice(0, 4));
    }
};


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
                      {visibleFacets.map((filterName) => (
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
                      <div>
                        {Object.keys(filters).length > 4 && (
                            <button onClick={toggleMoreFacets} className='filters-toggle'>
                                <FontAwesomeIcon icon={faSliders} /> <p>{showMoreFacets ? 'Less facets' : 'More facets'}</p>
                            </button>
                        )}
                      </div>
                  </div>
                  <div className='study-keyword-search'>
                    <span class="text-search search-title">Search by text <FontAwesomeIcon icon={faQuestionCircle} /></span>
                    <SearchBox />
                  </div>

                </div>

                <div className='table-pagination'>
                  <Pagination
                    pagination={pagination}
                    onPageChange={onPageChange}
                  />
                </div>
                
                <div className='table-results'>
                  <ResultsTable data={results} />
                </div>
              </div>
              
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
