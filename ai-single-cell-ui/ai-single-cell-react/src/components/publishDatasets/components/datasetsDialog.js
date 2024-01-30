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
    const [globalSearchTerm, setGlobalSearchTerm] = useState('');

    const [activeFilterCategory, setActiveFilterCategory] = useState(null);
    const [appliedFilters, setAppliedFilters] = useState([]);

    // Function to fetch data from the API
    const fetchData = async (currentPage, currentFilters, searchQuery) => {

      try {
        const response = await fetch(`${SERVER_URL}/api/datasets/search?q=${searchQuery}&page=${currentPage}` , {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ filters: currentFilters }),
        });
        const data = await response.json();
        setFilters(data.facets);
        setResults(data.results);
        setPagination(data.pagination);
      } catch (error) {
        console.error('Error fetching data:', error);
      }
    };   

    const handleApplyFilters = async () => {
      fetchData(pagination.page, activeFilters, globalSearchTerm);

      // Update the list of applied filters
      const filtersList = Object.entries(activeFilters).map(([category, values]) => {
        return values.map(value => ({ category, value }));
      }).flat();
      setAppliedFilters(filtersList);
    };

    useEffect(() => {   
      fetchData(pagination.page, activeFilters, globalSearchTerm);
    }, []); // Refetch when activeFilters change

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
      fetchData(newPage, activeFilters, globalSearchTerm);
    };

    const handleSearchSubmit = (event) => {
      event.preventDefault();
      fetchData(pagination.page, activeFilters, globalSearchTerm);
      console.log("Search Handled");
    };

    const handleRemoveFilter = (category, value) => {
      setActiveFilters(prevFilters => {
        const newFilters = { ...prevFilters };
        newFilters[category] = newFilters[category].filter(v => v !== value);
  
        if (newFilters[category].length === 0) {
          delete newFilters[category];
        }
  
        return newFilters;
      });
  
      // Remove filter from the list of applied filters
      setAppliedFilters(prevFilters => prevFilters.filter(filter => !(filter.category === category && filter.value === value)));

      fetchData(pagination.page, activeFilters, globalSearchTerm);

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
                              onApplyFilters={handleApplyFilters}
                          />
                      ))}
                      <div className='filters-toggle-div'>
                        {Object.keys(filters).length > 4 && (
                            <button onClick={toggleMoreFacets} className='filters-toggle'>
                                <FontAwesomeIcon icon={faSliders} /> <p>{showMoreFacets ? 'Less facets' : 'More facets'}</p>
                            </button>
                        )}
                      </div>
                  </div>
                  <div className='study-keyword-search'>
                    <span class="text-search search-title">Search by text <FontAwesomeIcon icon={faQuestionCircle} /></span>
                    <div>
                      <form onSubmit={handleSearchSubmit}>
                        <input
                            type="text"
                            autoComplete="off"
                            className="w-full dark:bg-gray-950 pl-8 form-input-alt h-9 pr-3 focus:shadow-xl"
                            placeholder="Search..."
                            value={globalSearchTerm}
                            onChange={(e) => setGlobalSearchTerm(e.target.value)}
                        />
                        
                        {/* <svg className="absolute left-2.5 text-gray-400 top-1/2 transform -translate-y-1/2" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 32 32">
                            <path d="M30 28.59L22.45 21A11 11 0 1 0 21 22.45L28.59 30zM5 14a9 9 0 1 1 9 9a9 9 0 0 1-9-9z" fill="currentColor"></path>
                        </svg>     */}

                      </form>
                    </div>
                  </div>

                </div>
                <div className='applied-filters-container'>
                  {appliedFilters.length > 0 && (
                          <div className="applied-filters">
                            <p>Applied Filters:</p>
                            {appliedFilters.map((filter, index) => (
                              <div key={index} className="applied-filter">
                                {filter.category}: {filter.value}
                                <span  className="cross-icon" onClick={() => handleRemoveFilter(filter.category, filter.value)}>&times;</span>
                              </div>
                            ))}
                          </div>
                  )}
                </div>

                <div className='table-pagination'>
                  <Pagination
                    pagination={pagination}
                    onPageChange={onPageChange}
                  />
                </div>
                
                <div className='total-results-count'>
                  <p>{pagination.totalCount} results found!</p>
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
