import React, { useState,useEffect, useMemo } from 'react';
import FilterComponent from './filtersComponent';
import {SERVER_URL} from '../../../constants/declarations';
import ResultsTable from './tableResultsComponent';
import Pagination from './tablePaginationComponent';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import {  faQuestionCircle, faSliders } from '@fortawesome/free-solid-svg-icons';
import SearchBox from '../../Header/searchBar';
import CreatableSelect from 'react-select/creatable';


const DatasetSelectionDialog = ({onSelect, multiple, onClose , isVisible, taskData }) => {

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };

    const dialogStyle1= {
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

    const [appliedFilters, setAppliedFilters] = useState([]);
    const [showview,setShowView]= useState(false)

    const handleShowView=(x)=>{
      setShowView(x);
    }

    const handlecloseView=()=>{
      setShowView(false);
    }

    const onSelectDataset = (dataset) => {
      // Get a copy of the current selected datasets from the state or props
      const currentSelectedDatasets = { ...taskData.task_builder.selectedDatasets };
      const datasetId = dataset.Id; // Make sure 'Id' is the correct field for dataset ID
    
      if (currentSelectedDatasets[datasetId]) {
          // Dataset is currently selected, deselect it
          delete currentSelectedDatasets[datasetId];
      } else {
          // Dataset is not selected, select it
          currentSelectedDatasets[datasetId] = dataset;
      }
    
      // Call onSelect with the updated selected datasets
      onSelect(currentSelectedDatasets);
    };

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
      setShowMoreFacets(false);
      fetchData(1, activeFilters, globalSearchTerm);

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
  
  useEffect(()=>{
console.log(taskData);
  })


  const toggleMoreFacets = () => {
    setShowMoreFacets(!showMoreFacets);
    // Show all facets if showMoreFacets is true, else show only the first four
    if (!showMoreFacets) {
        setVisibleFacets(Object.keys(filters));
    } else {
        setVisibleFacets(Object.keys(filters).slice(0, 4));
    }
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
      setShowMoreFacets(false);
      fetchData(newPage, activeFilters, globalSearchTerm);
    };

    const handleSearchSubmit = (event) => {
      event.preventDefault();
      setShowMoreFacets(false);
      fetchData(1 , activeFilters, globalSearchTerm);
      console.log("Search Handled");
    };

    const handleRemoveFilter = (category, value) => {
      setShowMoreFacets(false);
      setActiveFilters(prevFilters => {
        const newFilters = { ...prevFilters };
        newFilters[category] = newFilters[category].filter(v => v !== value);
  
        if (newFilters[category].length === 0) {
          delete newFilters[category];
        }
  
        fetchData(1, newFilters, globalSearchTerm);
        return newFilters;
      });
  
      // Remove filter from the list of applied filters
      setAppliedFilters(prevFilters => prevFilters.filter(filter => !(filter.category === category && filter.value === value)));

    };

    return (
      <div>
        <div style={dialogStyle} className="dialog-container">
            <div className="dialog-backdrop" onClick={onClose} />
            <div className="dialog">
              <div>
                
                <div className='filters-and-search-container'>
                  <div className='metadata-search-wrap filters-container'>
                  <span className="metadata-search search-title">Search by filters--- <FontAwesomeIcon icon={faQuestionCircle}/></span>
                      {visibleFacets.map((filterName) => (
                          <FilterComponent
                              key={filterName}
                              name={filterName}
                              options={filters[filterName]}
                              activeFilters={activeFilters}
                              onFilterChange={handleFilterChange}
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
                    <span className="text-search search-title">Search by text <FontAwesomeIcon icon={faQuestionCircle} /></span>
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
                
                <div className='table-results'>
                     <ResultsTable data={results} onSelectDataset={onSelectDataset} selectedDatasets={taskData.task_builder.selectedDatasets} multiple={multiple} pagination={pagination}  />
                </div>
                <div className='dialog-close'>
                    <button onClick={onClose}>Close</button>
                </div>
              </div>
              
            </div>
        </div>
        {/*--------------------------------– new form here *------------------------------------------*/}
        {false && (
          <div style={dialogStyle1} className="dialog-container">
          <div className="dialog">

          <div className="my-form-container">
        
        <div>
        <h2 className="form-title">My Form</h2>
        <form  className="form">
          {/* Dataset */}
          <div className="form-field">
            <label className="form-label">Dataset:</label>
            <input
              type="text"
              name="Dataset"
            //   required
            />
          </div>

          

          {/* Downloads */}
          <div className="form-field">
            <label className="form-label">Downloads:</label>
            <input
              type="text"
              required
              name="Downloads"
             
            />

          </div>

          <div className="form-field">
            <label className="form-label">Title:</label>
            <input
              type="text"
              name="Title"
              
              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Author:</label>
            <input
              type="text"
              name="Author"
              required
              
            />
          </div>

          <div className="form-field">
            <label className="form-label">Reference (paper):</label>
            <input
              type="text"
              name="Reference (paper)"
              
              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Abstract:</label>
            <textarea
              name="Abstract"
           
              className="form-input"
            />
          </div>

          {/* DOI */}
          <div className="form-field">
            <label className="form-label">DOI:</label>
            <input
              type="text"
              name="DOI"
              
              placeholder="http://"
              className="form-input"
            />
          </div>


          {/* Species (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Species:</label>
            <CreatableSelect
              name="Species"
            //   value={formData.Species}
              isClearable
              isSearchable
              required
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Species', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Species', inputValue)}
            //   options={options.Species} // Set options to the fetched options
            />
          </div>

          {/* "Sample Type" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Sample Type:</label>
            <CreatableSelect
              name="Sample Type"
            //   value={formData['Sample Type']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Sample Type', selectedOption)} // Use handleSelectChange             
               onCreateOption={(inputValue) => this.handleCreateOption('Sample Type', inputValue)}
            //   options={options['Sample Type']} // Set options to the fetched options
            />
          </div>


          {/* "Anatomical Entity" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Anatomical Entity:</label>
            <CreatableSelect
              name="Anatomical Entity"
            //   value={formData['Anatomical Entity']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Anatomical Entity', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Anatomical Entity', inputValue)}
            //   options={options['Anatomical Entity']} // Set options to the fetched options
            />
          </div>

          {/* "Organ Part" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Organ Part:</label>
            <CreatableSelect
              name="Organ Part"
            //   value={formData['Organ Part']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Organ Part', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Organ Part', inputValue)}
            //   options={options['Organ Part']} // Set options to the fetched options
            />
]                  </div>

          {/* "Model Organ" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Model Organ:</label>
            <CreatableSelect
              name="Model Organ"
            //   value={formData['Model Organ']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Model Organ', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Model Organ', inputValue)}
            //   options={options['Model Organ']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          {/* "Selected Cell Types" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Selected Cell Types:</label>
            <CreatableSelect
              name="Selected Cell Types"
            //   value={formData['Selected Cell Types']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Selected Cell Types', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Selected Cell Types', inputValue)}
            //   options={options['Selected Cell Types']} // Set options to the fetched options
            />
          </div>



          {/* "Library Construction Method" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Library Construction Method:</label>
            <CreatableSelect
              name="Library Construction Method"
            //   value={formData['Library Construction Method']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Library Construction Method', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Library Construction Method', inputValue)}
            //   options={options['Library Construction Method']} // Set options to the fetched options
              className="form-input"
            />
          </div>


          {/* "Nucleic Acid Source" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Nucleic Acid Source:</label>
            <CreatableSelect
              name="Nucleic Acid Source"
            //   value={formData['Nucleic Acid Source']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Nucleic Acid Source', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Nucleic Acid Source', inputValue)}
            //   options={options['Nucleic Acid Source']} // Set options to the fetched options
              className="form-input"
            />
          </div>


          <div className="form-field">
            <label className="form-label">Paired End:</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="Paired End"
                  value= "true"
                //   checked={formData["Paired End"] === "true"}
                //   onChange={this.handleChange}
                  className="form-input"
                />
                True
              </label>
              <label className="form-label">
                <input
                  type="radio"
                  name="Paired End"
                  value="false"
                //   checked={formData["Paired End"] === "false"}
                //   onChange={this.handleChange}
                  className="form-input"
                />
                False
              </label>
            </div>
          </div>

          <div className="form-field">
            <label className="form-label">Analysis Protocol:</label>
            <input
              type="text"
              name="Analysis Protocol"
            //   value={formData['Analysis Protocol']}
            //   onChange={this.handleChange}
              className="form-input"
            />
          </div>

          {/* "Disease Status (Specimen)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Specimen):</label>
            <CreatableSelect
              name="Disease Status (Specimen)"
            //   value={formData['Disease Status (Specimen)']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Disease Status (Specimen)', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Specimen)', inputValue)}
            //   options={options['Disease Status (Specimen)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Specimen)'] && <div className="error-tooltip">{errors['Disease Status (Specimen)']}</div>} */}
          </div>


          {/* "Disease Status (Donor)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Donor):</label>
            <CreatableSelect
              name="Disease Status (Donor)"
            //   value={formData['Disease Status (Donor)']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Disease Status (Donor)', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Donor)', inputValue)}
            //   options={options['Disease Status (Donor)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Donor)'] && <div className="error-tooltip">{errors['Disease Status (Donor)']}</div>} */}
          </div>

          {/* "Development Stage" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Development Stage:</label>
            <CreatableSelect
              name="Development Stage"
            //   value={formData['Development Stage']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Development Stage', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Development Stage', inputValue)}
            //   options={options['Development Stage']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Donor Count:</label>
            <input
              type="number"
              name="Donor Count"
            //   value={formData["Donor Count"]}
            //   onChange={this.handleChange}
              className="form-input"
            />
          </div>


         

          {/* "Source" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Source:</label>
            <CreatableSelect
              name="Source"
            //   value={formData['Source']}
            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Source', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Source', inputValue)}
            //   options={options['Source']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          
          {/* Source Key */}
          <div className="form-field">
            <label className="form-label">Source Key:</label>
            <input
              type="text"
              name="Source Key"
            //   value={formData['Source Key']}
            //   onChange={this.handleChange}
              placeholder="Enter ..."
              className="form-input"
            />
            {/* {errors['Source Key'] && <p className="error">{errors['Source Key']}</p>} */}
          </div>

          <div className="form-field">
            <label>Submission Date:</label>
            <input
              type="date"
              required
              name="Submission Date"
            //   value={formData["Submission Date"]}
            //   onChange={this.handleChange}
            />
            {/* {errors['Submission Date'] && <div className="error-tooltip">{errors['Submission Date']}</div>} */}
          </div>

          <div>
          <div className='dialog-close'>
                    <button onClick={handlecloseView} >Close</button>
                </div>
                
          </div>

          

        </form>
        </div>
      </div>
           

           
            </div>
          </div>
        )} 

        </div>

    );
  };

export default DatasetSelectionDialog;
