import React, {useState} from 'react';


const FilterComponent = ({ name, options, activeFilters, onFilterChange, isVisible, onCategoryChange  }) => {

    // const [isVisible, setIsVisible] = useState(false);

    // const toggleVisibility = () => setIsVisible(!isVisible);

    const isActive = (filterValue) => {
        return activeFilters[name] && activeFilters[name].includes(filterValue);
    };

    return (
        <div className='facet'>
            <div className='filter-category'>
                <p onClick={() => onCategoryChange(name)}>{name}</p>
            </div>
            {isVisible && (
                <div className='filters-box-searchable'>
                    <div className='facet-container'>
                        <div className='single-facet'>
                            <div>
                                <div className='filters-searchbar'>
                                    <input autocomplete="false" placeholder="Search for a filter" type="text" id="filters-search-bar-filters-box-searchable-organ" class="form-control" value="" fdprocessedid="vdtrz3"/>
                                </div>

                                <div className='filters-header'>
                                    <h4>Available filters:</h4>
                                </div>

                                <div className='filters-options'>
                                    <ul>
                                        {options.map(option => (
                                            <li key={option._id}>
                                                <label>
                                                    <input
                                                        type="checkbox"
                                                        checked={activeFilters[name]?.includes(option._id)}
                                                        onChange={() => onFilterChange(name, option._id)}
                                                    />
                                                    {option._id} ({option.count})
                                                </label>
                                            </li>
                                        ))}
                                    </ul>
                                </div>
                                
                            </div>
                        </div>
                    </div>

                    <div className='filters-footer'>
                        <button className="apply-filters-button">Apply</button>
                    </div>
                </div>
            )}
        </div>
    );

  };


  export default FilterComponent;
