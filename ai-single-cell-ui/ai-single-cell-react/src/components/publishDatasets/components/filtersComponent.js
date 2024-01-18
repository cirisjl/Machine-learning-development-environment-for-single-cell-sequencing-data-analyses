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
                                    <form>
                                        <input
                                            type="text"
                                            autoComplete="off"
                                            className="w-full dark:bg-gray-950 pl-8 form-input-alt h-9 pr-3 focus:shadow-xl"
                                            // placeholder={placeHolder}
                                            // spellCheck="false"
                                            // value={searchTerm}
                                            // onChange={handleInputChange}
                                        />
                                        
                                        {/* <svg className="absolute left-2.5 text-gray-400 top-1/2 transform -translate-y-1/2" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 32 32">
                                            <path d="M30 28.59L22.45 21A11 11 0 1 0 21 22.45L28.59 30zM5 14a9 9 0 1 1 9 9a9 9 0 0 1-9-9z" fill="currentColor"></path>
                                        </svg>     */}

                                        </form>
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
