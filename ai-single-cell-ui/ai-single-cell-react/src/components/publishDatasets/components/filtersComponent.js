import React, {useState} from 'react';


const FilterComponent = ({ name, options, activeFilters, onFilterChange, isVisible, onCategoryChange  }) => {

    // const [isVisible, setIsVisible] = useState(false);

    // const toggleVisibility = () => setIsVisible(!isVisible);

    const isActive = (filterValue) => {
        return activeFilters[name] && activeFilters[name].includes(filterValue);
    };

    return (
        <div>
            <h3 onClick={() => onCategoryChange(name)}>{name}</h3>
            {isVisible && (
                <div>
                    <ul>
                        {options.map(option => (
                            <li key={option.value}>
                                <label>
                                    <input
                                        type="checkbox"
                                        checked={activeFilters[name]?.includes(option.value)}
                                        onChange={() => onFilterChange(name, option.value)}
                                    />
                                    {option.label}
                                </label>
                            </li>
                        ))}
                    </ul>
                    <button className="apply-filters-button">Apply</button>
                </div>
            )}
        </div>
    );

  };


  export default FilterComponent;
