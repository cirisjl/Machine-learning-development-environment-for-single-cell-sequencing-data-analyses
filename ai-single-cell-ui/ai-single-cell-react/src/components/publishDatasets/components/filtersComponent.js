import React from 'react';


const FilterComponent = ({ name, options, activeFilters, onFilterChange }) => {

    const [isVisible, setIsVisible] = useState(false);

    const toggleVisibility = () => setIsVisible(!isVisible);

    const isActive = (filterValue) => {
        return activeFilters[name] && activeFilters[name].includes(filterValue);
    };

    return (
        <div>
            <h3 onClick={toggleVisibility}>{name}</h3>
            {isVisible && (
                <ul>
                    {options.map(option => (
                        <li key={option._id}>
                            <button
                                type="button"
                                onClick={() => onFilterChange(name, option._id)}
                                className={isActive(option._id) ? 'active' : ''}
                            >
                                {option._id} ({option.count})
                            </button>
                        </li>
                    ))}
                </ul>
            )}
        </div>
    );

  };


  export default FilterComponent;
