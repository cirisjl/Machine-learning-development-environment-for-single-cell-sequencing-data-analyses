import React, { useState, useEffect } from 'react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';

function LeftNav(props) {
  const [categories, setCategories] = useState([]);

  const TOOLS_FILTERS = "http://localhost:3001";

  const handleFilterSelection = props.handleFilterSelection


  useEffect(() => {
    fetch(TOOLS_FILTERS + "/api/tools/leftnav")
      .then(response => response.json())
      .then(data => {
        const updatedCategories = data.map((category, index) => {
          return {
            ...category,
            expanded: index < 2 // Set the first 4 categories to be expanded by default
          };
        });
        setCategories(updatedCategories);
      })
      .catch(error => console.log(error));
  }, []);

  const toggleCategory = categoryId => {
    console.log("Filter clicked");
    setCategories(
      categories.map(category =>
        category.category_id === categoryId
          ? { ...category, expanded: !category.expanded }
          : category
      )
    );
  };

  return (
    <nav>
      <div className='tool-panel-section'>
      {categories.map(category => (
        <div key={category.category_id} className='category-level'>
          <h3 className='category-level-header' onClick={() => toggleCategory(category.category_id)}>
            {category.category_name}
            <span className="category-icon">
            <FontAwesomeIcon
              icon={category.expanded ? faAngleDown : faAngleRight}
            />
            </span>
            </h3>
          {category.expanded && (
            <ul
            className={`category-filters filter-level-ul ${
              category.expanded ? 'expanded' : 'collapsed'
            }`}
          >
            {category.filters.map(filter => (
              <li key={filter} className='filter-level-li' onClick={() => handleFilterSelection(category.category_name.toLowerCase().replace(/\s/g, '_') + "_" + filter.toLowerCase())}>{filter}</li>
            ))}
          </ul>
            )}
        </div>
      ))}
      </div>
    </nav>
  );
}

export default LeftNav;
