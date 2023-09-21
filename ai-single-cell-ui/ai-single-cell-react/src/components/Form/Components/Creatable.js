import React, { useState, useEffect } from 'react';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';


function MyCreatableSelect({ fieldName, options }) {
  const [selectedOption, setSelectedOption] = useState(null);

  const filteredOptions = options[fieldName] || [];

 // Initialize the options 
 const [options, setOptions] = useState(filteredOptions);

  const handleChange = (newValue) => {
    setSelectedOption(newValue);
  };

  const handleCreateOption = (inputValue) => {

    // Update the options state to include the new option
    setOptions([...options, inputValue]);

    // Set the selected option to the newly created option
    setSelectedOption(inputValue);
  };

  return (
    <CreatableSelect
      isClearable
      isSearchable
      onChange={handleChange}
      onCreateOption={handleCreateOption}
      options={options.map((option) => ({ value: option, label: option }))}
      value={selectedOption}
    />
  );
}

export default MyCreatableSelect;
