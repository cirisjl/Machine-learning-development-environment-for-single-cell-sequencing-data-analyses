import React, { useState, useEffect } from 'react';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';

function MyCreatableSelect({ fieldName, options }) {
  console.log("mycreatbale component");
  console.log(fieldName);
  console.log(options);
  const [selectedOption, setSelectedOption] = useState(null);

  const filteredOptions = options[fieldName] || [];

  console.log(filteredOptions);

  // Initialize the options 
  const [filteredOptionsState, setFilteredOptionsState] = useState(filteredOptions);

  const handleChange = (newValue) => {
    setSelectedOption(newValue);
  };

  const handleCreateOption = (inputValue) => {
    setSelectedOption(inputValue);
    // Update the options state to include the new option
    setFilteredOptionsState([...filteredOptionsState, inputValue]);

    // Set the selected option to the newly created option
    setSelectedOption(inputValue);
  };

  return (
    <CreatableSelect
      isClearable
      isSearchable
      onChange={handleChange}
      onCreateOption={handleCreateOption}
      options={filteredOptionsState.map((option) => ({ value: option, label: option }))}
      value={selectedOption}
    />
  );
}

export default MyCreatableSelect;
