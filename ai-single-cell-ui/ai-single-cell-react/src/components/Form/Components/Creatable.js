import React, { useState, useEffect } from 'react';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';


function MyCreatableSelect({ fieldName }) {
  const [selectedOption, setSelectedOption] = useState(null);
 // Initialize the options with some initial values
 const [options, setOptions] = useState([
    { value: 'option1', label: 'Option 1' },
    { value: 'option2', label: 'Option 2' },
  ]);

  useEffect(() => {
    // Fetch options from MongoDB based on the fieldName
    const fetchOptions = async () => {
      try {
        let optionsURL = `${SERVER_URL}/mongoDB/api/options`
        const response = await fetch(optionsURL);
        const data = await response.json();
        setOptions(data.options.map(option => ({ value: option, label: option })));
      } catch (error) {
        console.error(`Error fetching options for ${fieldName}:`, error);
      }
    };

    fetchOptions();
  }, [fieldName]);

  const handleChange = (newValue) => {
    setSelectedOption(newValue);
  };

  const handleCreateOption = (inputValue) => {
    // Create a new option object
    const newOption = { value: inputValue, label: inputValue };

    // Update the options state to include the new option
    setOptions([...options, newOption]);

    // Set the selected option to the newly created option
    setSelectedOption(newOption);
  };

  return (
    <CreatableSelect
      isClearable
      isSearchable
      onChange={handleChange}
      onCreateOption={handleCreateOption}
      options={options}
      value={selectedOption}
    />
  );
}

export default MyCreatableSelect;
