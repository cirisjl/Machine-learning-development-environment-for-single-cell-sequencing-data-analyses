import React, { useState, useEffect } from 'react';
import CreatableSelect from 'react-select/creatable';

const createOption = (label) => ({
  label,
  value: label.toLowerCase().replace(/\W/g, ''),
});

function MyCreatableSelect({ fieldName, options, formData, setFormData}) {

  const [isLoading, setIsLoading] = useState(false);
  const [value, setValue] = useState(null);

  const filteredOptions = options[fieldName] || [];

  console.log(filteredOptions);

  // Initialize the options 
  const [filteredOptionsState, setFilteredOptionsState] = useState(
    filteredOptions.map((option) => createOption(option))
  );

  const handleChange = (inputValue) => {
    setValue(inputValue);

    // Update the formData state with the selected option
    setFormData({
      ...formData,
      [fieldName]: inputValue,
    });
  }

  const handleCreateOption = (inputValue) => {
    setIsLoading(true);

    setTimeout(() => {
      const newOption = createOption(inputValue);
      setIsLoading(false);
      setFilteredOptionsState((prev) => [...prev, newOption]);
      setValue(newOption);
    }, 1000);

    // Update the formData state with the newly created option
    setFormData({
      ...formData,
      [fieldName]: inputValue,
    });
  };

  return (
    <CreatableSelect
      isClearable
      isDisabled={isLoading}
      isSearchable
      onChange={handleChange}
      // onChange={(newValue) => setValue(newValue)}
      onCreateOption={handleCreateOption}
      options={filteredOptionsState}
      value={value}
    />
  );
}

export default MyCreatableSelect;
