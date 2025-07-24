import React, { useState } from 'react';
import CreatableSelect from 'react-select';

const SelectComponent = ({value, onChange,options}) => {


    const selectOptions = options.opts.map(option => ({
        label: option, // Use the string value as both label and value
        value: option
      }));

  const handleChange = selectedOption => {
    onChange(selectedOption ? selectedOption.value : null);
  };

  return (
    <div>
      <CreatableSelect 
        value={selectOptions.find(option => option.value === value) || null}
        onChange={handleChange}
        options={selectOptions}
        isClearable={options.clearable}
        placeholder={options.placeholder}
        isCreatable={options.creatable} 
        isSearchable={options.searchable}
        />
    </div>
  );
};

export default SelectComponent;
