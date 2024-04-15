import React from 'react';
import Box from '@mui/material/Box';
import Slider from '@mui/material/Slider';

// Custom widget component
const GeneRangeSlider = ({value, onChange}) => {
  const handleChange = (event, newValue) => {
    onChange(newValue); // newValue will be an array [min, max]
  };

  // // This function is optional, only necessary if you want to handle the commit of the change
  // const handleDragChange = (event, newValue) => {
  //   onChange(newValue); // newValue will be an array [min, max]
  // };

  // Ensure we have a default value
  const defaultValue = value || [200, 20000];

  return (
    <Box sx={{ m: 2 }}>
      <Slider
        value={defaultValue}
        onChange={handleChange}
        // onChangeCommitted={handleDragChange} // Use onChangeCommitted for smoother drag behavior
        valueLabelDisplay="auto"
        min={0}
        max={20000}
        step={25}
        marks={[
          { value: 200, label: '200*' },
          { value: 1000, label: '1000' },
          { value: 5000, label: '5000' },
          { value: 10000, label: '10000' },
          { value: 15000, label: '15000' },
          { value: 20000, label: '20000' }
        ]}
      />
    </Box>
  );
};

export default GeneRangeSlider;
