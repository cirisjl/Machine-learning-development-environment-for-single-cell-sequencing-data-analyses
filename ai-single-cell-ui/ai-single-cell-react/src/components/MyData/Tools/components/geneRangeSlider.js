import React from 'react';
import Box from '@mui/material/Box';
import Typography from '@mui/material/Typography';
import Slider from '@mui/material/Slider';

// Custom widget component
const GeneRangeSlider = ({ value, onChange }) => {
  const handleChange = (event, newValue) => {
    onChange(newValue); // newValue will be an array [min, max]
  };

  const handleDragChange = (event, newValue) => {
    onChange(newValue); // newValue will be an array [min, max]
  };

  return (
    <Box sx={{ m: 2 }}>
      <Typography gutterBottom>Min Genes - Max Genes: <b>[{value[0]} - {value[1]}]</b></Typography>
      <Slider
        value={value}
        onChange={handleChange}
        onChangeCommitted={handleDragChange} // Use onChangeCommitted for smoother drag behavior
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
