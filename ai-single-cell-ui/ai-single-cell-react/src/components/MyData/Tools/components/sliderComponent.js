import React from 'react';
import Typography from '@mui/material/Typography';
import Slider from '@mui/material/Slider';

const RangeSlider = ({ title, value, onChange,options }) => {
  const handleChange = (event, newValue) => {
    onChange(newValue);
  };

  const handleDragChange = (event, newValue) => {
    onChange(newValue); // newValue will be an array [min, max]
  };

  return (
    <div>
      {/* <Typography gutterBottom>{options.title}: <b>{value}</b></Typography> */}
      <Slider
        value={value}
        onChange={handleChange}
        onChangeCommitted={handleDragChange} // Use onChangeCommitted for smoother drag behavior
        valueLabelDisplay="auto"
        min={options.min}
        max={options.max}
        step={options.step}
        marks={options.marks}
      />
    </div>
  );
};

export default RangeSlider;
