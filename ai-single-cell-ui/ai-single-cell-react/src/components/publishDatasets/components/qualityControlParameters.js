import React, { useState } from 'react';
import Box from '@material-ui/core/Box';
import Typography from '@material-ui/core/Typography';
import Slider from '@mui/material/Slider';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import {ExpansionPanel, ExpansionPanelSummary, ExpansionPanelDetails, Button } from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import { makeStyles } from '@material-ui/core/styles';
import Switch from "react-switch";

const QualityControlParameters = ({values, setValues, defaultValues, shouldHideForSeurat}) => {

  const handleSliderChange = (event) => {
    const { name, value } = event.target;
    setValues((prevValues) => ({
      ...prevValues,
      [name]: value,
      use_default: false,
    }));
  };

  const useStyles = makeStyles((theme) => ({
    root: {
      width: '100%',
      marginBottom: theme.spacing(2), // Adds space between expansion panels
    },
    panelSummary: {
      '& .MuiTypography-h6': {
        fontSize: '15px', // Targeting the content directly
      },
      '& .MuiTypography-body1': {
        fontSize: '15px', // Ensuring all typography inside summary also gets the font size
      }
    },
    panelDetails: {
      display: 'block', // Ensures the contents take the full width and are blocked
      '& .MuiTypography-h6': {
        fontSize: '15px', // Targeting the content directly
      },
      '& .MuiTypography-body1': {
        fontSize: '15px', // Ensuring all typography inside summary also gets the font size
      }
    },
    valueLabel: {
      // Increase the size of the value label to fit longer text
      '& .MuiSlider-valueLabel': {
        width: 'auto', // Allow the width to auto-adjust to content
        minWidth: '40px', // Ensure a minimum width to avoid too narrow labels
        borderRadius: '4px', // Adjust for a more rectangular shape
        padding: '0 8px', // Add some padding horizontally to accommodate wider numbers
  
        // Ensure the font size is 14px within the value label
        '& .MuiTypography-root': {
          fontSize: '14px',
        },
      },
    },
    customSwitch: {
      '& .MuiSwitch-switchBase.Mui-checked': {
        color: '#1976d2',
        '&:hover': {
          backgroundColor: 'rgba(25, 118, 210, 0.04)', // Lighter shade for hover, adjust opacity as needed
        },
      },
      '& .MuiSwitch-switchBase.Mui-checked + .MuiSwitch-track': {
        backgroundColor: '#1976d2',
      },
    },
  }));

  const classes = useStyles();


  const handleSwitchChange = (name) => (checked) => {
    setValues((prevValues) => ({
      ...prevValues,
      [name]: checked,
      ...(name === 'use_default' && checked ? defaultValues : {}),
    }));
  };

  const resetToDefaults = () => {
    setValues(defaultValues);
  };

  // Constructing dual-slider for min_genes and max_genes as they share the same slider
  const handleGeneSliderChange = (event, newValue) => {
    setValues((prevValues) => ({
      ...prevValues,
      min_genes: newValue[0],
      max_genes: newValue[1],
      use_default: false,
    }));
  };

  return (
    <div>
     <Typography variant="h6" gutterBottom style={{ fontWeight: 'bold' }}>
        Advanced Quality Control Parameters
      </Typography>
      
      <ExpansionPanel className={classes.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={classes.panelSummary}>
            {/* QC Parameters */}
            <Typography variant="h6" gutterBottom>QC Parameters</Typography>
        </ExpansionPanelSummary>
        <ExpansionPanelDetails className={classes.panelDetails}>
             <FormGroup>
                <Box sx={{ m: 2 }}>
                  <Typography variant="caption" display="block" gutterBottom>
                    <b>Note:</b> An asterisk (*) indicates the default value.
                  </Typography>
                </Box>

                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>Min Genes - Max Genes: <b>[{values.min_genes} - {values.max_genes}]</b></Typography>
                    <Slider
                    value={[values.min_genes, values.max_genes]}
                    onChange={handleGeneSliderChange}
                    valueLabelDisplay="auto"
                    min={0}
                    max={20000}
                    step={25}
                    name="min_max_genes"
                    marks={[
                      { value: 200, label: '200*' },
                      { value: 1000, label: '1000' },
                      { value: 5000, label: '5000' },
                      { value: 10000, label: '10000' },
                      { value: 15000, label: '15000' },
                      { value: 20000, label: '20000' },
                    ]}
                    />
                </Box>

                <Box sx={{ m: 2 }}>

                    <Typography gutterBottom>Min Cells: <b>{values.min_cells}</b></Typography>
                    <Slider
                    value={values.min_cells}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'min_cells', value: val } })}
                    valueLabelDisplay="auto"
                    min={1}
                    max={200}
                    step={1}
                    name="min_cells"
                    marks={[
                      { value: 2, label: '2*' },
                      { value: 10, label: '10' },
                      { value: 50, label: '50' },
                      { value: 100, label: '100' },
                      { value: 200, label: '200' },
                    ]}
                    />
                </Box>
              
                {!shouldHideForSeurat && (
                  <Box sx={{ m: 2 }}>
                      <Typography gutterBottom>Target Sum: <b>{values.target_sum.toExponential()}</b></Typography>
                      <Slider
                      value={values.target_sum}
                      onChange={(e, val) => handleSliderChange({ target: { name: 'target_sum', value: val } })}
                      //   valueLabelDisplay="auto"
                      min={0}
                      max={1e6}
                      step={1e4}
                      name="target_sum"
                      marks={[
                        // { value: 0, label: '0' },
                        { value: 1e4, label: '1e4*' },
                        { value: 1e5, label: '1e5' },
                        { value: 1e6, label: '1e6' },
                      ]}
                      />
                  </Box>
                )}

                <Box sx={{ m: 2 }}>

                    <Typography gutterBottom>Highly Variable Genes (n_top_genes): <b>{values.n_top_genes}</b></Typography>
                    <Slider
                    value={values.n_top_genes}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'n_top_genes', value: val } })}
                    valueLabelDisplay="auto"
                    min={100}
                    max={10000}
                    step={25}
                    name="n_top_genes"
                    marks={[
                      { value: 100, label: '100' },
                      { value: 500, label: '500' },
                      { value: 1000, label: '1000' },
                      { value: 2000, label: '2000*' },
                      { value: 5000, label: '5000' },
                      { value: 10000, label: '10000'}]}
                    />
                </Box>

                  <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>
                      Expected Doublet Rate: <b>{`${(values.doublet_rate * 100).toFixed(2)}%`}</b>
                    </Typography>
                    <Typography variant="caption" display="block" gutterBottom>
                      <b>Note:</b> A rate of <b>0%</b> means not to classify doublets.
                    </Typography>
                    <Slider
                      value={values.doublet_rate}
                      min={0}
                      max={0.5}
                      step={0.001}
                      onChange={(e, val) => handleSliderChange({ target: { name: 'doublet_rate', value: val } })}
                      valueLabelDisplay="auto"
                      valueLabelFormat={(value) => `${(value * 100).toFixed(2)}%`}
                      marks={[
                        { value: 0, label: '0%*' }, // Default
                        // { value: 0.008, label: '0.8%' },
                        // { value: 0.023, label: '2.3%' },
                        // { value: 0.038, label: '3.8%' },
                        // { value: 0.046, label: '4.6%' },
                        // { value: 0.061, label: '6.1%' },
                        { value: 0.08, label: '8%' }, 
                        { value: 0.125, label: '12.5%' },
                        { value: 0.2, label: '20%' },
                        { value: 0.5, label: '50%' },
                      ]}
                      name="doublet_rate"
                      className={classes.valueLabel}
                    />
                  </Box>

                <Box sx={{ m: 2 }}>
                    <div>
                      <label htmlFor="material-switch">
                        <p>{`Regress Cell Cycle: ${values.regress_cell_cycle ? 'Yes' : 'No'}`}</p>
                        <Switch
                          checked={values.regress_cell_cycle}
                          onChange={handleSwitchChange("regress_cell_cycle")}
                          onColor="#86d3ff"
                          onHandleColor="#2693e6"
                          handleDiameter={30}
                          uncheckedIcon={false}
                          checkedIcon={false}
                          boxShadow="0px 1px 5px rgba(0, 0, 0, 0.6)"
                          activeBoxShadow="0px 0px 1px 10px rgba(0, 0, 0, 0.2)"
                          height={20}
                          width={48}
                          className="react-switch"
                          id="regress_cell_cycle"
                        />
                      </label>
                    </div>
                </Box>
            </FormGroup>
            </ExpansionPanelDetails>
      </ExpansionPanel>


      <ExpansionPanel className={classes.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={classes.panelSummary}>
            {/* Projection Parameters */}
            <Typography variant="h6" gutterBottom>Projection Parameters</Typography>
        </ExpansionPanelSummary>
        
        <ExpansionPanelDetails className={classes.panelDetails}>
            <FormGroup>
                <Box sx={{ m: 2 }}>
                  <Typography variant="caption" display="block" gutterBottom>
                  <b>Note:</b> An asterisk (*) indicates the default value.
                  </Typography>
                </Box>
                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>n_neighbors: <b>{values.n_neighbors}</b></Typography>
                    <Slider
                    value={values.n_neighbors}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'n_neighbors', value: val } })}
                    valueLabelDisplay="auto"
                    min={2}
                    max={100}
                    step={1}
                    name="n_neighbors"
                    marks={[
                      { value: 2, label: '2' },
                      { value: 5, label: '5' },
                      { value: 10, label: '10' },
                      { value: 15, label: '15*' },
                      { value: 20, label: '20' },
                      { value: 50, label: '50' },
                      { value: 100, label: '100' },
                    ]}
                    />
                </Box>

                {!shouldHideForSeurat && (
                  <Box sx={{ m: 2 }}>
                      <Typography gutterBottom>n_pcs: <b>{values.n_pcs}</b></Typography>
                      <Slider
                      value={values.n_pcs}
                      onChange={(e, val) => handleSliderChange({ target: { name: 'n_pcs', value: val } })}
                      valueLabelDisplay="auto"
                      min={0}
                      max={200}
                      step={1}
                      name="n_pcs"
                      marks={[
                        { value: 0, label: '0*' },
                        { value: 5, label: '5' },
                        { value: 10, label: '10' },
                        { value: 20, label: '20' },
                        { value: 40, label: '40' },
                        { value: 50, label: '50' },
                        { value: 125, label: '125' },
                        { value: 200, label: '200' },
                      ]}
                      />
                  </Box>
                )}
            </FormGroup>
        </ExpansionPanelDetails>
      </ExpansionPanel> 


      <ExpansionPanel className={classes.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={classes.panelSummary}>
          <Typography variant="h6" gutterBottom>Clustering Parameters</Typography>
        </ExpansionPanelSummary>

        <ExpansionPanelDetails className={classes.panelDetails}>
            <FormGroup>
                <Box sx={{ m: 2 }}>
                  <Typography variant="caption" display="block" gutterBottom>
                  <b>Note:</b> An asterisk (*) indicates the default value.
                  </Typography>
                </Box>
                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>Resolution: <b>{values.resolution}</b></Typography>
                    <Slider
                    value={values.resolution}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'resolution', value: val } })}
                    valueLabelDisplay="auto"
                    min={0}
                    max={5}
                    step={0.05}
                    name="resolution"
                    marks={[
                      // { value: 0, label: '0' },
                      { value: 0.1, label: '0.1' },
                      // { value: 0.25, label: '0.25' },
                      { value: 0.5, label: '0.5' },
                      { value: 1, label: '1*' },
                      { value: 2.5, label: '2.5' },
                      { value: 5, label: '5' },
                    ]}
                    />
                </Box>
            </FormGroup>
        </ExpansionPanelDetails>
      </ExpansionPanel>

      <div style={{ marginTop: '10px' }}>
        <div>
          <label htmlFor="material-switch">
            <p>Use Default Values</p>
            <Switch
              checked={values.use_default}
              onChange={handleSwitchChange("use_default")}
              onColor="#86d3ff"
              onHandleColor="#2693e6"
              handleDiameter={30}
              // uncheckedIcon={false}
              // checkedIcon={false}
              boxShadow="0px 1px 5px rgba(0, 0, 0, 0.6)"
              activeBoxShadow="0px 0px 1px 10px rgba(0, 0, 0, 0.2)"
              height={20}
              width={48}
              className="react-switch"
              id="use_default"
            />
          </label>
        </div>
      </div>
    </div>
  );
};

export default QualityControlParameters;
