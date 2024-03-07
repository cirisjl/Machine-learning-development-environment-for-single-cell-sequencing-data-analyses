import React, { useState } from 'react';
import Box from '@material-ui/core/Box';
import Typography from '@material-ui/core/Typography';
import Slider from '@material-ui/core/Slider';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Switch from '@material-ui/core/Switch';
import {ExpansionPanel, ExpansionPanelSummary, ExpansionPanelDetails, Button } from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import { makeStyles } from '@material-ui/core/styles';


const QualityControlParameters = ({values, setValues, defaultValues}) => {

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
      fontSize: '14px', // Adjusts font size for panel titles
      '& .MuiExpansionPanelSummary-content': {
        margin: theme.spacing(1), // Optionally adjust spacing if needed
      },
    },
    panelDetails: {
      display: 'block',
      fontSize: '12px', // Adjusts font size for panel content
    },
    formGroup: {
      margin: theme.spacing(2), 
    },
    formControlLabel: { // Optionally if you want to target labels specifically
      '& .MuiTypography-root': {
        fontSize: '12px', // Adjusts font size for form control labels
      },
    },
  }));

  const handleSwitchChange = (event) => {
    const { name, checked } = event.target;
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
     <Typography variant="h4" gutterBottom>
        Advanced Quality Control Parameters
      </Typography>
      
      <ExpansionPanel className={useStyles.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={useStyles.panelSummary}>
            {/* QC Parameters */}
            <Typography variant="h6" gutterBottom>QC Parameters</Typography>
        </ExpansionPanelSummary>
        <ExpansionPanelDetails className={useStyles.panelDetails}>
             <FormGroup className={useStyles.formGroup}>
                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>Min Genes - Max Genes: [{values.min_genes} - {values.max_genes}]</Typography>
                    <Slider
                    value={[values.min_genes, values.max_genes]}
                    onChange={handleGeneSliderChange}
                    valueLabelDisplay="auto"
                    min={0}
                    max={20000}
                    step={25}
                    name="min_max_genes"
                    marks={[{}]}
                    />
                </Box>

                <Box sx={{ m: 2 }}>

                    <Typography gutterBottom>Min Cells: {values.min_cells}</Typography>
                    <Slider
                    value={values.min_cells}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'min_cells', value: val } })}
                    valueLabelDisplay="auto"
                    min={1}
                    max={200}
                    step={1}
                    name="min_cells"
                    />
                </Box>

                <Box sx={{ m: 2 }}>

                    <Typography gutterBottom>Target Sum: {values.target_sum}</Typography>
                    <Slider
                    value={values.target_sum}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'target_sum', value: val } })}
                    //   valueLabelDisplay="auto"
                    min={0}
                    max={1e6}
                    step={1e4}
                    name="target_sum"
                    />
                </Box>

                <Box sx={{ m: 2 }}>

                    <Typography gutterBottom>Highly Variable Genes (n_top_genes): {values.n_top_genes}</Typography>
                    <Slider
                    value={values.n_top_genes}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'n_top_genes', value: val } })}
                    valueLabelDisplay="auto"
                    min={100}
                    max={10000}
                    step={25}
                    name="n_top_genes"
                    />
                </Box>

                <Box sx={{ m: 2 }}>
                  <Typography gutterBottom>
                    Doublet Rate: {`${(values.doublet_rate * 100).toFixed(2)}%`}
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
                      { value: 0, label: '0%' },
                      { value: 0.008, label: '0.8%' },
                      { value: 0.023, label: '2.3%' },
                      { value: 0.038, label: '3.8%' },
                      { value: 0.046, label: '4.6%' },
                      { value: 0.061, label: '6.1%' },
                      { value: 0.08, label: '8%*' }, // Default
                      { value: 0.125, label: '12.5%' },
                      { value: 0.2, label: '20%' },
                      { value: 0.5, label: '50%' },
                    ]}
                    name="doublet_rate"
                  />
                </Box>

                <Box sx={{ m: 2 }}>

                    <FormControlLabel
                        control={<Switch checked={values.regress_cell_cycle} onChange={handleSwitchChange} name="regress_cell_cycle" />}
                        label={`Regress Cell Cycle: ${values.regress_cell_cycle ? 'Yes' : 'No'}`}
                        className={useStyles.formControlLabel}
                    />
                </Box>
            </FormGroup>
            </ExpansionPanelDetails>
      </ExpansionPanel>


      <ExpansionPanel className={useStyles.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={useStyles.panelSummary}>
            {/* Projection Parameters */}
            <Typography variant="h6" gutterBottom>Projection Parameters</Typography>
        </ExpansionPanelSummary>
        
        <ExpansionPanelDetails className={useStyles.panelDetails}>
            <FormGroup className={useStyles.formGroup}>
                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>n_neighbors: {values.n_neighbors}</Typography>
                    <Slider
                    value={values.n_neighbors}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'n_neighbors', value: val } })}
                    valueLabelDisplay="auto"
                    min={2}
                    max={100}
                    step={1}
                    name="n_neighbors"
                    />
                </Box>

                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>n_pcs: {values.n_pcs}</Typography>
                    <Slider
                    value={values.n_pcs}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'n_pcs', value: val } })}
                    valueLabelDisplay="auto"
                    min={0}
                    max={200}
                    step={1}
                    name="n_pcs"
                    />
                </Box>
            </FormGroup>
        </ExpansionPanelDetails>
      </ExpansionPanel> 


      <ExpansionPanel className={useStyles.root}>
        <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} className={useStyles.panelSummary}>
          <Typography variant="h6" gutterBottom>Clustering Parameters</Typography>
        </ExpansionPanelSummary>

        <ExpansionPanelDetails className={useStyles.panelDetails}>
            <FormGroup className={useStyles.formGroup}>

                <Box sx={{ m: 2 }}>
                    <Typography gutterBottom>Resolution: {values.resolution}</Typography>
                    <Slider
                    value={values.resolution}
                    onChange={(e, val) => handleSliderChange({ target: { name: 'resolution', value: val } })}
                    valueLabelDisplay="auto"
                    min={0}
                    max={5}
                    step={0.05}
                    name="resolution"
                    />
                </Box>
            </FormGroup>
        </ExpansionPanelDetails>
      </ExpansionPanel>

        {/* Switches */}
        <FormControlLabel
            control={<Switch checked={values.use_default} onChange={handleSwitchChange} name="use_default" />}
            label="Use Default Values"
            className={useStyles.formControlLabel}
        />
    </div>
  );
};

export default QualityControlParameters;
