import React, { useState, useEffect } from 'react';
import InputDataComponent from '../../Tools/inputDataCollection';
import Schema from '../../../../schema/react-json-schema/Workflows/clusteringWorkflow.json'
import { uiSchema } from '../../../../schema/UI-schema/Workflows/clusteringWorkflow';
import Form from '@rjsf/core';
import validator from '@rjsf/validator-ajv8';
import SelectComponent from '../../Tools/components/selectComponent';
import GeneRangeSlider from '../../Tools/components/geneRangeSlider';
import MultiSelectComponent from '../../Tools/components/multiselectComponent';
import RangeSlider from '../../Tools/components/sliderComponent';
import SwitchComponent from '../../Tools/components/switchComponent';
import UseDefaultSwitch from '../../Tools/components/useDefaultSwitch';
import ClusterLabelInput from '../../Tools/components/customInputComponent';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';

export function ClusteringWorkFlowComponent(props) {
    const [selectedDatasets, setSelectedDatasets] = useState({});
    const [formErrors, setFormErrors] = useState("");
    const [value, setValue] = useState("");
    const [formData, setFormData] = useState({});


    const handleChange = (event) => {
        setValue(event.target.value);
      };

    const onDeleteDataset = (id) => {
    const currentSelectedDatasets = { ...selectedDatasets};
    
    if (currentSelectedDatasets[id]) {
        delete currentSelectedDatasets[id];
    }
    setSelectedDatasets(currentSelectedDatasets);
    };

    const onSelectDataset = (dataset) => {
        let datasetId = dataset.Id; 
        let currentSelectedDatasets = { ...selectedDatasets};
    
        if(currentSelectedDatasets[datasetId]) {
          delete currentSelectedDatasets[datasetId];
        } else {
          if(props.selectedWorkflow !== "integration") {
            currentSelectedDatasets = {};
          }
          currentSelectedDatasets[datasetId] = dataset;
        }
      setSelectedDatasets(currentSelectedDatasets)
    };

    const widgets = {
        SelectComponent: SelectComponent,
        geneRangeSlider: GeneRangeSlider,
        MultiSelectComponent: MultiSelectComponent,
        toggle: (props) => (
          <Toggle
            checked={props.value}
            onChange={(e) => props.onChange(e.target.checked)}
          />
        ),
        GeneRangeSlider: GeneRangeSlider,
        RangeSlider: RangeSlider,
        SwitchComponent: SwitchComponent,
        UseDefaultSwitch: UseDefaultSwitch,
        ClusterLabelInput: ClusterLabelInput
      };

    const handleSubmit = ({ formData }) => {
        console.log("Submitted data:", formData);
    };
    
  return (
    <div className='workflow-container common-class-tools-and-workflows'>
      <div className="separator heading">
        <div className="stripe"></div> 
          <h2 className="h-sm font-weight-bold">
            Tool Parameters 
          </h2> 
        <div className="stripe"></div>
      </div>

      <div>
        <InputDataComponent formErrors={formErrors} filterCategory={props.selectedWorkflow} selectedDatasets={selectedDatasets}
            onSelectDataset={onSelectDataset} onDeleteDataset={onDeleteDataset}/>
      </div>

      <div>
        <Form
            schema = {Schema}
            uiSchema={uiSchema}
            widgets={widgets}
            formData={formData}
            onChange={({ formData }) => setFormData(formData)}
            onSubmit={handleSubmit}
            onError={(errors) => console.log("Form errors:", errors)}
            validator={validator}
        />
      </div>

    </div>
  )
};
