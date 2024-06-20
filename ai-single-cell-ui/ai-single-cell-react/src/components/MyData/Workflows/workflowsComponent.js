import React, { useState, useEffect } from 'react';
import { getCookie} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import schema from '../../../schema/react-json-schema/Workflows/clusteringUsingRaceIDSchema.json';
import Form from 'react-jsonschema-form';
import styled from 'styled-components';
import { Container, Button } from "reactstrap";
import Switch from 'react-switch';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import InputDataComponent from '../Tools/inputDataCollection';
import { ClusteringWorkFlowComponent } from './components/clusteringWorkflow';




export function WorkflowsComponent(props) {

  const selectedWorkflow = props.selectedWorkflow;
  console.log(selectedWorkflow);

  return (
    <div>
      {selectedWorkflow.toLowerCase() === "clustering" && <ClusteringWorkFlowComponent selectedWorkflow={selectedWorkflow}/>}
    </div>
  )
};
