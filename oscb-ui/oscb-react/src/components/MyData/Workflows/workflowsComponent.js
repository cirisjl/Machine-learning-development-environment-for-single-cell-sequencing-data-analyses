import React from 'react';
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
