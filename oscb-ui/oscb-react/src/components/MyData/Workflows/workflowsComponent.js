import React from 'react';
import { ClusteringWorkFlowComponent } from './components/clusteringWorkflow';
import { IntegrationWorkFlowComponent } from './components/integrationWorkflow';
import { AnnotationWorkFlowComponent } from './components/annotationWorkflow';
import { PhoenixWorkFlowComponent } from './components/phoenixWorkflow';
import { MycoSNPWorkFlowComponent } from './components/mycosnpWorkflow'; 

export function WorkflowsComponent(props) {

  const selectedWorkflow = props.selectedWorkflow;
  console.log(selectedWorkflow);

  return (
    <div>
      {selectedWorkflow.toLowerCase() === "clustering" && <ClusteringWorkFlowComponent selectedWorkflow={selectedWorkflow}/>}
      {selectedWorkflow.toLowerCase() === "integration" && <IntegrationWorkFlowComponent selectedWorkflow={selectedWorkflow} />}
      {selectedWorkflow.toLowerCase() === "annotation" && <AnnotationWorkFlowComponent selectedWorkflow={selectedWorkflow} />}
      {selectedWorkflow.toLowerCase() === "phoenix" && <PhoenixWorkFlowComponent selectedWorkflow={selectedWorkflow} />}
      {selectedWorkflow.toLowerCase() === "mycosnp" && <MycoSNPWorkFlowComponent selectedWorkflow={selectedWorkflow} />}
    </div>
  )
};
