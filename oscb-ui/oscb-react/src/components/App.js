import './App.css';
import {
  BrowserRouter,
  Route, 
  Routes
} from 'react-router-dom'
import React, { useState } from 'react';

// import Layouts
import RootLayout from './../layouts/rootLayout'

// import Pages
import GetStarted from './../pages/getStarted'
import Competitions from './../pages/competitions'
import Updates from './../pages/updates'
import Benchmarks from './../pages/benchmarks'
import Leaderboards from './../pages/leaderboards'
import MyData from '../pages/MyData/mydata'
import Team from './../pages/team'
import Tutorial from './../pages/tutorial'
import PreviewDatasets from '../pages/MyData/previewDatasets'
import UploadData from './MyData/uploadData';
import Login from '../pages/login/login';
import SignUp from '../pages/login/signup';
import RoutingTemplate from '../pages/login/loginRouting';
import WorkflowsPageComponent from '../pages/MyData/Workflows/workflowsPageComponent';
import ToolsComponentPage from '../pages/MyData/Tools/toolsComponentPage';
import MyTasks from '../pages/myTasks';
import ResultFiles from '../pages/ResultFiles';
import FlaskDashboard from './MyData/dashboard';
import NewApp from './Form/Components/component2';
import MyForm from './Form/Components/customComponent';
import AccessDenied from './AccessDeniedPage';
import ManageOptions from './Form/Components/editablePageOptions';
// import PublishDataset from './publishDatasets/publishDataset';
import FlowControl from './publishDatasets/flowControl';
import TaskResultsComponent from './Benchmarks/taskResultsComponent';
import UploadDataset from './MyData/UploadData/uploadDataset';
import SessionReminder from './Session/sessionManager';
import { SessionProvider } from './Session/context/sessionContext'; 
import QualityControlParameters from './publishDatasets/components/qualityControlParameters';
import TaskDetailsComponent from './MyData/MyTasks/taskDetailsComponent';
import TreeTableComponent from './common_components/treeTableComponent';
import BenchmarksViewDetailsComponent from './Benchmarks/components/benchmarksViewDetailsComponent';

function App() {
  
  return (
    <>
    <div className="session-manager-component"><SessionReminder/> </div>
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<RootLayout />}>
          <Route path='getStarted'   element={<GetStarted/>} />
          <Route path="updates"      element={<Updates/>} />
          <Route path="competitions" element={<Competitions/>}/>
          <Route path="benchmarks"   element={<Benchmarks/>}/>
          <Route path="benchmarks/uploads"   element={<FlowControl/>}/>
          <Route path="benchmarks/clustering"   element={<TaskResultsComponent task_type="Clustering"/>}/>
          <Route path="benchmarks/viewDetails"   element={<BenchmarksViewDetailsComponent/>}/>
          <Route path="benchmarks/imputation" element={<TaskResultsComponent task_type="Imputation" />} />
          <Route path="benchmarks/maker-gene-identification" element={<TaskResultsComponent task_type="Marker Gene Identification" />} />
          <Route path="benchmarks/trajectory" element={<TaskResultsComponent task_type="Trajectory" />} />
          <Route path="benchmarks/cell-cell-communication" element={<TaskResultsComponent task_type="Cell-Cell Communication" />} />
          <Route path="benchmarks/multiomics-data-integration" element={<TaskResultsComponent task_type="Multiomics Data Integration" />} />
          <Route path="benchmarks/gene-regulatory-relations" element={<TaskResultsComponent task_type="Gene Regulatory Relations" />} />
          <Route path="benchmarks/genes-over-time" element={<TaskResultsComponent task_type="Genes Over Time" />} />
          <Route path="benchmarks/genes-over-condition" element={<TaskResultsComponent task_type="Genes Over Condition" />} />
          <Route path="benchmarks/cell-type" element={<TaskResultsComponent task_type="Cell Type" />} />
          <Route path="leaderboards" element={<Leaderboards/>}/>
          <Route path="mydata"       element={<MyData/>}></Route>
          <Route path="mydata/upload-data"       element={<UploadDataset/>}></Route>
          <Route path="mydata/update-dataset"       element={<UploadData/>}></Route>
          <Route path="mydata/preview-datasets" element={<PreviewDatasets/>}></Route>
          <Route path="mydata/taskDetails"       element={<TaskDetailsComponent/>}></Route>
          <Route path="mydata/workflows" element={<WorkflowsPageComponent/>}></Route>
          <Route path="mydata/tools" element={<ToolsComponentPage/>}></Route>
          <Route path="team"         element={<Team/>}/>
          <Route path="dashboard"         element={<FlaskDashboard/>}/>
          <Route path="tutorial" element={<Tutorial/>}/>
          <Route path="login"         element={<Login/>}/>
          <Route path="signup"         element={<SignUp/>}/>
          <Route path="routing"         element={<RoutingTemplate/>}/>
          <Route path="myTasks"         element={<MyTasks/>}/>
          <Route path="resultfiles"         element={<ResultFiles/>}/>
          <Route path="new"         element={<NewApp/>}/>
          <Route path="custom"         element={<MyForm/>}/>
          <Route path="accessDenied"         element={<AccessDenied/>}/>
          <Route path="manageOptions"         element={<ManageOptions/>}/>
          <Route path="testapi"         element={<TreeTableComponent/>}/>
          {/* <Route path="publishDataset"         element={<PublishDataset/>}/> */}

        </Route>
      </Routes>
    </BrowserRouter>
    </>
  );
}

export default App
