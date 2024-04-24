import React, { useState } from 'react';
import PublishDataset from './publishDataset';
import TaskBuilder from './taskBuilder';

const FlowControl = () => {

  const [flow, setFlow] = useState('');

  const [activeTask, setActiveTask] = useState(4); // Initialize with the first task

  const [taskStatus, setTaskStatus] = useState({
    1: false, // Task 1 is initially not completed
    2: false,
    3: false,
    4: false,
    5: false,
    6: false,
    // Add other tasks here
  });


  const [taskData, setTaskData] = useState({
    upload: {
      files:[]
    },
    // validation: {
    //     inputFiles: [],
    //     seuratFile: {
    //       default_assay: '',
    //       assay_names: [],
    //       file: ''
    //     } ,
    //     selectedAssayName: '',
    //     fileMappings:[],
    //     token: '',
    //     displayAssayNames: false
    //   },
    quality_control: {
      qc_results: [],
      file_paths: [],
      output: '',
      seurat_meta: {
        default_assay: '',
        assay_names: [],
        file: '',
        displayAssayNames: false
      } ,
      shouldHideForSeurat: false,
      token: '',
      selectedAssayName:''
    },
    metadata: {
      formData: {
        Dataset: '',
        Downloads: '',
        Title: '',
        Author: '',
        'Reference (paper)':'',
        Abstract: '',
        DOI: '',
        Species: '',
        'Sample Type': '',
        'Anatomical Entity': '',
        'Organ Part': '',
        'Model Organ': '',
        'Selected Cell Types': '',
        'Library Construction Method': '',
        'Nucleic Acid Source': '',
        'Paired End': false,
        'Analysis Protocol': '',
        'Disease Status (Specimen)': '',
        'Disease Status (Donor)': '',
        'Development Stage': '',
        'Donor Count': 0,
        'Source': '',
        'Source Key': '',
        'Submission Date': '', // Set your initial date placeholder here    
       },
      taskOptions: [],
      options: {
        Task: [], 
        Author: '',
        Species: [],
        'Sample Type':[],
        'Anatomical Entity': [],
        'Organ Part': [],
        'Model Organ': [],
        'Selected Cell Types': [],
        'Library Construction Method': [],
        'Nucleic Acid Source': [],
        'Disease Status (Specimen)': [],
        'Disease Status (Donor)': [],
        'Development Stage': [],
        'Cell Count Estimate': [],
        'Source': []
      },
      newOptions: []
    },
    task_builder: {
      selectedDatasets: {},
      task_type: '',
      task_id:'',
      task_label: [],
      table_data: [],
      task_data_split:[],
      task_states: {
        trainFraction: 0.8,
        validationFraction: 0.1,
        testFraction: 0.1,
        dataSplitPerformed: false,
        archivePath: ''
      }
    },
    benchmarks: {
      benchmarks_results: []
    },
    Review: {}
  });


  const startFromBeginning = () => {
    setFlow('upload');
    setActiveTask(1);
  };

  const startFromFourthStep = () => {
    setFlow('taskBuilder');
    setActiveTask(4);
  };

  return (
    <div>
        {flow === '' && 
            <div className='flow-messaging'>
                <h3>If you want to create a new dataset and then build tasks, click on "Create a New Dataset". If you want to use existing datasets to build tasks, click on "Jump to Task Builder".</h3>
                <button onClick={startFromBeginning}>Create a new Dataset</button>
                <button onClick={startFromFourthStep}>Jump to Task builder</button>
            </div>
        }
      
        <div>
            {flow === 'upload' && <div><PublishDataset taskStatus={taskStatus} setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} activeTask={activeTask} setActiveTask={setActiveTask} flow={flow} setFlow={setFlow}/></div>}
            {flow === 'taskBuilder' && <div><TaskBuilder taskStatus={taskStatus} setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} activeTask={activeTask} setActiveTask={setActiveTask} flow={flow} setFlow={setFlow}/></div>}
        </div>
    </div>
  );
};

export default FlowControl;
