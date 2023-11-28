import RightRail from "../RightNavigation/rightRail";
import LeftNav from "./components/leftNav";
import MiddleContent from "./components/mainContent";
import React, { useState } from 'react';
import './publishDatasets.css'; // Import a CSS file for styles

export default function PublishDataset() {

    const [taskStatus, setTaskStatus] = useState({
        1: false, // Task 1 is initially not completed
        2: false,
        3: false,
        4: false,
        5: false,
        6: false,
        7: false
        // Add other tasks here
      });
    

      const [taskData, setTaskData] = useState({
        upload: {},
        validation: {
            inputFiles: [],
            seuratFiles: [],
            selectedSeuratFile: null,
            fileMappings:[]
          },
        quality_control: {
          qc_results: []
        },
        metadata: {
          formData: {
            Dataset: '',
            Task: '',
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
            'Cell Count Estimate': 0,
            'Source': '',
            'Source Key': '',
            'Submission Date': 'YYYY-MM-DD', // Set your initial date placeholder here    
            'Id':'' 
           },
          taskOptions: [],
          options: {
            Task: [], 
            Author: [],
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
          task_type: '',
          task_id:'',
          task_label: [],
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

    const [activeTask, setActiveTask] = useState(1); // Initialize with the first task

    return(
        <div className="page-container">
            <div className="left-nav">
            <LeftNav activeTask={activeTask} setActiveTask={setActiveTask} taskStatus={taskStatus} taskData={taskData} setTaskData={setTaskData} />
            </div>
            <div className="main-content">
                <MiddleContent activeTask={activeTask} setActiveTask={setActiveTask} setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} taskStatus={taskStatus}/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}