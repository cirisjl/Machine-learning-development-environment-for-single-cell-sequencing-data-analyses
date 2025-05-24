import React , { useState, useEffect ,useRef}from 'react';
import Select from 'react-select';
import { Card, CardContent, Typography, Slider, Button } from '@material-ui/core';
import { CELERY_BACKEND_API, NODE_API_URL, WEB_SOCKET_URL} from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';
import axios from 'axios';
import ReactPlotly from './reactPlotly';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import {faFile} from "@fortawesome/free-solid-svg-icons";
import DatasetSelectionDialog from './datasetsDialog';
import TableComponent from './labelTableComponent';
import useWebSockets from './useWebSockets';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [ isError, setIsError ] = useState(false);
  const [loading, setLoading] = useState({});
  const [jobId, setjobId] = useState("");
  const webSocketStatus = useRef(null);

  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const [selectionMode, setSelectionMode] = useState(''); // or 'multiple'

  const taskOptions = [
    { label: "Clustering", value: "CL" },
    { label: "Imputation", value: "IM" },
    { label: "Marker gene identification", value: "MGI" },
    { label: "Trajectory", value: "TR" },
    { label: "Cell-cell communication", value: "CCC" },
    { label: "Multi-omic data integration", value: "MDI" },
    { label: "Gene regulatory relations", value: "GRR" },
    { label: "Cell type identification", value: "CTI" },
    { label: "Spatial", value: "SP" }
  ];

  const handleOpenDialog = (mode) => {
    if (selectionMode !== mode) {
      // Reset selectedDatasets in taskData when the selection mode changes
      setTaskData(prevTaskData => ({
        ...prevTaskData,
        task_builder: {
          ...prevTaskData.task_builder,
          selectedDatasets: {}, // Resetting selectedDatasets when the mode changes
        },
      }));
    }
    setSelectionMode(mode);
    setIsDialogOpen(true);
  };

  const updateTaskDataWithResults = (data) => {
    // Update taskData state with results
    // For example, setting archive_path
    setTaskData(prevTaskData => ({
        ...prevTaskData,
        task_builder: {
            ...prevTaskData.task_builder,
            selectedDatasets: {
                ...prevTaskData.task_builder.selectedDatasets,
                [data.datasetId]: {
                    ...prevTaskData.task_builder.selectedDatasets[data.datasetId],
                    dataSplit: {
                        ...prevTaskData.task_builder.selectedDatasets[data.datasetId].dataSplit,
                        adataPath: data.adata_path,
                        dataSplitPerformed: true,
                    }
                }
            }
        }
    }));
  };
  const handleStatusMessage = (event) => {
    try {
      const data = JSON.parse(event.data);
      console.log("Task Response");
      console.log(data);
      if (data.task_status && (data.task_status === "SUCCESS" || data.task_status === "FAILURE")) {
        if(data.task_status === "SUCCESS"){
            updateTaskDataWithResults(data.task_result); // Update your state with results
            setLoading(prevLoading => ({ ...prevLoading, [data.task_result.datasetId]: false })); 
        }
        // Close WebSocket connection for this jobId
        closeWebSocket(jobId);
    }
    } catch (error) {
      console.error("Error parsing status message:", error);
      // setLoading(prevLoading => ({ ...prevLoading, [datasetId]: false })); 
    }
  };

  const { closeWebSocket } = useWebSockets(jobId, handleStatusMessage, WEB_SOCKET_URL);

  const handleSelectDatasets = async (newSelectedDatasets) => {
    // Initialize additional parameters for new datasets
    Object.keys(newSelectedDatasets).forEach(key => {
      if (!taskData.task_builder.selectedDatasets[key]) {
        newSelectedDatasets[key] = {
          ...newSelectedDatasets[key],
          taskType: null,
          taskLabel: '',
          dataSplit: {
            trainFraction: 1,
            validationFraction: 0,
            testFraction: 0,
            dataSplitPerformed: false,
            archivePath: ''
          }
        };
      }
    });

    setTaskData(prevTaskData => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        selectedDatasets: newSelectedDatasets
      },
    }));
  };

  const handleCloseDialog = () => {
    setIsDialogOpen(false);
  };

  const handleDataSplit = async (adata_path, datasetId) => {
    try {
      const dataset = taskData.task_builder.selectedDatasets[datasetId];
      if(!dataset.taskType) {
        setMessage('You must select taskType before the data split!');
        setHasMessage(true);
        setIsError(true);
        return;
      }
      setLoading(prevLoading => ({ ...prevLoading, [datasetId]: true })); // Set loading to true for the specific dataset
  
      const userData = {
        benchmarksId: dataset.taskType.value + "-" + datasetId,
        datasetId: datasetId,
        userId: dataset.Owner,
        adata_path: adata_path || '',
        train_fraction: dataset.dataSplit.trainFraction,
        validation_fraction: dataset.dataSplit.validationFraction,
        test_fraction: dataset.dataSplit.testFraction,
      };
  
      const totalFraction = userData.train_fraction + userData.validation_fraction + userData.test_fraction;
  
      if (totalFraction !== 1) {
        setHasMessage(true);
        setMessage("The sum of train, validation, and test fractions must equal 1.");
        setLoading(false);
        setIsError(true);
        return;
      }
  
      console.log(userData);
  
      // Make the API call
      const response = await fetch(`${CELERY_BACKEND_API}/benchmarks/data-split`, {        
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(userData),
      });
  
      if (response.ok) {
        const result = await response.json();
        const jobId = result.job_id;
        setjobId(jobId);
  
        // Update only the specific dataset's dataSplit parameters
        setTaskData(prevTaskData => ({
          ...prevTaskData,
          task_builder: {
            ...prevTaskData.task_builder,
            selectedDatasets: {
              ...prevTaskData.task_builder.selectedDatasets,
              [datasetId]: {
                ...prevTaskData.task_builder.selectedDatasets[datasetId],
                dataSplit: {
                  ...prevTaskData.task_builder.selectedDatasets[datasetId].dataSplit,
                  dataSplitPerformed: true,
                  adata_path: result.adata_path,
                }
              }
            }
          }
        }));
      } else {
        const error = await response.json();
        console.error(error.error); // Handle the error
        setLoading(prevLoading => ({ ...prevLoading, [datasetId]: false })); // Set loading to false for the specific dataset
      }
    } catch (error) {
      console.error('Error:', error);
      setLoading(prevLoading => ({ ...prevLoading, [datasetId]: false })); // Set loading to false for the specific dataset
    } 
  };
  

  const handleTaskCompletion = () => {
    // Check if all datasets have a task type, label and data split performed
    const allDatasetsValid = Object.values(taskData.task_builder.selectedDatasets).every(dataset => 
      dataset.taskType && 
      dataset.taskLabel && 
      dataset.dataSplit.dataSplitPerformed
    );
  
    if (allDatasetsValid) {
      setTaskStatus(prevTaskStatus => ({
        ...prevTaskStatus,
        5: true, // Mark Task 5 as completed
      }));
  
      setActiveTask(5);
    } else {
      setMessage('Please ensure that the task type, labels, and data split for each dataset are valid.');
      setHasMessage(true);
      setIsError(true);
    }
  };
  

  const handleTaskTypeChange = (datasetId, newTaskType) => {
    setTaskData(prevTaskData => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        selectedDatasets: {
          ...prevTaskData.task_builder.selectedDatasets,
          [datasetId]: {
            ...prevTaskData.task_builder.selectedDatasets[datasetId],
            taskType: newTaskType
          }
        }
      }
    }));
  };


  const handleDataSplitChange = (datasetId, parameter, value) => {
    setTaskData(prevTaskData => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        selectedDatasets: {
          ...prevTaskData.task_builder.selectedDatasets,
          [datasetId]: {
            ...prevTaskData.task_builder.selectedDatasets[datasetId],
            dataSplit: {
              ...prevTaskData.task_builder.selectedDatasets[datasetId].dataSplit,
              [parameter]: value
            }
          }
        }
      }
    }));
  };
  
  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

  const handleLabelChange = (datasetId, selectedOption) => {
    setTaskData(prevTaskData => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        selectedDatasets: {
          ...prevTaskData.task_builder.selectedDatasets,
          [datasetId]: {
            ...prevTaskData.task_builder.selectedDatasets[datasetId],
            taskLabel: selectedOption
          }
        }
      }
    }));
  };

  const onSelectDataset = (dataset) => {
    let datasetId = dataset.Id; 
    let currentSelectedDatasets = { ...taskData.task_builder.selectedDatasets };

    if(currentSelectedDatasets[datasetId]) {
      delete currentSelectedDatasets[datasetId];
    } else {
      if(selectionMode === "single") {
        currentSelectedDatasets = {};
      }
      currentSelectedDatasets[datasetId] = dataset;
    }

    // Call onSelect with the updated selected datasets
    handleSelectDatasets(currentSelectedDatasets);
  };

    // Function to handle selection of sub-items
const onSelectSubItem = (mainItem, subItem) => {
  const mainItemId = mainItem.Id;
  let currentSelectedDatasets = { ...taskData.task_builder.selectedDatasets };

  // Check if the main item is already selected
  if (currentSelectedDatasets[mainItemId]) {
      // If sub-item is already selected, deselect it
      if (currentSelectedDatasets[mainItemId].selectedSubItem?.process_id  === subItem.process_id ) {
          delete currentSelectedDatasets[mainItemId];
      } else {
          // Update the selected main item with the selected sub-item
          currentSelectedDatasets[mainItemId] = {
              ...mainItem,
              selectedSubItem: subItem
          };
      }
  } else {
      // Select the main item and the sub-item
      currentSelectedDatasets = {
          [mainItemId]: {
              ...mainItem,
              selectedSubItem: subItem
          }
      };
  }

    Object.keys(currentSelectedDatasets).forEach(key => {
      currentSelectedDatasets[key] = {
          ...currentSelectedDatasets[key],
          taskType: null,
          taskLabel: '',
          dataSplit: {
            trainFraction: 0.8,
            validationFraction: 0.1,
            testFraction: 0.1,
            dataSplitPerformed: false,
            archivePath: ''
          }
        };
    });
    setTaskData(prevTaskData => ({
      ...prevTaskData,
      task_builder: {
          ...prevTaskData.task_builder,
          selectedDatasets: currentSelectedDatasets
      },
  }));

  };
  

  return (
    <div className='task-builder-task'>
      <div>
        <div>
          <div>
            <button onClick={() => handleOpenDialog('single')}>Select Single Dataset</button>
            <button onClick={() => handleOpenDialog('multiple')}>Select Multiple Datasets</button>
          </div>
          {isDialogOpen && (
            <DatasetSelectionDialog
              onSelect={onSelectDataset}
              multiple={selectionMode === 'multiple'}
              onClose={handleCloseDialog}
              isVisible={isDialogOpen !== false}
              selectedDatasets={taskData.task_builder.selectedDatasets}
              onSelectSubItem = {onSelectSubItem}
            />
          )}
        </div>
      </div>
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError} />}
      <div>

    {Object.entries(taskData.task_builder.selectedDatasets).length > 0 && (
      <div className="metadata-section">
        <Typography variant="h6" component="h6">Select task type and choose the label for each dataset accordingly.</Typography>
          {Object.entries(taskData.task_builder.selectedDatasets).map(([key, dataset], index) => (
            <Card key={index} className="metadata-item">
              <CardContent>
              <Typography variant="h6" component="h6">Dataset {index + 1} : {dataset.Id}</Typography>
              <Typography variant="body2" component="p">
                <label>
                  <p>Please Choose the Task Type:</p>
                    <Select
                      value={dataset.taskType}
                      options={taskOptions}
                      onChange={(selectedOption) => handleTaskTypeChange(key, selectedOption)}
                      />
                </label>
              </Typography>

              <Typography variant="body2" component="p">
                  Please Choose the Label:
                </Typography>
                <Select
                  options={Object.keys(JSON.parse(dataset.cell_metadata_head)).map((key) => ({
                    label: key,
                    value: key,
                  }))}
                  onChange={(selectedOption) => handleLabelChange(key, selectedOption)}
                  value={dataset.taskLabel}
                />
              {/* 
                  {dataset.tablePlot && (
                    <>
                      <h2>Table: </h2>
                      <ReactPlotly plot_data={dataset.tablePlot} />
                    </>
                  )} */}

                <div className="label-table-container">
                  <TableComponent cellMetadataObs={JSON.parse(dataset.cell_metadata_head)} />
                </div>
                
                <Typography variant="body2" component="p" style={{ marginTop: '20px' }}>
                  Data Split Parameters
                </Typography>
               
                <Typography variant="body2" component="p">

               {/* Slider input for Train Fraction */}
               <label>
                      <p>Train Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={dataset.dataSplit.trainFraction}
                        onChange={(e) => handleDataSplitChange(key, 'trainFraction', parseFloat(e.target.value))}

                      />
                        {dataset.dataSplit.trainFraction}
                    </label>

                    {/* Slider input for Validation Fraction */}
                    <label>
                     <p> Validation Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={dataset.dataSplit.validationFraction}
                        onChange={(e) => handleDataSplitChange(key, 'validationFraction', parseFloat(e.target.value))}

                      />
                       {dataset.dataSplit.validationFraction}
                    </label>

                    {/* Slider input for Test Fraction */}
                    <label>
                      <p>Test Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={dataset.dataSplit.testFraction}
                        onChange={(e) => handleDataSplitChange(key, 'testFraction', parseFloat(e.target.value))}

                      />
                      {dataset.dataSplit.testFraction}
                    </label>

                
                <Button 
                  onClick={() => handleDataSplit(dataset.adata_path, key)} 
                  disabled={dataset.dataSplit.dataSplitPerformed || loading[key]} // Use dataset-specific loading state
                  variant="contained" 
                  color="primary"
                  style={{ marginTop: '20px' }}
                >
                {loading[key] ? 'Processing, please wait...' : 'Perform Data Split'} 
                </Button>

                </Typography>

                {dataset.dataSplit.dataSplitPerformed && (
                  <Typography variant="body2" component="p">
                    <b>Archive Path: </b>{dataset.dataSplit.archivePath}
                  </Typography>
                )}
              </CardContent>
            </Card>
        ))}
      </div>
    )}

      <div className='navigation-buttons'>
            {/* <div className="previous">
              <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                Previous
              </button>
            </div> */}
            <div className="next-upon-success">
              <button type="submit" className="btn btn-info button" onClick={handleTaskCompletion}>
                Next
              </button>
            </div>
      </div>

    </div>
    </div>
  );
}

export default TaskBuilderTaskComponent;
