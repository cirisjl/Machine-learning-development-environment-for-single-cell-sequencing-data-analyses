import React , { useState, useEffect }from 'react';
import Select from 'react-select';
import { Card, CardContent, Typography, Slider, Button } from '@material-ui/core';
import { CELERY_BACKEND_API, SERVER_URL} from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';
import axios from 'axios';
import ReactPlotly from './reactPlotly';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import {faFile} from "@fortawesome/free-solid-svg-icons";
import DatasetSelectionDialog from './datasetsDialog';
import TableComponent from './labelTableComponent';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [loading, setLoading] = useState({});

  const [datasets, setDatasets] = useState([]);

  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const [selectionMode, setSelectionMode] = useState(''); // or 'multiple'

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
  const handleSelectDatasets = async (newSelectedDatasets) => {
    // Initialize additional parameters for new datasets
    Object.keys(newSelectedDatasets).forEach(key => {
      if (!taskData.task_builder.selectedDatasets[key]) {
        newSelectedDatasets[key] = {
          ...newSelectedDatasets[key],
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
      }
    });

  //   // Extract file paths for new datasets to make an API call
  // const filePaths = Object.values(newSelectedDatasets).map(dataset => dataset.adata_path).filter(Boolean);

  // if (filePaths.length > 0) {
  //   try {
  //     // Make an API call to get the table plots for new datasets
  //     const response = await fetch(`${CELERY_BACKEND_API}/convert/api/getTablePlot`, {
  //       method: 'POST',
  //       headers: {
  //         'Content-Type': 'application/json',
  //       },
  //       body: JSON.stringify(filePaths),
  //     });

  //     if (!response.ok) {
  //       throw new Error('Network response was not ok');
  //     }

  //     const tablePlots = await response.json();

  //     // Assuming tablePlots is an array of results corresponding to filePaths,
  //     // Attach each tablePlot result to its respective dataset
  //     filePaths.forEach((filePath, index) => {
  //       const key = Object.keys(newSelectedDatasets).find(key => newSelectedDatasets[key].adata_path === filePath);
  //       if (key) {
  //         newSelectedDatasets[key].tablePlot = tablePlots[index];
  //       }
  //     });

  //   } catch (error) {
  //     console.error('There was a problem with the fetch operation:', error);
  //   }
  // }
  
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
      setLoading(prevLoading => ({ ...prevLoading, [datasetId]: true })); // Set loading to true for the specific dataset
      const dataset = taskData.task_builder.selectedDatasets[datasetId];
  
      const userData = {
        data: adata_path || '',
        train_fraction: dataset.dataSplit.trainFraction,
        validation_fraction: dataset.dataSplit.validationFraction,
        test_fraction: dataset.dataSplit.testFraction,
      };
  
      const totalFraction = userData.train_fraction + userData.validation_fraction + userData.test_fraction;
  
      if (totalFraction !== 1) {
        setHasMessage(true);
        setMessage("The sum of train, validation, and test fractions must equal 1.");
        setLoading(false);
        return;
      }
  
      console.log(userData);
  
      // Make the API call
      const response = await fetch(`${CELERY_BACKEND_API}/convert/api/data-split`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(userData),
      });
  
      if (response.ok) {
        const result = await response.json();
        console.log(result); // Handle the result as needed
  
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
                  archivePath: result.archive_path,
                }
              }
            }
          }
        }));
      } else {
        const error = await response.json();
        console.error(error.error); // Handle the error
      }
    } catch (error) {
      console.error('Error:', error);
    } finally {
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
  
      setActiveTask(6);
    } else {
      setMessage('Please ensure that the task type, labels, and data split for each dataset are valid.');
      setHasMessage(true);
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
              onSelect={handleSelectDatasets}
              multiple={selectionMode === 'multiple'}
              onClose={handleCloseDialog}
              isVisible={isDialogOpen !== false}
              taskData={taskData}
            />
          )}
        </div>
      </div>
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} />}
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
                      options={dataset.taskOptions}
                      onChange={(selectedOption) => handleTaskTypeChange(key, selectedOption)}
                      />
                </label>
              </Typography>

              <Typography variant="body2" component="p">
                  Please Choose the Label:
                </Typography>
                <Select
                  options={Object.keys(JSON.parse(dataset.cell_metadata_obs)).map((key) => ({
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
                  <TableComponent cellMetadataObs={JSON.parse(dataset.cell_metadata_obs)} />
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
