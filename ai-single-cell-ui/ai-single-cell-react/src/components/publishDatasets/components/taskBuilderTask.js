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

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [loading, setLoading] = useState(false);

  const [datasets, setDatasets] = useState([]);

  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const [selectionMode, setSelectionMode] = useState(''); // or 'multiple'
  const [selectedDatasets, setSelectedDatasets] = useState({});


  const handleOpenDialog = (mode) => {
    if(selectionMode !== mode) {
      setSelectedDatasets({});
    }
    setSelectionMode(mode);
    setIsDialogOpen(true);
  };

  const handleSelectDatasets = (datasets) => {
    console.log(datasets); // Replace this with what you want to do with the selected datasets
  };

  const handleCloseDialog = () => {
    setIsDialogOpen(false);
  };

  const handleDataSplit = async (index) => {
    try {

      setLoading(true); // Set loading to true when data split is initiated

      const dataPath = taskData.quality_control.qc_results[index].adata_path || ''
      // User input for data split fractions
      const userData = {
        data: dataPath,
        train_fraction: taskData.task_builder.task_states.trainFraction,
        validation_fraction: taskData.task_builder.task_states.validationFraction,
        test_fraction: taskData.task_builder.task_states.testFraction,
      };

      console.log(userData);

      // Make the API call
      const response = await fetch(`${CELERY_BACKEND_API}/convert/api/data-split`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(userData),
      });

      // Check if the response is successful
      if (response.ok) {
        const result = await response.json();
        console.log(result); // Handle the result as needed
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          task_builder: {
            ...prevTaskData.task_builder,
            task_states: {
              ...prevTaskData.task_builder.task_states,
              dataSplitPerformed: true,
              archivePath: result.archive_path,
            },
          },
        }));
      } else {
        const error = await response.json();
        console.error(error.error); // Handle the error
      }
    } catch (error) {
      console.error('Error:', error);
    } finally {
      setLoading(false); // Set loading to false after data split (success or failure)
    }
  };

  const handleTaskCompletion = () => {
    const isValid =
            taskData.task_builder.task_type &&
            taskData.task_builder.task_label.length > 0 &&
            taskData.task_builder.task_states.dataSplitPerformed
            
    if(isValid) {
        // After Task 5 is successfully completed, update the task status
        setTaskStatus((prevTaskStatus) => ({
          ...prevTaskStatus,
          5: true, // Mark Task 5 as completed
        }));
  
        //The current task is finished, so make the next task active
        setActiveTask(6);
    } else {
      // Display an error message or throw an error
      setMessage('Please ensure that the task type, labels, dataset split, and data are valid.');
      setHasMessage(true);
    }
  };

  const handleTaskChange = (selectedOption) => {
    // Update the task_type in task_builder
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        task_type: selectedOption,
      },
    }));
  };

  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

  const handleLabelChange = (selectedOption, index) => {

    // Update the task_label in task_builder for the specific result
    setTaskData((prevTaskData) => {
      const updatedLabels = [...prevTaskData.task_builder.task_label];
      updatedLabels[index] = selectedOption;

      return {
        ...prevTaskData,
        task_builder: {
          ...prevTaskData.task_builder,
          task_label: updatedLabels,
        },
      };
    });
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
              selectedDatasets={selectedDatasets}
              setSelectedDatasets={setSelectedDatasets}
            />
          )}
        </div>
      </div>
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} />}
      <div>
        <div className="task-section">
          <label>
            <p>Please Choose the Task Type:</p>
              <Select
                value={taskData.task_builder.task_type}
                options={taskData.metadata.taskOptions}
                onChange={handleTaskChange}
            />
          </label>
          {taskData.task_builder.task_id && (
            <div className="task-id-section">
              Task ID: {taskData.task_builder.task_id}
            </div>
          )}
        </div>

    {Object.entries(selectedDatasets).length > 0 && (
      <div className="metadata-section">
          {Object.entries(selectedDatasets).map(([key, dataset], index) => (
            <Card key={index} className="metadata-item">
              <CardContent>
                <Typography variant="h6" component="h2">
                  Please Choose the Label:
                </Typography>
                <Select
                  options={Object.keys(JSON.parse(dataset.cell_metadata_obs)).map((key) => ({
                    label: key,
                    value: key,
                  }))}
                  onChange={(selectedOption) => handleLabelChange(selectedOption, index)}
                  value={taskData.task_builder.task_label[index]}
                />
                
                <Typography variant="h6" component="h2" style={{ marginTop: '20px' }}>
                  Data Split Parameters
                </Typography>
               
               {/* Slider input for Train Fraction */}
               <label>
                      <p>Train Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={taskData.task_builder.task_states.trainFraction}
                        onChange={(e) => {
                          setTaskData((prevTaskData) => ({
                            ...prevTaskData,
                            task_builder: {
                              ...prevTaskData.task_builder,
                              task_states: {
                                ...prevTaskData.task_builder.task_states,
                                dataSplitPerformed: false,
                                trainFraction: parseFloat(e.target.value),
                              },
                            },
                          }));
                        }}
                      />
                      {taskData.task_builder.task_states.trainFraction}
                    </label>

                    {/* Slider input for Validation Fraction */}
                    <label>
                     <p> Validation Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={taskData.task_builder.task_states.validationFraction}
                        onChange={(e) => {
                          setTaskData((prevTaskData) => ({
                            ...prevTaskData,
                            task_builder: {
                              ...prevTaskData.task_builder,
                              task_states: {
                                ...prevTaskData.task_builder.task_states,
                                dataSplitPerformed: false,
                                validationFraction: parseFloat(e.target.value),
                              },
                            },
                          }));
                        }}
                      />
                      {taskData.task_builder.task_states.validationFraction}
                    </label>

                    {/* Slider input for Test Fraction */}
                    <label>
                      <p>Test Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={taskData.task_builder.task_states.testFraction}
                        onChange={(e) => {
                          setTaskData((prevTaskData) => ({
                            ...prevTaskData,
                            task_builder: {
                              ...prevTaskData.task_builder,
                              task_states: {
                                ...prevTaskData.task_builder.task_states,
                                dataSplitPerformed: false,
                                testFraction: parseFloat(e.target.value),
                              },
                            },
                          }));
                        }}
                      />
                      {taskData.task_builder.task_states.testFraction}
                    </label>

                
                <Button 
                  onClick={() => handleDataSplit(index)} 
                  disabled={taskData.task_builder.task_states.dataSplitPerformed || loading}
                  variant="contained" 
                  color="primary"
                  style={{ marginTop: '20px' }}
                >
                  {loading ? 'Processing, please wait...' : 'Perform Data Split'}
                </Button>

                {taskData.task_builder.task_states.dataSplitPerformed && (
                  <Typography variant="body2" component="p">
                    <b>Archive Path: </b>{taskData.task_builder.task_states.archivePath}
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
