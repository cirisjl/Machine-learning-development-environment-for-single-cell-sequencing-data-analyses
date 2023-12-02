import React , { useState, useEffect }from 'react';
import Select from 'react-select';
import { CELERY_BACKEND_API} from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';
import axios from 'axios';
import ReactPlotly from './reactPlotly';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [loading, setLoading] = useState(false);

  // useEffect(() => {
  //   if(taskData.task_builder.status !== 'completed') {
  //     let adata_paths = []
  //     taskData.quality_control.qc_results
  //             .filter((result) => result.metadata && result.metadata.cell_metadata_obs)
  //             .map((results, index) => (
  //         adata_paths.push(results.adata_path)
  //     ))

  //     axios.post(`${CELERY_BACKEND_API}/convert/api/getTablePlot`, adata_paths)
  //     .then(response => {
  //       const data = response.data;
       
  //       // Updating taskData state
  //       setTaskData(prevTaskData => ({
  //         ...prevTaskData,
  //         task_builder: {
  //           ...prevTaskData.task_builder,
  //           table_data: data,
  //         },
  //       }));
  //     })
  //     .catch(error => {
  //       setMessage('Error while getting the table data: ' + error.response.data.detail);
  //       setHasMessage(true);
  //     });
  //   }
  // }, []);

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
          {taskData.quality_control.qc_results && taskData.quality_control.qc_results.length > 0 && (
          <div className="metadata-section">
            {taskData.quality_control.qc_results
              .filter((result) => result.metadata && result.metadata.cell_metadata_obs)
              .map((metadata, index) => (
                <div key={index} className="metadata-item">
                  {/* <div>
                    <p>Please Analyse the table before choosing the label</p>
                    {taskData.task_builder.table_data.length > 0 && taskData.task_builder.table_data[index] && (
                      <ReactPlotly plot_data={taskData.task_builder.table_data[index]} />                   
                    )}
                  </div> */}
                  <div>
                    <label>
                      <p>Please Choose the Label:</p>
                      <Select
                        options={Object.keys(metadata.metadata.cell_metadata_obs).map((key) => ({
                          label: key,
                          value: key,
                        }))}
                        onChange={(selectedOption) => handleLabelChange(selectedOption, index)}
                        value={taskData.task_builder.task_label[index]}
                        />
                    </label>
                  </div>

                  <div className="split-parameters">
                    <h3> Data Split Parameters</h3>

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

                    {/* Button to perform data split */}
                    <button onClick={() => handleDataSplit(index)} disabled={taskData.task_builder.task_states.dataSplitPerformed || loading}>
                      {loading ? 'Processing, please wait...' : 'Perform Data Split'}
                    </button>

                    {taskData.task_builder.task_states.dataSplitPerformed && <p><b>Archive Path: </b>{taskData.task_builder.task_states.archivePath}</p>}

                  </div>
                </div>
            ))}
    </div>
    )}

      <div className='navigation-buttons'>
            <div className="previous">
              <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                Previous
              </button>
            </div>
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
