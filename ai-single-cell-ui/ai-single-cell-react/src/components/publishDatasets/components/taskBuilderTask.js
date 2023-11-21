import React , { useState, useEffect }from 'react';
import Select from 'react-select';
import { CELERY_BACKEND_API} from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  // const [data_path, setDataPath] = useState('');
  const [trainFraction, setTrainFraction] = useState(0.8);
  const [validationFraction, setValidationFraction] = useState(0.1);
  const [testFraction, setTestFraction] = useState(0.1);
  const [dataSplitPerformed, setDataSplitPerformed] = useState(false);
  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [archivePath, setArchivePath] = useState('');


  const handleDataSplit = async (index) => {
    try {

      const dataPath = taskData.quality_control.qc_results[index].adata_path || ''
      // User input for data split fractions
      const userData = {
        data: dataPath,
        train_fraction: trainFraction,
        validation_fraction: validationFraction,
        test_fraction: testFraction,
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
        setArchivePath(result.archive_path); // Set the archive path
        setDataSplitPerformed(true); // Set the state to indicate data split performed
      } else {
        const error = await response.json();
        console.error(error.error); // Handle the error
      }
    } catch (error) {
      console.error('Error:', error);
    }
  };

  const handleTaskCompletion = () => {
    const isValid =
            taskData.task_builder.task_type &&
            taskData.task_builder.task_label.length > 0 &&
            dataSplitPerformed
            
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
                        value={trainFraction}
                        onChange={(e) => {
                          setTrainFraction(parseFloat(e.target.value));
                          setDataSplitPerformed(false);
                        }}
                      />
                      {trainFraction}
                    </label>

                    {/* Slider input for Validation Fraction */}
                    <label>
                     <p> Validation Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={validationFraction}
                        onChange={(e) => {
                          setValidationFraction(parseFloat(e.target.value));
                          setDataSplitPerformed(false);
                        }}
                      />
                      {validationFraction}
                    </label>

                    {/* Slider input for Test Fraction */}
                    <label>
                      <p>Test Fraction:</p>
                      <input
                        type="range"
                        min={0}
                        max={1}
                        step={0.01}
                        value={testFraction}
                        onChange={(e) => {
                          setTestFraction(parseFloat(e.target.value));
                          setDataSplitPerformed(false);
                        }}
                      />
                      {testFraction}
                    </label>

                    {/* Button to perform data split */}
                    <button onClick={() => handleDataSplit(index)} disabled={dataSplitPerformed}>Perform Data Split</button>

                    {dataSplitPerformed && <p><b>Archive Path: </b>{archivePath}</p>}

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
