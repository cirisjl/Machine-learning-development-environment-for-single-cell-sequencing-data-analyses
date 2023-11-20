import React , { useState }from 'react';
import Select from 'react-select';
import { CELERY_BACKEND_API} from '../../../constants/declarations';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  // const [data_path, setDataPath] = useState('');
  const [trainFraction, setTrainFraction] = useState(0.8);
  const [validationFraction, setValidationFraction] = useState(0.1);
  const [testFraction, setTestFraction] = useState(0.1);


  const handleDataSplit = async (index) => {
    try {
      // User input for data split fractions
      const userData = {
        data: taskData.quality_control.qc_results[index].metadata.cell_metadata_obs.data_path,
        train_fraction: trainFraction,
        validation_fraction: validationFraction,
        test_fraction: testFraction,
      };

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
      } else {
        const error = await response.json();
        console.error(error.error); // Handle the error
      }
    } catch (error) {
      console.error('Error:', error);
    }
  };

  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      5: true, // Mark Task 5 as completed
    }));

    //The current task is finished, so make the next task active
    setActiveTask(6);
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

  // const handleLabelChange = (selectedOption) => {
  //   // Update the task_type in task_builder
  //   setTaskData((prevTaskData) => ({
  //     ...prevTaskData,
  //     task_builder: {
  //       ...prevTaskData.task_builder,
  //       task_label: selectedOption,
  //     },
  //   }));
  // };

  const handleLabelChange = (selectedOption, index) => {
    // Ensure that qc_results[index] and cell_metadata_obs are defined
    if (
      taskData.quality_control.qc_results[index] &&
      taskData.quality_control.qc_results[index].metadata &&
      taskData.quality_control.qc_results[index].metadata.cell_metadata_obs
    ) {
      const cellMetadataObs = taskData.quality_control.qc_results[index].metadata.cell_metadata_obs;
  
      // Now, you can access cell_metadata_obs and perform any actions needed
      console.log(cellMetadataObs);
  
      // Update the task_label in task_builder for the specific result
      setTaskData((prevTaskData) => {
        const updatedLabels = [...prevTaskData.task_builder.task_labels];
        updatedLabels[index] = selectedOption.value;
  
        return {
          ...prevTaskData,
          task_builder: {
            ...prevTaskData.task_builder,
            task_labels: updatedLabels,
          },
        };
      });
    }
  };
  

  return (
    <div className='task-builder-task'>
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
                        onChange={(e) => setTrainFraction(parseFloat(e.target.value))}
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
                        onChange={(e) => setValidationFraction(parseFloat(e.target.value))}
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
                        onChange={(e) => setTestFraction(parseFloat(e.target.value))}
                      />
                      {testFraction}
                    </label>

                    {/* Button to perform data split */}
                    <button onClick={() => handleDataSplit(index)}>Perform Data Split</button>
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
