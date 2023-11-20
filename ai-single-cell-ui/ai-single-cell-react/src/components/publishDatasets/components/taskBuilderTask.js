import React from 'react';
import Select from 'react-select';

function TaskBuilderTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {
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

  const handleLabelChange = (selectedOption) => {
    // Update the task_type in task_builder
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      task_builder: {
        ...prevTaskData.task_builder,
        task_label: selectedOption,
      },
    }));
  };

  return (
    <div className='task-builder-task'>
      <div>
        <label>
          Please Choose the Task Type:
            <Select
              value={taskData.task_builder.task_type}
              options={taskData.metadata.taskOptions}
              onChange={handleTaskChange}
          />
        </label>
        <br />
          {taskData.task_builder.task_id && (
            <div>
              Task ID: {taskData.task_builder.task_id}
            </div>
          )}
        <br />
          {taskData.qc_results &&
            taskData.qc_results
              .filter((result) => result.metadata)
              .map((result, index) => (
                <div key={index}>
                  <label>
                    Please Choose the Label:
                    <Select
                      options={Object.keys(result.cell_metadata_obs).map((key) => ({
                        label: key,
                        value: key,
                      }))}
                       onChange={handleLabelChange}
                       value={taskData.task_builder.task_label}
                    />
                  </label>

                  <button>Perform Data Split</button>

                </div>
            ))}
    </div>
      {/* Task 1 content here */}
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
  );
}

export default TaskBuilderTaskComponent;
