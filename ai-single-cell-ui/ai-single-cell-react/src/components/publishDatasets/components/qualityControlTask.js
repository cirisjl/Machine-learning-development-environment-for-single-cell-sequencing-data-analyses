import React from 'react';

function QualityControlTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {
  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      3: true, // Mark Task 1 as completed
    }));

    //The current task is finished, so make the next task active
    setActiveTask(4);
  };

  return (
    <div>
      {/* Task 1 content here */}
      <button onClick={handleTaskCompletion}>QualityControlTaskComponent button</button>
      <div className='previous'>
        <button type="submit" class="btn btn-info" onClick={() => setActiveTask(activeTask - 1)} >Previous</button>
      </div>
      <div className='next-upon-success'>
        <button type="submit" class="btn btn-info" onClick={handleTaskCompletion} >Next</button>
      </div>
    </div>
  );
}

export default QualityControlTaskComponent;
