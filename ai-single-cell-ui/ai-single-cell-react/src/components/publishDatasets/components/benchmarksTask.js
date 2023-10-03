import React from 'react';

function BenchmarksTaskComponent({ setTaskStatus }) {
  const handleTask1Completion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      1: true, // Mark Task 1 as completed
    }));
  };

  return (
    <div>
      {/* Task 1 content here */}
      <button onClick={handleTask1Completion}>GetMetaDataComponent button</button>
    </div>
  );
}

export default BenchmarksTaskComponent;
