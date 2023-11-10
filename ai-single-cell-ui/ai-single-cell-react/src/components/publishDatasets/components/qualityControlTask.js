import React from 'react';
import { useState, useEffect } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API} from '../../../constants/declarations';

function QualityControlTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {
  
  useEffect(() => {

    // Check the fileMappings state to determine if quality control should be run
    if (taskData.validation.fileMappings.length > 0) {
      try {
        // Make an API call to run quality control
        const runQualityControl = async () => {
          try {
            // Make an API call to run quality control
            const response = await axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/run/quality_control`, taskData.validation.fileMappings);

            const qualityControlResults = response.data.qc_results;
          
            // Update the qc_results state with the quality control results
            setTaskData((prevTaskData) => ({
              ...prevTaskData,
              quality_control: {
                ...prevTaskData.quality_control,
                qc_results: qualityControlResults,
              },
            }));

          } catch (error) {
            console.error('Error running quality control:', error);
          }
        };

        runQualityControl();
      } catch (error) {
        console.error('Error checking fileMappings:', error);
      }
    } else {
      // No files to run quality control on
    }
  }, []); // Empty dependency array ensures this effect runs only once when the component mounts

  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      3: true, // Mark Task 3 as completed
    }));

    //The current task is finished, so make the next task active
    setActiveTask(4);
  };

  return (
    <div>
      {/* Task 1 content here */}
      <button onClick={handleTaskCompletion}>QualityControlTaskComponent button</button>
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

export default QualityControlTaskComponent;
