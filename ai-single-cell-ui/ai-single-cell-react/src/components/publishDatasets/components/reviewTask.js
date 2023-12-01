import React, { useState } from 'react';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';

function ReviewTaskComponent({setTaskStatus, taskData, setTaskData, setActiveTask, activeTask}) {
  
  const [sectionsVisibility, setSectionsVisibility] = useState({
    inputData: true,
    qcPlots: true,
    metadata: true,
    taskBuilder: true,
    benchmarks: true,
  });

  const toggleSectionVisibility = (section) => {
    setSectionsVisibility((prevVisibility) => ({
      ...prevVisibility,
      [section]: !prevVisibility[section],
    }));
  };
  
  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 7 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      7: true, // Mark Task 7 as completed
    }));

    console.log("All tasks completed");
    console.log(taskData);
  };

  return (
    <div className='review-task'>

      <div className='section'>
        <div className='section-heading' onClick={() => toggleSectionVisibility('inputData')}>
          <h3>Input Data</h3>
          <span className="category-icon">
            <FontAwesomeIcon
              icon={sectionsVisibility.inputData ? faAngleDown : faAngleRight}
            />
          </span>
        </div>
        <div className='section-content' style={{ display: sectionsVisibility.inputData ? 'block' : 'none' }}>
          <h3> Content of input section</h3>
        </div>
      </div>
      
      <div className='qc-plots-section'>
        {/* Add button to toggle section visibility */}
        <button onClick={() => toggleSectionVisibility('qcPlots')}>Quality Control</button>

        {/* Render sections based on visibility state */}
        <div style={{ display: sectionsVisibility.qcPlots ? 'block' : 'none' }}>
          
        </div>
      </div>
      <div className='metadata-section'>
        {/* Add button to toggle section visibility */}
        <button onClick={() => toggleSectionVisibility('metadata')}>Metadata</button>

        {/* Render sections based on visibility state */}
        <div style={{ display: sectionsVisibility.metadata ? 'block' : 'none' }}>
          
        </div>
      </div>
      <div className='task-builder-section'>
        {/* Add button to toggle section visibility */}
        <button onClick={() => toggleSectionVisibility('taskBuilder')}>Task Builder</button>

        {/* Render sections based on visibility state */}
        <div style={{ display: sectionsVisibility.taskBuilder ? 'block' : 'none' }}>
          
        </div>
      </div>
      <div className='benchmarks-section'>
        {/* Add button to toggle section visibility */}
        <button onClick={() => toggleSectionVisibility('benchmarks')}>Benchmarks</button>

        {/* Render sections based on visibility state */}
        <div style={{ display: sectionsVisibility.benchmarks ? 'block' : 'none' }}>
          
        </div>
      </div>
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

export default ReviewTaskComponent;
