import React, { useState } from 'react';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import ReactPlotly from './reactPlotly';
import BenchmarksPlots from './benchmarksPlots';

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
          <h3>List of Datasets Uploaded</h3>
          <ul>
            {taskData.validation.inputFiles.map((dataset, index) => (
              <li key={index}>{dataset}</li>
            ))}
          </ul>
        </div>
      </div>
      <div className='section'>
        <div className='section-heading' onClick={() => toggleSectionVisibility('qcPlots')}>
          <h3>Quality Control</h3>
          <span className="category-icon">
            <FontAwesomeIcon
              icon={sectionsVisibility.qcPlots ? faAngleDown : faAngleRight}
            />
          </span>
        </div>
        <div className='section-content' style={{ display: sectionsVisibility.qcPlots ? 'block' : 'none' }}>
          <p>Quality Control Plots</p>
          {taskData.quality_control.qc_results &&
            taskData.quality_control.qc_results.map((result, index) => (
              <React.Fragment key={index}>
                    {result.umap_plot && (
                      <>
                        <h2>UMAP Plot</h2>
                        <ReactPlotly plot_data={result.umap_plot} />
                      </>
                    )}
                    {result.violin_plot && (
                      <>
                        <h2>Violin Plot</h2>
                        <ReactPlotly plot_data={result.violin_plot} />
                      </>
                    )}
                    {result.scatter_plot && (
                      <>
                        <h2>Scatter Plot</h2>
                        <ReactPlotly plot_data={result.scatter_plot} />
                      </>
                    )}
                    {result.highest_expr_genes_plot && (
                      <>
                        <h2>Highest expression Genes Plot</h2>
                        <ReactPlotly plot_data={result.highest_expr_genes_plot} />
                      </>
                    )}
                  </React.Fragment>
          ))}
        </div>
      </div>
      <div className='section'>
        <div className='section-heading' onClick={() => toggleSectionVisibility('metadata')}>
          <h3>Metadata</h3>
          <span className="category-icon">
            <FontAwesomeIcon
              icon={sectionsVisibility.metadata ? faAngleDown : faAngleRight}
            />
          </span>
        </div>
        <div className='section-content' style={{ display: sectionsVisibility.metadata ? 'block' : 'none' }}>
          <h3>Metadata Information</h3>
          {/* <ul>
            {Object.entries(taskData.metadata.formData).map(([key, value]) => (
              <li key={key}>
                <strong>{key}:</strong> {value}
              </li>
            ))}
          </ul> */}
        </div>
      </div>
      <div className='section'>
        <div className='section-heading' onClick={() => toggleSectionVisibility('taskBuilder')}>
          <h3>Task Builder</h3>
          <span className="category-icon">
            <FontAwesomeIcon
              icon={sectionsVisibility.taskBuilder ? faAngleDown : faAngleRight}
            />
          </span>
        </div>
        <div className='section-content' style={{ display: sectionsVisibility.taskBuilder ? 'block' : 'none' }}>
          {/* <strong>Task Type:</strong> {taskData.task_builder.task_type}
          <strong>Task Labels:</strong> */}
          {/* <ul>
            {taskData.task_builder.task_label.map((labelItem, index) => (
              <li key={index}>
                {labelItem.label} - {labelItem.value}
              </li>
            ))}
          </ul> */}
          {/* <strong>Train Fraction:</strong> {taskData.task_builder.task_states.trainFraction}
          <strong>Validation Fraction:</strong> {taskData.task_builder.task_states.validationFraction}
          <strong>Test Fraction:</strong> {taskData.task_builder.task_states.testFraction} */}
        </div>
      </div>
      <div className='section'>
        <div className='section-heading' onClick={() => toggleSectionVisibility('benchmarks')}>
          <h3>Benchmarks</h3>
          <span className="category-icon">
            <FontAwesomeIcon
              icon={sectionsVisibility.benchmarks ? faAngleDown : faAngleRight}
            />
          </span>
        </div>
        <div className='section-content' style={{ display: sectionsVisibility.benchmarks ? 'block' : 'none' }}>
          <p>Benchmarks Plots</p>
          {taskData.benchmarks &&
              taskData.benchmarks.benchmarks_results &&
              taskData.benchmarks.benchmarks_results.map((result, index) => (
                <React.Fragment key={index}>
                  {/* Assuming BenchmarksPlot is a component that you want to render */}
                  <BenchmarksPlots
                    barPlot={result.bar_plot}
                    linePlot={result.line_plot}
                  />
              </React.Fragment>
          ))}
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
