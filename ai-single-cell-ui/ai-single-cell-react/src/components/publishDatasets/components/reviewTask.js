import React, { useState } from 'react';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import ReactPlotly from './reactPlotly';
import BenchmarksPlots from './benchmarksPlots';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';

function ReviewTaskComponent({setTaskStatus, taskData, setTaskData, setActiveTask, activeTask}) {
  
  const [sectionsVisibility, setSectionsVisibility] = useState({
    inputData: true,
    qcPlots: false,
    metadata: false,
    taskBuilder: false,
    benchmarks: false,
  });

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);

  const toggleSectionVisibility = (section) => {
    setSectionsVisibility((prevVisibility) => ({
      ...prevVisibility,
      [section]: !prevVisibility[section],
    }));
  };
  
  const handleTaskCompletion = () => {

    let {formData} = taskData.metadata;

    // populate required fields

    formData.Task = taskData.task_builder.task_type
    formData['Cell Count Estimate'] = taskData.quality_control.qc_results[0]?.metadata?.nCells || 0;

    // construct ID 
    const task_abbv = formData.Task.value;
    const species = formData.Species.value;
    const tissue = formData['Anatomical Entity'].label;
    const cellCount = formData['Cell Count Estimate'];
    const author = formData['Author'];
    const submissionDate = formData['Submission Date'];
    const year = submissionDate ? new Date(submissionDate).getFullYear().toString() : '';

    // Check if cellCount is greater than 1000
    const useCellCount = cellCount && parseInt(cellCount) > 1000;

    const constructedID = `${task_abbv}-${species}-${tissue}${useCellCount ? `-${cellCount}` : ''}-${author}-${year}`;

    formData.Id = constructedID;

    //Add label
    formData.Label = taskData.task_builder.task_label[0].label

    //Add genes and cells
    formData.Cells = JSON.stringify(taskData.quality_control.qc_results[0]?.metadata?.cells);
    formData.Genes = JSON.stringify(taskData.quality_control.qc_results[0]?.metadata?.genes);

    //Add Datasplit metadata
    formData['Data Split'] = {
      "trainFraction": taskData.task_builder.task_states.trainFraction,
      "testFraction": taskData.task_builder.task_states.testFraction,
      "validationFraction": taskData.task_builder.task_states.validationFraction,
      "archivePath": taskData.task_builder.task_states.archivePath
    }

    formData['QC_Plots'] = {
      "scatter_plot": taskData.quality_control.qc_results[0]?.scatter_plot,
      "umap_plot": taskData.quality_control.qc_results[0]?.umap_plot,
      "violin_plot": taskData.quality_control.qc_results[0]?.violin_plot,
      "highest_expr_genes_plot": taskData.quality_control.qc_results[0]?.highest_expr_genes_plot
    }

    formData.Status = "Review"
    console.log(formData);

    axios.post(`${SERVER_URL}/mongoDB/api/submitDatasetMetadata`, formData)
    .then(response => {
      console.log(response);
      console.log('Form data submitted successfully');
      setMessage("Form data submitted successfully")
      setHasMessage(true)
    })
    .catch(error => {
      console.error('Error submitting form data:', error.response.data.error);
      setMessage('Error submitting form data:', error.response.data.error)
      setHasMessage(true)
    });

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
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} />}
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
          <ul>
          {Object.entries(taskData.metadata.formData).map(([key, value]) => (
              <li key={key}>
                <strong>{key}:</strong> {typeof value === 'object' ? value.label : value}
              </li>
            ))}
          </ul>
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
          <ul>
            <li><strong>Task Type:</strong> {typeof taskData.task_builder.task_type === 'object' ? taskData.task_builder.task_type.label : taskData.task_builder.task_type}</li>
            <li><strong>Task Labels:</strong>
            <ul>
              {Array.isArray(taskData.task_builder.task_label) &&
                taskData.task_builder.task_label.map((labelItem, index) => (
                  <li key={index}>
                    {labelItem && typeof labelItem === 'object' && (
                      <>
                        {labelItem.label}
                      </>
                    )}
                  </li>
                ))}
            </ul></li>
            <li><strong>Train Fraction:</strong> {taskData.task_builder.task_states.trainFraction}</li>
            <li><strong>Validation Fraction:</strong> {taskData.task_builder.task_states.validationFraction}</li>
            <li><strong>Test Fraction:</strong> {taskData.task_builder.task_states.testFraction}</li>
          </ul>
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
