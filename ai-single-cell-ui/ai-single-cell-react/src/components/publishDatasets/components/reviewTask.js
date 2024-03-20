import React, { useState } from 'react';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import ReactPlotly from './reactPlotly';
import BenchmarksPlots from './benchmarksPlots';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';
import AlertMessageComponent from './alertMessageComponent';
import { Card, CardContent, Typography} from '@material-ui/core';


function ReviewTaskComponent({setTaskStatus, taskData, setTaskData, setActiveTask, activeTask}) {
  
  const [sectionsVisibility, setSectionsVisibility] = useState({
    taskBuilder: false,
    benchmarks: false,
  });

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [ isError, setIsError ] = useState(false);

  const toggleSectionVisibility = (section) => {
    setSectionsVisibility((prevVisibility) => ({
      ...prevVisibility,
      [section]: !prevVisibility[section],
    }));
  };
  
  const handleTaskCompletion = () => {


    let documents = [];

    Object.entries(taskData.task_builder.selectedDatasets).forEach(([datasetId, datasetDetails]) => {
        let document = {};
        
          // Find the corresponding benchmark result for this dataset
        const benchmarkResult = taskData.benchmarks.benchmarks_results.find(result => result.datasetId === datasetId);

        if(benchmarkResult) {
          // construct ID for each document
          const taskType = datasetDetails.taskType.value; 
          const constructedID = `${taskType}-${datasetId}`;
        
          // populate the document with details from task builder and benchmarks results
          document = {
            Id: constructedID,
            TaskType: datasetDetails.taskType,
            TaskLabel: datasetDetails.taskLabel,
            DatasetId: datasetId,
            TrainFraction: datasetDetails.dataSplit.trainFraction,
            ValidationFraction: datasetDetails.dataSplit.validationFraction,
            TestFraction: datasetDetails.dataSplit.testFraction,
            ArchivePath: datasetDetails.dataSplit.archivePath,
            BenchmarkResults: {
              metrics: benchmarkResult.metrics,
              Scanpy: JSON.stringify(benchmarkResult.Scanpy),
              Seurat: JSON.stringify(benchmarkResult.Seurat),
              bar_plot: JSON.stringify(benchmarkResult.bar_plot),
              line_plot: JSON.stringify(benchmarkResult.line_plot)
            },            // ... Add any other necessary details
          };
        
          documents.push(document);
        } else {
          console.error(`No benchmark result found for datasetId: ${datasetId}`);
        }

    });
    

    axios.post(`${SERVER_URL}/mongoDB/api/submitTaskMetadata`, documents)
    .then(response => {
      const { data } = response;
  
      // Check if the response is an array (indicating multiple documents were processed)
      if (Array.isArray(data)) {
        const allSuccess = data.every(result => result.status === 'success');
        const messages = data.map(result => `${result.Id}: ${result.message}`).join('\n');
  
        if (allSuccess) {
          console.log('All documents submitted successfully:', messages);
          setMessage("All documents submitted successfully.");
        } else {
          console.warn('Some documents were not submitted successfully:', messages);
          setMessage("Some documents encountered errors. Check the console for details.");
          setIsError(true);
        }
      } else {
        // Handle single document submission response
        console.log('Document submitted successfully:', data.message);
        setMessage(data.message);
      }
  
      setHasMessage(true);
    })
    .catch(error => {
      console.error('Error submitting documents:', error.response ? error.response.data.error : error.message);
      setMessage(`Error submitting documents: ${error.response ? error.response.data.error : error.message}`);
      setHasMessage(true);
      setIsError(true);
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
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}
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
          {Object.entries(taskData.task_builder.selectedDatasets).map(([datasetId, datasetDetails], index) => (
            <Card key={index} className="benchmarks-results">
              <CardContent>
              <div key={datasetId} className="dataset-details">
                <Typography variant="body2">Dataset {index + 1}: {datasetDetails.Id} 
                  <ul>
                    <li><strong>Task Type:</strong> {datasetDetails.taskType ? datasetDetails.taskType.label : 'Not set'}</li>
                    <li><strong>Task Label:</strong> {datasetDetails.taskLabel ? datasetDetails.taskLabel.label : 'Not set'}</li>
                    <li><strong>Train Fraction:</strong> {datasetDetails.dataSplit.trainFraction}</li>
                    <li><strong>Validation Fraction:</strong> {datasetDetails.dataSplit.validationFraction}</li>
                    <li><strong>Test Fraction:</strong> {datasetDetails.dataSplit.testFraction}</li>
                  </ul>
                </Typography>
              </div>
              </CardContent>
            </Card>
          ))}
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
          {taskData.benchmarks &&
              taskData.benchmarks.benchmarks_results &&
              taskData.benchmarks.benchmarks_results.map((result, index) => (
                <React.Fragment key={index}>
                  <Card key={index} className="benchmarks-results">
                    <CardContent>
                      <Typography variant="body2">Benchmark Results for {result.datasetId}</Typography>
                      <BenchmarksPlots
                        barPlot={result.bar_plot}
                        linePlot={result.line_plot}
                      />
                    </CardContent>
                  </Card>
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
