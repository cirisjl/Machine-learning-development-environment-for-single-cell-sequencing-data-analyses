import React, {useEffect, useState} from 'react';
import {CELERY_BACKEND_API} from '../../../constants/declarations';
import {  ScaleLoader } from 'react-spinners';
import AlertMessageComponent from './alertMessageComponent';
import BenchmarksPlots from './benchmarksPlots';
import { Card, CardContent, Typography} from '@material-ui/core';


function BenchmarksTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [loading, setLoading] = useState(false);
  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);

  const handleTaskCompletion = () => {

    // Update the fileMappings state with the new list
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      benchmarks: {
        ...prevTaskData.benchmarks,
        status: 'completed'
      },
    }));

    // After Task 6 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      6: true, // Mark Task 6 as completed
    }));
    //The current task is finished, so make the next task active
    setActiveTask(7);
  };

  useEffect(() => {
    if(taskData.benchmarks.status !== 'completed') {

    if (taskData.task_builder && taskData.task_builder.selectedDatasets) {
      setLoading(true);

      const selectedDatasets  = taskData.task_builder.selectedDatasets;

      const postBody = {
        task_type: 'Clustering',
        data: Object.entries(selectedDatasets).map(([key, dataset], index) => ({
          adata_path: dataset.adata_path,
          task_label: dataset.taskLabel.label || '',
          datasetId: dataset.Id
        })),
      };

      fetch(`${CELERY_BACKEND_API}/convert/publishDatasets/benchmarks`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(postBody),
      })
        .then((response) => response.json())
        .then((data) => {
          if (Array.isArray(data)) {
            const benchmarksResults = data;

            // Update the benchmarks section in taskData with the received results
            setTaskData((prevTaskData) => ({
              ...prevTaskData,
              benchmarks: {
                benchmarks_results: benchmarksResults,
              },
            }));
          } else {
            console.error('Invalid response format');
            setMessage('Invalid response format');
            setHasMessage(true);
          }
        })
        .catch((error) => {
          console.error('Error during API call:', error);
          setMessage(`Error during API call: ${error}`);
          setHasMessage(true);
        })
        .finally(() => {
          setLoading(false);
        });
    }
  }
  }, []);


  return (
    <div className='benchmarks-task'>
      
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} />}

      {loading ? (
              <div className="spinner-container">
                <ScaleLoader color="#36d7b7" loading={loading} />
              </div>
            ) : (
              <>
               {/* Iterate over benchmarks_results and call BenchmarksPlot */}
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
          </>
        )}
    </div>
  );
}

export default BenchmarksTaskComponent;
