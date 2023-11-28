import React, {useEffect, useState} from 'react';
import {CELERY_BACKEND_API} from '../../../constants/declarations';
import {  ScaleLoader } from 'react-spinners';
import AlertMessageComponent from './alertMessageComponent';
import BenchmarksPlots from './benchmarksPlots';

function BenchmarksTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [loading, setLoading] = useState(false);
  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);

  const handleTaskCompletion = () => {
    // After Task 6 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      6: true, // Mark Task 6 as completed
    }));
    //The current task is finished, so make the next task active
    setActiveTask(7);
  };

  useEffect(() => {
    if (taskData.quality_control && taskData.quality_control.qc_results) {
      setLoading(true);

      const { qc_results } = taskData.quality_control;

      const postBody = {
        task_type: taskData.task_builder.task_type.label,
        data: qc_results.map((result, index) => ({
          adata_path: result.adata_path,
          task_label: taskData.task_builder.task_label[index].label || '',
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
  }, [taskData]);


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
            taskData.benchmarks.benchmarks_results.map((result, index) => {
              // Assuming BenchmarksPlot is a component that you want to render
              return (
                <BenchmarksPlots
                  barPlot={result.bar_plot}
                  linePlot={result.line_plot}
                />
              );
            })}
            
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
