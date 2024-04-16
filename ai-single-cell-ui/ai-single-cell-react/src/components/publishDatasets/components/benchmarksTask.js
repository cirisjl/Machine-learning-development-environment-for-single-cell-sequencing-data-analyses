import React, {useEffect, useState} from 'react';
import {CELERY_BACKEND_API} from '../../../constants/declarations';
import {  ScaleLoader } from 'react-spinners';
import AlertMessageComponent from './alertMessageComponent';
import BenchmarksPlots from './benchmarksPlots';
import useWebSocket from '../../MyData/MyTasks/useWebSocket';
import { Typography,Paper, Grid, Card, CardContent,CardHeader } from '@mui/material';


function BenchmarksTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [loading, setLoading] = useState(false);
  const [ message, setMessage ] = useState('');
  const [ isError, setIsError ] = useState(false);
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [wsLogs, setWsLogs] = useState('');
  const [currentStatus, setCurrentStatus] = useState(null); // Set to null initially
  const [taskId, setTaskId] = useState('');
  const [celeryTaskResults, setCeleryTaskResults] = useState({});

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
    setActiveTask(6);
  };

  
  const handleStatusMessage = (event) => {
    try {
      const data = JSON.parse(event.data);
      console.log("Task Response");
      console.log(data);
      if (data.task_status) {
        setCurrentStatus(data.task_status);
        if(data.task_status === "SUCCESS" || data.task_status === "FAILURE"){
          setCeleryTaskResults(data);
        }
      }
    } catch (error) {
      console.error("Error parsing status message:", error);
    }
  };

  const handleLogMessage = (event) => {
    setWsLogs(event.data);
    // Auto-scroll to the bottom of the logs
    const logsElement = document.getElementById("_live_logs");
    if (logsElement) {
      logsElement.scrollTop = logsElement.scrollHeight;
    }
  };

      // A utility function to safely sanitize logs before using dangerouslySetInnerHTML
      const createMarkup = (logs) => {
        return { __html: logs };
      };

      const { closeWebSockets } = useWebSocket(taskId, handleStatusMessage, handleLogMessage);


  useEffect(() => {
    if(taskData.benchmarks.status !== 'completed') {

    if (taskData.task_builder && taskData.task_builder.selectedDatasets) {
      setLoading(true);

      const selectedDatasets  = taskData.task_builder.selectedDatasets;

      let body = Object.entries(selectedDatasets).map(([key, dataset], index) => ({
        benchmarksId: dataset.taskType.value + "-" + dataset.Id,
        datasetId: dataset.Id,
        userID: dataset.Owner,
        task_type: dataset.taskType.label,
        adata_path: dataset.adata_path,
        label: dataset.taskLabel.label || '',
      }));
      const postBody = body[0];

      fetch(`${CELERY_BACKEND_API}/api/benchmarks/create`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(postBody),
      })
        .then((response) => response.json())
        .then((data) => {
          const taskId = data.task_id;
          setTaskId(taskId);
        })
        .catch((error) => {
          console.error('Error during API call:', error);
          setMessage(`Error during API call: ${error}`);
          setHasMessage(true);
          setIsError(true);
          setLoading(false);
        })
    }
  }
  }, []);

  useEffect(() => {
    if(currentStatus === "SUCCESS" || currentStatus === "FAILURE") {
      closeWebSockets(); // Close WebSockets when task is done
      if(currentStatus === "SUCCESS") {
        // fetchProcessResults(celeryTaskResults.task_result.process_ids);
      }
      setLoading(false);
      setHasMessage(true);
      setMessage("Benchmarks task Success or failed");
      setIsError(true);
    }
  }, [currentStatus]); // Empty dependency array ensures this runs on mount and unmount only



  return (
    <div className='benchmarks-task'>
      
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}

      <Grid item xs={12} sx={{paddingTop: '10px'}}>
            <Card raised>
              <CardHeader title="Live Logs" />
              <CardContent>
                <Paper 
                  sx={{ 
                    maxHeight: 300, 
                    overflow: 'auto', 
                    '&::-webkit-scrollbar': { width: '0.4em' },
                    '&::-webkit-scrollbar-thumb': { 
                      backgroundColor: 'rgba(0,0,0,.1)',
                      borderRadius: '4px',
                    }
                  }} 
                  id="_live_logs"
                >
                  <Typography 
                    variant="body2" 
                    component="div" 
                    sx={{ fontFamily: 'monospace' }}
                    dangerouslySetInnerHTML={createMarkup(wsLogs || 'No Live logs...')}
                  />
                </Paper>
              </CardContent>
            </Card>
        </Grid>

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
