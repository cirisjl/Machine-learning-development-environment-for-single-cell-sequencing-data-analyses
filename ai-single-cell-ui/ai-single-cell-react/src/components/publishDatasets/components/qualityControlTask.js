import React from 'react';
import { useState, useEffect } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API} from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from './reactPlotly';

function QualityControlTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {
  
  const [loading, setLoading] = useState(false);

  useEffect(() => {

    // Check the fileMappings state to determine if quality control should be run
    if (taskData.validation.fileMappings.length > 0) {
      setLoading(true);
      try {
        // Make an API call to run quality control
        const runQualityControl = async () => {
          try {
            // Make an API call to run quality control
            const response = await axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/run/quality_control`, taskData.validation.fileMappings);

            const qualityControlResults = response.data;
          
            // Update the qc_results state with the quality control results
            setTaskData((prevTaskData) => ({
              ...prevTaskData,
              quality_control: {
                ...prevTaskData.quality_control,
                qc_results: qualityControlResults,
              },
            }));
            setLoading(false);

          } catch (error) {
            console.error('Error running quality control:', error);
            setLoading(false);
          }
        };

        runQualityControl();
      } catch (error) {
        console.error('Error checking fileMappings:', error);
        setLoading(false);
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
    <div className='quality-control-task'>
      {/* Task 1 content here */}
      {loading ? (
        <div className="spinner-container">
          <ScaleLoader color="#36d7b7" loading={loading} />
        </div>
      ) : (
      <div>
        <div className="App">
        {/* {taskData.quality_control.qc_results &&
              taskData.quality_control.qc_results.map((result, index) => (
                result.umap_plot && <ReactPlotly plot_data={result.umap_plot} />
        ))} */}
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
      )}
    </div>
  );
}

export default QualityControlTaskComponent;
