import React from 'react';
import { useState, useEffect } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API, STORAGE} from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from './reactPlotly';
import {isUserAuth, getCookie} from '../../../utils/utilFunctions';
import QualityControlParameters from './qualityControlParameters';
import {Button } from '@material-ui/core';
import { useNavigate } from 'react-router-dom';

const defaultValues = {
  min_genes: 200,
  max_genes: 20000, // No limit
  min_cells: 2,
  target_sum: 1e4,
  n_top_genes: 2000,
  n_neighbors: 15,
  n_pcs: 0, // None
  resolution: 1,
  regress_cell_cycle: false,
  use_default: true,
  doublet_rate: 0.08
};

function QualityControlTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {
  
  const [loading, setLoading] = useState(false);
  const [values, setValues] = useState(defaultValues);
  const navigate = useNavigate();

  const runQualityControl = async() => {

    setLoading(true);

    let file_paths = taskData.quality_control.file_paths;
    let pathToUse;
  
    if (file_paths.length > 1) {
        const firstFilePath = file_paths[0]; 
        pathToUse = firstFilePath.substring(0, firstFilePath.lastIndexOf('/'));
    } else if (file_paths.length === 1) {
        // Only one file, use its complete path
        pathToUse = file_paths[0];
    }

    let inputRequest = {
      fileDetails: pathToUse,
      min_genes : values.min_genes,
      max_genes : values.max_genes,
      min_cells : values.min_cells,
      target_sum : values.target_sum,
      n_top_genes : values.n_top_genes,
      n_neighbors : values.n_neighbors,
      n_pcs : values.n_pcs,
      resolution : values.resolution,
      regress_cell_cycle : values.regress_cell_cycle,
      use_default : values.use_default,
      doublet_rate: values.doublet_rate
    }

    try {
      const response = await axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/run/quality_control`, inputRequest);
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
      console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
      setLoading(false);
    }

  };


  useEffect(() => {
    isUserAuth(getCookie('jwtToken'))
    .then((authData) => {
      if (authData.isAdmin) {
        if(taskData.quality_control.status !== 'completed') {
          let files = taskData.upload.files;
          let newDirectoryPath = taskData.upload.newDirectoryPath;

          let inputFiles = [];
          for (let file of files) {
            let path = STORAGE + "/" + authData.username + "/" + newDirectoryPath + "/" + file;
            inputFiles.push(path);
          }
          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            quality_control: {
              ...prevTaskData.quality_control,
              file_paths: inputFiles,
              token: authData.username
            },
          }));   
        }
      }  else {
        console.warn("Unauthorized - you must be an admin to access this page");
        navigate("/accessDenied");
      }
    })
    .catch((error) => {
      console.error(error);
    });
  }, []); 

  // useEffect(() => {
  //   if(taskData.quality_control.status !== 'completed') {
  //   // Check the fileMappings state to determine if quality control should be run
  //   if (taskData.validation.fileMappings.length > 0) {
  //     setLoading(true);
  //     try {
  //       // Make an API call to run quality control
  //       const runQualityControl = async () => {
  //         try {
  //           // Make an API call to run quality control
  //           const response = await axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/run/quality_control`, taskData.validation.fileMappings);

  //           const qualityControlResults = response.data;
          
  //           // Update the qc_results state with the quality control results
  //           setTaskData((prevTaskData) => ({
  //             ...prevTaskData,
  //             quality_control: {
  //               ...prevTaskData.quality_control,
  //               qc_results: qualityControlResults,
  //             },
  //           }));
  //           setLoading(false);

  //         } catch (error) {
  //           console.error('Error running quality control:', error);
  //           setLoading(false);
  //         }
  //       };

  //       runQualityControl();
  //     } catch (error) {
  //       console.error('Error checking fileMappings:', error);
  //       setLoading(false);
  //     }
  //   } else {
  //     // No files to run quality control on
  //   }
  // }
  // }, []); // Empty dependency array ensures this effect runs only once when the component mounts

  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

  const handleTaskCompletion = () => {
      // Update the fileMappings state with the new list
      setTaskData((prevTaskData) => ({
        ...prevTaskData,
        quality_control: {
          ...prevTaskData.quality_control,
          status: 'completed'
        },
      }));

    // After Task 3 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      2: true, // Mark Task 3 as completed
    }));

    //The current task is finished, so make the next task active
    setActiveTask(3);
  };

  return (
    <div className='quality-control-task'>
      <div>
        <QualityControlParameters values={values} setValues={setValues} defaultValues={defaultValues}/>
        <div style={{ display: 'flex', justifyContent: 'center', marginTop: '20px' }}>
          <Button onClick={runQualityControl} variant="contained" color="primary">
            Run Quality Control
          </Button>
        </div>
      </div>

      {loading ? (
        <div className="spinner-container">
          <ScaleLoader color="#36d7b7" loading={loading} />
        </div>
      ) : (
      <div>
        <div className="App">
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
                  {/* {result.highest_expr_genes_plot && (
                    <>
                      <h2>Highest expression Genes Plot</h2>
                      <ReactPlotly plot_data={result.highest_expr_genes_plot} />
                    </>
                  )} */}
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
