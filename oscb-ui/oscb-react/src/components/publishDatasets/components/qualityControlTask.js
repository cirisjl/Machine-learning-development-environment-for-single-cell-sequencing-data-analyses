import React from 'react';
import { useState, useEffect, useRef  } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API, STORAGE, defaultValues} from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from './reactPlotly';
import {isUserAuth, getCookie, plotUmapObs} from '../../../utils/utilFunctions';
import QualityControlParameters from './qualityControlParameters';
import { Button, makeStyles } from '@material-ui/core';
import { useNavigate } from 'react-router-dom';
import AlertMessageComponent from './alertMessageComponent';
import ReactSelect from 'react-select';
import FormControl from '@mui/material/FormControl';
import Radio from '@mui/material/Radio';
import RadioGroup from '@mui/material/RadioGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import useWebSocket from '../../MyData/MyTasks/useWebSocket';
import LogComponent from '../../common_components/liveLogs';
import {Select, MenuItem, InputLabel } from '@mui/material';

function QualityControlTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [plotDimension, setPlotDimension] = useState('2D');
 
  const useStyles = makeStyles((theme) => ({
    customButton: {
      backgroundColor: '#007BFF',
      color: '#fff', // Change text color to white for better contrast
      textTransform: 'capitalize', // Capitalize the first letter of each word
      '&:hover': {
        backgroundColor: theme.palette.primary.dark, // Darken button on hover, adjust color as needed
      },
    },
  }));

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [ isError, setIsError ] = useState(false);
  const classes = useStyles(); // Use the custom styles
  const [loading, setLoading] = useState(false);
  const [values, setValues] = useState(taskData.quality_control.qc_params);
  const navigate = useNavigate();
  const [wsLogs, setWsLogs] = useState('');
  const [currentStatus, setCurrentStatus] = useState(null); // Set to null initially
  const [jobId, setjobId] = useState('');
  const [celeryTaskResults, setCeleryTaskResults] = useState({});
  const [clusteringPlotType, setClusteringPlotType] = useState('');
  const [plotData, setPlotData] = useState(null); // State to store the fetched plot data
  const [loadingPlot, setLoadingPlot] = useState(false); // State to handle loading spinner


  const fetchPlotData = async (plotType, cell_metadata, umap, umap_3d) => {
      setLoadingPlot(true); // Set loading to true before making the API call
  
      const selectedCellType = null
  
      try {

        let data = getUmapPlotData(cell_metadata, umap, umap_3d, plotType, selectedCellType)
  
        if(data.umap_plot && data.umap_plot_3d) {
          setPlotData({umap_plot: data.umap_plot, umap_plot_3d: data.umap_plot_3d});
        }
      } catch (error) {
        console.error('Error fetching plot data:', error);
        alert(`Error fetching plot data: ${error}`);
      } finally {
        setLoadingPlot(false); 
      }
    };

  const getUmapPlotData = (cell_metadata, umap, umap_3d, clustering_plot_type, annotation) => {
    const umap_plot = plotUmapObs(cell_metadata, umap, clustering_plot_type, [], annotation, 2);
    const umap_plot_3d = plotUmapObs(cell_metadata, umap_3d, clustering_plot_type, [], annotation, 3);
  
    return { umap_plot, umap_plot_3d };
  };

  const extractDir =  (inputFile) => {
    const fileLocParts = inputFile.split('/');
    fileLocParts.pop(); // Remove the file name from the array
    const output = fileLocParts.join('/'); // Join the remaining parts with '/'
    return output;
};

const handleStatusMessage = (event) => {
  try {
    const data = JSON.parse(event.data);
    if (data.task_status) {
      setCurrentStatus(data.task_status);
      if(data.task_status === "SUCCESS" || data.task_status === "FAILURE"){
        setCeleryTaskResults(data);
      }
    }
  } catch (error) {
    setLoading(false);
    console.error("Error parsing status message:", error);
  }
};

const handleLogMessage = (event) => {
  setWsLogs((prevLogs) => prevLogs + event.data);
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

  const { closeWebSockets } = useWebSocket(jobId, handleStatusMessage, handleLogMessage, setLoading);

  const fetchProcessResults = async (processIds) => {
    if (!processIds.length) return;

    try {
      const response = await axios.post(`${CELERY_BACKEND_API}/getPreProcessResults`, { process_ids: processIds });
      console.log('Process Results:', response.data);
      setTaskData((prevTaskData) => ({
        ...prevTaskData,
        quality_control: {
          ...prevTaskData.quality_control,
          qc_results: response.data,
        },
      }));
      setLoading(false);
    } catch (error) {
      console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
      setLoading(false);
      setHasMessage(true);
      setMessage("Failed to retrieve pre processed results from MongoDB");
      setIsError(true);
    }
  };

  const runQualityControl = async() => {
    setLoading(true);

    try {
      let file_paths = taskData.quality_control.file_paths;
      let pathToUse;
    
      if (file_paths.length > 1) {
          const firstFilePath = file_paths[0]; 
          pathToUse = firstFilePath.substring(0, firstFilePath.lastIndexOf('/')+1);
      } else if (file_paths.length === 1) {
          // Only one file, use its complete path
          pathToUse = file_paths[0];
      }

      let inputRequest = {
        dataset: taskData.upload.title,
        input: pathToUse,
        output: pathToUse + "/Results",
        userID: taskData.quality_control.token,
        qc_params : {
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
      }

 
      const response = await axios.post(`${CELERY_BACKEND_API}/tools/qc`, inputRequest);
      const taskInfo = response.data;
      const jobId = taskInfo.job_id;
      setjobId(jobId);

    } catch (error) {
      console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
      setLoading(false);
      setHasMessage(true);
      setMessage("Failed to execute the quality control task.")
      setIsError(true);
    }
  };

  useEffect(() => {
    if(currentStatus === "SUCCESS" || currentStatus === "FAILURE") {
      if(currentStatus === "SUCCESS" && celeryTaskResults.task_result.ddl_assay_names) {
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            seurat_meta: {
              ...prevTaskData.quality_control.seurat_meta,
              default_assay: celeryTaskResults.task_result.default_assay,
              assay_names: celeryTaskResults.task_result.assay_names,
              file: celeryTaskResults.task_result.inputfile,
              displayAssayNames: celeryTaskResults.task_result.ddl_assay_names
            },
          }
        }));
      }
      else if(currentStatus === "SUCCESS" && celeryTaskResults.task_result.process_ids) {
        fetchProcessResults(celeryTaskResults.task_result.process_ids);
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
             nCells: celeryTaskResults.task_result.nCells,
            seurat_meta: {
              ...prevTaskData.quality_control.seurat_meta,
              displayAssayNames: false
            },
          },
        }));
        setMessage("quality control task is Successful");
        setHasMessage(true);
        setIsError(false);
      } else if(currentStatus === "FAILURE"){
        setMessage("quality control task is Failed");
        setHasMessage(true);
        setIsError(true);
      }
      setLoading(false);
      closeWebSockets();
    }
  }, [currentStatus]); // Empty dependency array ensures this runs on mount and unmount only

  useEffect(() => {
    isUserAuth(getCookie('jwtToken'))
    .then((authData) => {
      if (authData.isAdmin) {
        if(taskData.quality_control.status !== 'completed') {
          let files = taskData.upload.files;
          let newDirectoryPath = taskData.upload.newDirectoryPath;

          let inputFiles = [];
          for (let file of files) {
            // let path = STORAGE + "/" + authData.username + "/" + newDirectoryPath + "/" + file;
            let path = STORAGE + "/Benchmarks/" + newDirectoryPath + "/" + file;
            inputFiles.push(path);
          }

          let shouldHideForSeurat = false;
          if(inputFiles.length === 1 && (inputFiles[0].toLowerCase().endsWith('h5seurat') || inputFiles[0].toLowerCase().endsWith('rds') || inputFiles[0].toLowerCase().endsWith('robj'))) {
            shouldHideForSeurat = true;
          } 

          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            quality_control: {
              ...prevTaskData.quality_control,
              // nCells: celeryTaskResults.task_result.nCells,
              file_paths: inputFiles,
              token: authData.username,
              shouldHideForSeurat: shouldHideForSeurat
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

  const handleAssaySelection = selectedOption => {
    setTaskData(prevTaskData => ({
        ...prevTaskData,
        quality_control: {
            ...prevTaskData.quality_control,
            selectedAssayName: (selectedOption ? selectedOption.value : null)
        },
    }));
};

const handleAssaySelectionSubmit = async () => {

  try {
    if(!taskData.quality_control.selectedAssayName) {
      setHasMessage(true);
      setMessage("Please select the default assay");
      setIsError(true);
      return;
    }

    setLoading(true);

    let inputRequest = {
      dataset: taskData.upload.title,
      input: taskData.quality_control.seurat_meta.file,
      output: extractDir(taskData.quality_control.seurat_meta.file) + "/Results",
      userID: taskData.quality_control.token,
      qc_params : {
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
        doublet_rate: values.doublet_rate,
        assay: taskData.quality_control.selectedAssayName
      } 
    }

  
      const response = await axios.post(`${CELERY_BACKEND_API}/tools/qc`, inputRequest);

      const taskInfo = response.data;

      const jobId = taskInfo.job_id;
      setjobId(jobId);

    } catch (error) {
      console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
      setLoading(false);
      setHasMessage(true);
      setMessage("Failed to execute the quality control task.")
      setIsError(true);
    }
};

  const handleTaskCompletion = () => {

      if(taskData.quality_control.qc_results.length > 0) {
        // Update the fileMappings state with the new list
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            status: 'completed',
            qc_params: values
          },
        }));

      // After Task 3 is successfully completed, update the task status
      setTaskStatus((prevTaskStatus) => ({
        ...prevTaskStatus,
        2: true, // Mark Task 3 as completed
      }));

      //The current task is finished, so make the next task active
      setActiveTask(3);
    } else {
      setHasMessage(true);
      setMessage("Run Quality Control Task Before moving to the next step.");
      setIsError(true);
    }
  };

  return (
    <div className='quality-control-task'>

      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}

      <div>
        <QualityControlParameters values={values} setValues={setValues} defaultValues={defaultValues} shouldHideForSeurat={taskData.quality_control.shouldHideForSeurat}/>
        <div style={{ display: 'flex', justifyContent: 'center', marginTop: '30px', marginLeft: '10px' }}>
        <Button onClick={runQualityControl} variant="contained" className={classes.customButton}>
          Run Quality Control
        </Button>
        </div>
      </div>

    <LogComponent wsLogs = {wsLogs}/>



      {taskData.quality_control.seurat_meta.displayAssayNames && (
            <div>
            <h3><b>File:</b> {taskData.quality_control.seurat_meta.file}</h3>
            <h3><b>Default Assay:</b> {taskData.quality_control.seurat_meta.default_assay}</h3>
            <p>Do you want to change the default assay?</p>
            <ReactSelect
            id="assaySelection"
            placeholder="Select an Assay"
            options={taskData.quality_control.seurat_meta.assay_names.map(name => ({ value: name, label: name }))}
            onChange={handleAssaySelection}
            />

      {taskData.quality_control.qc_results.length === 0 && (

            <div className='navigation-buttons'>
              <div className="previous">
                <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                  Previous
                </button>
              </div>
              <div className="next-upon-success">
                <button type="submit" className="btn btn-info button" onClick={handleAssaySelectionSubmit}>
                  Next
                </button>
              </div>
            </div>

      )}
        </div>
        )}

      {loading ? (
        <div className="spinner-container">
          <ScaleLoader color="#36d7b7" loading={loading} />
        </div>
      ) : (

      <div>
        <div align="center">
        {taskData.quality_control.qc_results &&
          taskData.quality_control.qc_results.map((result, index) => (
            <React.Fragment key={index}>
                  {result.umap_plot && (
                    <>
                    <h2>UMAP Plot</h2>
                    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', gap: '16px', marginTop: '16px' }}>
                      {/* Radio buttons */}
                      <div>
                        <label>
                          <input
                            type="radio"
                            value="2D"
                            checked={plotDimension === "2D"}
                            onChange={(e) => setPlotDimension(e.target.value)}
                          />
                          2D
                        </label>
                        <label style={{ marginLeft: '8px' }}>
                          <input
                            type="radio"
                            value="3D"
                            checked={plotDimension === "3D"}
                            onChange={(e) => setPlotDimension(e.target.value)}
                          />
                          3D
                        </label>
                      </div>

                      <select
                        value={clusteringPlotType}
                        onChange={(e) => {
                          const selectedPlotType = e.target.value;
                          setClusteringPlotType(selectedPlotType);
                          fetchPlotData(selectedPlotType, result.obs, result.umap, result.umap_3d);
                        }}
                        style={{
                          padding: '8px 12px',
                          borderRadius: '6px',
                          border: '1px solid #ccc',
                          backgroundColor: '#fff',
                          fontSize: '14px',
                          fontFamily: 'Arial, sans-serif',
                          outline: 'none',
                          transition: 'border-color 0.2s, box-shadow 0.2s',
                          boxShadow: '0 1px 2px rgba(0,0,0,0.1)',
                          cursor: 'pointer',
                          minWidth: '150px',
                        }}
                        onFocus={(e) => {
                          e.target.style.borderColor = '#1976d2';
                          e.target.style.boxShadow = '0 0 0 2px rgba(25, 118, 210, 0.2)';
                        }}
                        onBlur={(e) => {
                          e.target.style.borderColor = '#ccc';
                          e.target.style.boxShadow = '0 1px 2px rgba(0,0,0,0.1)';
                        }}
                      >
                        {Object.keys(result.cell_metadata).map((key) => (
                          <option key={key} value={key}>
                            {key}
                          </option>
                        ))}
                      </select>

                    </div>

                
                    {plotDimension === '2D' && (plotData?.umap_plot || result.umap_plot) && (
                      <>
                        <ReactPlotly plot_data={plotData?.umap_plot || result.umap_plot} />
                      </>
                    )}
                    {plotDimension === '3D' && (plotData?.umap_plot_3d || result.umap_plot_3d) && (
                      <>
                        <ReactPlotly plot_data={plotData?.umap_plot_3d || result.umap_plot_3d} />
                      </>
                    )}
                       
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
                      <h2>Highest Expression Genes Plot</h2>
                      <ReactPlotly plot_data={result.highest_expr_genes_plot} />
                    </>
                  )}
                </React.Fragment>
        ))}
        </div>
        {!taskData.quality_control.seurat_meta.displayAssayNames && (

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
        )}
      </div>
      )}

    </div>
  );
}

export default QualityControlTaskComponent;
