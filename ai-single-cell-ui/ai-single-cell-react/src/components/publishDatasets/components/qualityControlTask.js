import React from 'react';
import { useState, useEffect, useRef  } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API, STORAGE, defaultValues, WEB_SOCKET_URL, SERVER_URL} from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from './reactPlotly';
import {isUserAuth, getCookie} from '../../../utils/utilFunctions';
import QualityControlParameters from './qualityControlParameters';
import { Button, makeStyles } from '@material-ui/core';
import { useNavigate } from 'react-router-dom';
import AlertMessageComponent from './alertMessageComponent';
import ReactSelect from 'react-select';
import { v4 as uuid } from 'uuid';
import FormControl from '@mui/material/FormControl';
import FormLabel from '@mui/material/FormLabel';
import Radio from '@mui/material/Radio';
import RadioGroup from '@mui/material/RadioGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import { Typography,Paper, Grid, Card, CardContent,CardHeader } from '@mui/material';
import useWebSocket from '../../MyData/MyTasks/useWebSocket';
import LogComponent from '../../common_components/liveLogs';

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
  const [values, setValues] = useState(defaultValues);
  const navigate = useNavigate();
  const [wsLogs, setWsLogs] = useState('');
  const [currentStatus, setCurrentStatus] = useState(null); // Set to null initially
  const [taskId, setTaskId] = useState('');
  const [celeryTaskResults, setCeleryTaskResults] = useState({});


  const extractDir =  (inputFile) => {
    const fileLocParts = inputFile.split('/');
    fileLocParts.pop(); // Remove the file name from the array
    const output = fileLocParts.join('/'); // Join the remaining parts with '/'
    return output;
};

const handleStatusMessage = (event) => {
  try {
    const data = JSON.parse(event.data);
    console.log("insde handle status message");
    console.log(event.data);
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
  setWsLogs(event.data);
  console.log("insde handle log message");
  console.log(event.data);
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

  const { closeWebSockets } = useWebSocket(taskId, handleStatusMessage, handleLogMessage, setLoading);

  const fetchProcessResults = async (processIds) => {
    if (!processIds.length) return;

    try {
      const response = await axios.post(`${SERVER_URL}/benchmarks/api/getPreProcessResults`, { processIds });
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
    console.log(values);

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

 
      const response = await axios.post(`${CELERY_BACKEND_API}/api/tools/qc`, inputRequest);
      const taskInfo = response.data;

      const taskId = taskInfo.task_id;
      setTaskId(taskId);

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
      if(currentStatus === "SUCCESS" && celeryTaskResults.task_result.process_ids) {
        fetchProcessResults(celeryTaskResults.task_result.process_ids);
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
            let path = STORAGE + "/" + authData.username + "/" + newDirectoryPath + "/" + file;
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

  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

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

    if(!taskData.quality_control.selectedAssayName) {
      setHasMessage(true);
      setMessage("Please select the default assay");
      setIsError(true);
      return;
    }
    console.log(values);

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
        doublet_rate: values.doublet_rate
      } 
    }

    try {
      const response = await axios.post(`${CELERY_BACKEND_API}/api/tools/qc`, inputRequest);

      const qualityControlResults = response.data;

        // Update the qc_results state with the quality control results
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            qc_results: qualityControlResults,
          },
        }));

        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            seurat_meta: {
              ...prevTaskData.quality_control.seurat_meta,
              displayAssayNames: false
            },
          }
        }));
        setLoading(false);

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
{/* 
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
        </Grid> */}

{/* <Grid item xs={12} sx={{ paddingTop: '10px' }}>
  <Card raised>
    <CardHeader title="Live Logs" />
    <CardContent>
      <Paper 
        sx={{ 
          maxHeight: 300, 
          overflow: 'auto', 
          backgroundColor: 'rgba(211,211,211,0.5)', // Setting the gray background color
          '&::-webkit-scrollbar': { 
            width: '0.4em' 
          },
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
          sx={{ 
            fontFamily: 'monospace', 
            color: 'rgba(0, 0, 0, 0.87)', // Ensuring the text color is darker for better readability
            padding: '8px' // Adding some padding for visual relief
          }}
          dangerouslySetInnerHTML={createMarkup(wsLogs || 'No Live logs...')}
        />
      </Paper>
    </CardContent>
  </Card>
</Grid> */}
{/* <div style={{ padding: '10px', width: '100%' }}>
  <div style={{ 
    boxShadow: '0px 0px 2px 1px rgba(0,0,0,0.2)', 
    backgroundColor: '#232323', // Dark background for card
    color: 'white', // White text for card header
    marginBottom: '10px',
  }}>
    <h2 style={{ 
      backgroundColor: '#323232', // Slightly lighter header for contrast
      margin: 0,
      padding: '10px',
    }}>Live Logs</h2>
  </div>
  <div style={{ 
    maxHeight: '300px', 
    overflow: 'auto', 
    backgroundColor: '#2D2D2D', // Dark background for content
    color: 'white', // White text color for logs
    padding: '8px', // Padding inside the Paper equivalent
    fontSize: '0.875rem', // Default font size for Typography variant="body2"
    fontFamily: 'monospace', // Monospace font for logs
  }}>
    <div dangerouslySetInnerHTML={createMarkup(wsLogs || 'No Live logs...')} />
  </div>
</div> */}
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
        <div className="App">
        {taskData.quality_control.qc_results &&
          taskData.quality_control.qc_results.map((result, index) => (
            <React.Fragment key={index}>
                  {result.umap_plot && (
                    <>
                      <h2>UMAP Plot</h2>

                      <FormControl>
                        {/* <FormLabel id="demo-row-radio-buttons-group-label">Dimension</FormLabel> */}
                        <RadioGroup
                          row
                          aria-labelledby="demo-row-radio-buttons-group-label"
                          name="row-radio-buttons-group"
                          value={plotDimension}
                          onChange={(event) => setPlotDimension(event.target.value)}
                        >
                          <FormControlLabel value="2D" control={<Radio color="secondary"/>} label="2D" />
                          <FormControlLabel value="3D" control={<Radio color="secondary"/>} label="3D" />
                        </RadioGroup>
                      </FormControl>

                      {plotDimension === '2D' && result.umap_plot && (
                          <>
                            <ReactPlotly plot_data={result.umap_plot} />
                          </>
                        )}
                        {plotDimension === '3D' && result.umap_plot_3d && (
                          <>
                            <ReactPlotly plot_data={result.umap_plot_3d}/>
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
                      <h2>Highest expression Genes Plot</h2>
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
