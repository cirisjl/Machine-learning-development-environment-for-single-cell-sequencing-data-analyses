import React from 'react';
import { useState, useEffect, useRef  } from 'react';
import axios from 'axios';
import { CELERY_BACKEND_API, STORAGE} from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from './reactPlotly';
import {isUserAuth, getCookie} from '../../../utils/utilFunctions';
import QualityControlParameters from './qualityControlParameters';
import { Button, makeStyles } from '@material-ui/core';
import { useNavigate } from 'react-router-dom';
import AlertMessageComponent from './alertMessageComponent';
import ReactSelect from 'react-select';
import { v4 as uuid } from 'uuid';
import FormControl from '@material-ui/core/FormControl';
import FormLabel from '@material-ui/core/FormLabel';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Radio from '@material-ui/core/Radio';

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
  const webSocketInstance = useRef(null);

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

  const runQualityControl = async() => {
    console.log(values);

    setLoading(true);

    const unique_id = uuid()

    
    // Establish WebSocket connection right away
    const websocketURL = `ws://${process.env.REACT_APP_HOST_URL}:5000/log/${unique_id}`;  
    webSocketInstance.current = new WebSocket(websocketURL);

    webSocketInstance.current.onopen = () => {
      console.log('WebSocket Connected');
    };

    webSocketInstance.current.onmessage = (event) => {
      const message = event.data;
      console.log(message);
      setWsLogs(message);
      let logs = document.getElementById("logs");
        let log_data = event.data;
        logs.innerHTML = log_data;
    };

    webSocketInstance.current.onerror = (error) => {
      console.error('WebSocket Error:', error);
    };

    webSocketInstance.current.onclose = () => {
      console.log('WebSocket Disconnected');
    };


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
      doublet_rate: values.doublet_rate,
      unique_id: unique_id
    }

    try {
      const response = await axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/run/quality_control`, inputRequest);
      const qualityControlResults = response.data;

      const ddl_assay_names = qualityControlResults[0].ddl_assay_names;

      if(!ddl_assay_names) {
        // Update the qc_results state with the quality control results
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            qc_results: qualityControlResults,
          },
        }));
        setLoading(false);
      } else {

        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          quality_control: {
            ...prevTaskData.quality_control,
            seurat_meta: {
              ...prevTaskData.quality_control.seurat_meta,
              default_assay: qualityControlResults[0].default_assay,
              assay_names: qualityControlResults[0].assay_names,
              file: qualityControlResults[0].inputfile,
              displayAssayNames: ddl_assay_names
            },
          }
        }));
        setLoading(false);
        return;
      }
    } catch (error) {
      console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
      setLoading(false);
      setHasMessage(true);
      setMessage("Failed to execute the quality control task.")
      setIsError(true);
    }
  };

  useEffect(() => {
    // This cleanup function will be called on component unmount
    return () => {
      if (webSocketInstance.current) {
        webSocketInstance.current.close();
        console.log('WebSocket Disconnected');
      }
    };
  }, []); // Empty dependency array ensures this runs on mount and unmount only


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

    const unique_id = uuid();

    if(!webSocketInstance.current) {
     // Establish WebSocket connection right away
     const websocketURL = `ws://${process.env.REACT_APP_HOST_URL}:5000/log/${unique_id}`; 
     webSocketInstance.current = new WebSocket(websocketURL);
    }


    webSocketInstance.current.onopen = () => {
      console.log('WebSocket Connected');
    };

    webSocketInstance.current.onmessage = (event) => {
      const message = event.data;
      setWsLogs(message);
    };

    webSocketInstance.current.onerror = (error) => {
      console.error('WebSocket Error:', error);
    };

    webSocketInstance.current.onclose = () => {
      console.log('WebSocket Disconnected');
    };

    let inputRequest = {
      fileDetails: taskData.quality_control.seurat_meta.file,
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
      assay: taskData.quality_control.selectedAssayName,
      unique_id: unique_id
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

      <div class="flex items-center py-2 px-3">
        <div
          id="logs"
          class="block p-2.5 w-full text-sm text-gray-900 bg-gray-50 rounded-lg border border-gray-300 focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
        >
          reading logs...
        </div>
      </div>

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
                        <FormLabel id="demo-row-radio-buttons-group-label">Dimension</FormLabel>
                        <RadioGroup
                          row
                          aria-labelledby="demo-row-radio-buttons-group-label"
                          name="row-radio-buttons-group"
                          value={plotDimension}
                          onChange={(event) => setPlotDimension(event.target.value)}
                        >
                          <FormControlLabel value="2D" control={<Radio />} label="2D" />
                          <FormControlLabel value="3D" control={<Radio />} label="3D" />
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
                  {/* {result.highest_expr_genes_plot && (
                    <>
                      <h2>Highest expression Genes Plot</h2>
                      <ReactPlotly plot_data={result.highest_expr_genes_plot} />
                    </>
                  )} */}
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

{/* {wsLogs && (
         <div class="flex items-center py-2 px-3">
         <div
           id="logs"
           class="block p-2.5 w-full text-sm text-gray-900 bg-gray-50 rounded-lg border border-gray-300 focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
         >
           reading logs...
         </div>
       </div>
)} */}

    </div>
  );
}

export default QualityControlTaskComponent;
