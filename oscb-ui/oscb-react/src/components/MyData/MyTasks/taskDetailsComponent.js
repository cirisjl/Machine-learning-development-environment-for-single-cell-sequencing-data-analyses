import React, { useState, useEffect} from 'react';
import { useLocation,useNavigate } from 'react-router-dom';
import useWebSocket from './useWebSocket'; // Custom hook for WebSocket
import { 
  Container, Typography, Chip, Box, CircularProgress, Paper, Grid,TextField,Button, 
  Card, CardContent, Link, CardHeader 
} from '@mui/material';
import { green, red, yellow } from '@mui/material/colors';
import RightRail from '../../RightNavigation/rightRail';
import LogComponent from '../../common_components/liveLogs';
import axios from 'axios';
import { Octokit } from "@octokit/rest";
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import { ScaleLoader } from 'react-spinners';
import FormControl from '@mui/material/FormControl';
import Radio from '@mui/material/Radio';
import RadioGroup from '@mui/material/RadioGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import ReactPlotly from '../../publishDatasets/components/reactPlotly';
import { getCookie, plotUmapObs } from '../../../utils/utilFunctions';
//GitImports
import { CELERY_BACKEND_API, NODE_API_URL, owner, repo } from '../../../constants/declarations';
import {Select, MenuItem, InputLabel } from '@mui/material';


// Initialize Octokit with your GitHub personal access token
const octokit = new Octokit({ auth: process.env.REACT_APP_TOKEN });

let jwtToken = getCookie('jwtToken');


function StatusChip({ status }) {
  const getStatusColor = () => {
    switch (status?.toLowerCase()) { // Ensure status is defined
      case 'success': return green[500];
      case 'failure': return red[500];
      case 'started': return yellow[700];
      default: return yellow[700];
    }
  };

  return status ? (
    <Chip label={status.toUpperCase()} style={{ backgroundColor: getStatusColor(), color: '#fff' }} />
  ) : (
    <Chip label="In Progress" style={{ backgroundColor: yellow[700], color: '#fff' }} />
  );
}


function getFileNameFromURL(fileUrl){
  if (fileUrl) {
    try { 
      const filename = fileUrl.substring(fileUrl.lastIndexOf('/') + 1);
      return filename;
    } 
    catch (e) { 
      console.error(e); 
    }
  } else{
    return '';
  }
};


function downloadFile(fileUrl) {
  const apiUrl = `${NODE_API_URL}/download`;
  const pwd = "jobResults";

  if (fileUrl) {
    fetch(`${apiUrl}?fileUrl=${fileUrl}&authToken=${jwtToken}&pwd=${pwd}`)
      .then(response => {
        return response.blob();
      })
      .then(blob => {
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        const filename = fileUrl.substring(fileUrl.lastIndexOf('/') + 1);
        link.href = url;
        link.download = filename;

        document.body.appendChild(link);
        link.click();
        // Remove the link from the DOM
        // document.body.removeChild(link);
      })
      .catch(error => {
        console.error('Error downloading file:', error);
      });
  }
}


function TaskDetailsComponent() {
  const location = useLocation();
  const navigate = useNavigate();
  const { job_id, method, datasetURL, description, process, output, results, status } = location.state || {};
  const searchParams = new URLSearchParams(location.search);
  // const resultsPath = searchParams.get('results_path');
  // const taskTitle = searchParams.get('description');
  const [taskStatus, setTaskStatus] = useState(null); // Set to null initially
  const [taskOutput, setTaskOutput] = useState(null); // Set to null initially
  const [liveLogs, setLiveLogs] = useState('');
  const [loading, setLoading] = useState(true);
  const [toolResultsFromMongo, setToolResultsFromMongo] = useState([]);
  const [uName, setUName] = useState(null);
  const [uIat, setUIat] = useState(null);
  const [taskResult, setTaskResult] = useState("");
  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [ isError, setIsError ] = useState(false);
  const [plotDimension, setPlotDimension] = useState('2D');
  const [userComment, setUserComment] = useState(''); // State for user comment
  const [isSaving, setIsSaving] = useState(false); // State to indicate save operation
  const [isSent, setIsSent] = useState(false); // State to disable button after success
  const [commentSuccessMessage, setCommentSuccessMessage] = useState('');
  const [showErrorLog, setShowErrorLog] = useState(true); // State to show/hide the error log card
  let plotLoaded = false;
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
    // A utility function to safely sanitize logs before using dangerouslySetInnerHTML
    const createMarkup = (logs) => {
      return { __html: logs };
    };

    const cardStyle = {
      height: '100%', // Makes the cards take the full height of their container
      display: 'flex', // Allows child items to be flex items
      flexDirection: 'column', // Stacks child items vertically
    };
  
    const cardContentStyle = {
      flexGrow: 1, // Allows the content to expand and fill the space
      overflow: 'auto' // Adds scroll for overflow content
    };

    const fetchProcessResults = async (processIds) => {
      if (!processIds.length) return;
  
      try {
        const response = await axios.post(`${CELERY_BACKEND_API}/getPreProcessResults`, { process_ids: processIds });
        console.log('Process Results:', response.data);
        setToolResultsFromMongo(response.data);
        setLoading(false);
      } catch (error) {
        console.error('There was a problem with the axios operation:', error.response ? error.response.data : error.message);
        setLoading(false);
        setHasMessage(true);
        setMessage("Failed to retrieve pre processed results from MongoDB");
        setIsError(true);
      }
    };
    

  const handleStatusMessage = (event) => {
    try {
      const data = JSON.parse(event.data);
      if (status?.toLowerCase() === "success" || status?.toLowerCase() === "failure"){
        setTaskStatus(status);
        if (status?.toLowerCase() === "success" && !plotLoaded) {
          // if (results.process_ids && (process === "Quality Control" || process === "Normalization" || process === "Visualization")) {
          if (results.process_ids) {
            console.log("results: ", results)
            setTaskResult(results);
            fetchProcessResults(results.process_ids);
            plotLoaded = true;
          } else {
            setLoading(false);
          }
        }       
      } else if (data.task_status) {
        setTaskStatus(data.task_status);
        if (data.task_status?.toLowerCase() === "success" || data.task_status?.toLowerCase() === "failure"){
          if (data.task_status?.toLowerCase() === "success" && !plotLoaded) {
            if(data.task_result.process_ids) {
              console.log("data.task_result: ", data.task_result)
              fetchProcessResults(data.task_result.process_ids);
              plotLoaded = true;
            } else {
              setLoading(false);
            }

            if(data.task_result.output){
              setTaskOutput(data.task_result.output);
            }
          } else {
            setLoading(false);
          } 
        }
      }
    } catch (error) {
      setLoading(false);
      console.error("Error parsing status message:", error);
    }
  };

  useEffect(() => {
    async function fetchFiles() {
      if (status?.toLowerCase() === "success" && output) {
        setTaskOutput(output);
      }
      else {
        try {
          const taskInfoResponse = await fetch(`${CELERY_BACKEND_API}/task/${job_id}`);
          const taskInfoData = await taskInfoResponse.json();
          // console.log("taskInfoData", taskInfoData);
          setTaskResult(taskInfoData.task_result);
          // console.log("taskInfoData.task_result", taskInfoData.task_result);
          // if (taskInfoData.task_result.output) {
          //   setTaskOutput(taskInfoData.task_result.output);
          // }

          if (jwtToken) {
            fetch(NODE_API_URL + "/protected", { //to get username, id
              method: 'GET',
              credentials: 'include',
              headers: { 'Authorization': `Bearer ${jwtToken}` },
            })
              .then((response) => response.json())
              .then((data) => {

                if (data.authData !== null) {
                  // console.log("userdata: ", data.authData);
                  setUName(data.authData.username);
                  setUIat(data.authData.iat);
                }
              })
              .catch((error) => {
                console.error(error);
              })
          }
        } catch (error) {
          console.error('Error fetching task status:', error);
        }
      }
    }
    fetchFiles();
  }, [job_id]);


  const handleLogMessage = (event) => {
    setLiveLogs((prevLogs) => prevLogs + event.data);
    // Auto-scroll to the bottom of the logs
    const logsElement = document.getElementById("_live_logs");
    if (logsElement) {
      logsElement.scrollTop = logsElement.scrollHeight;
    }
  };
  
  const saveErrorLogData = async () => {
    try {
      const response = await axios.post(`${NODE_API_URL}/errorlogdata`, {
        name: uName,
        id: uIat,
        taskResult: taskResult,
        taskStatus: taskStatus,
        job_id: job_id,
        userComments: userComment
      });

      if (response.status === 200) {
        setCommentSuccessMessage('Feedback sent successfully.');
        setIsSaving(false);
        setIsSent(true);
        setShowErrorLog(false); // Hide the error log card
      } else {
        setCommentSuccessMessage('Failed to send feedback.');
        setIsSaving(false);
      }
    } catch (error) {
      console.error('Error saving comment:', error);
      setCommentSuccessMessage('Failed to send feedback.');
      setIsSaving(false);
    }
  };

  const createGitHubIssue = async () => {
    try {
      const response = await octokit.issues.create({
        //owner: 'SAYEERA', // Replace with your GitHub username
        //repo: 'issues-list', // Replace with your repository name
        owner: owner,
        repo: repo,
        title: `Issue for Job ID: ${job_id}`,
        body: `
            User Name:${uName}
            User ID: ${uIat}
            Job Result: ${taskResult}
            Job Status: ${taskStatus}
            Job ID: ${job_id}
            User Comments: ${userComment}
        `
      });
      console.log("Git response: ", response);
      if (response.status === 201) {
        console.log('GitHub issue created successfully:', response.data.html_url);
      } else {
        console.error('Failed to create GitHub issue:', response);
      }
    } catch (error) {
      console.error('Error creating GitHub issue:', error);
    }
  };

  const handleSaveComment = async () => {
    setIsSaving(true);
    setIsSent(false);
    await saveErrorLogData();
    await createGitHubIssue();
  };

    // Use the WebSocket hook
    useWebSocket(job_id, handleStatusMessage, handleLogMessage);

  return (

    <div className="task-details-container eighty-twenty-grid">

      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}

      <div className="main-content">
        <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
          <Box display="flex" justifyContent="center">
            <Typography variant="h4" gutterBottom component="div">
                Job Details for Job ID: {job_id || 'Loading ...'}
            </Typography>
          </Box>
          
          <Grid container spacing={3}>
            <Grid item xs={12} md={6}>
              <Card raised sx={cardStyle}>
                <CardHeader title="Dataset Information" />
                <CardContent sx={cardContentStyle}>
                  <Typography variant="subtitle1"><strong>Job Description:</strong></Typography>
                  <Typography variant="body1" gutterBottom>{description || 'Not available'}</Typography>
                  <Typography variant="subtitle1"><strong>Dataset:</strong></Typography>
                  <Typography variant="body1" gutterBottom>
                    { /* <Button onClick={downloadFile(datasetURL)}>
                      {getFileNameFromURL(datasetURL) || 'Not available'}
                    </Button> */ }
                    { process === "Integration" ?
                      ( datasetURL.map((inpput, index) => (
                          <a download onClick={() => { downloadFile(inpput) }} style={{ marginLeft: '10px', textAlign: 'center' }}>
                            {getFileNameFromURL(inpput) || 'Not available'}
                          </a>
                        ))
                      ) :
                      (<a download onClick={() => { downloadFile(datasetURL) }} style={{ marginLeft: '10px', textAlign: 'center' }}>
                        {getFileNameFromURL(datasetURL) || 'Not available'}
                      </a>)
                    }
                  </Typography>
                </CardContent>
              </Card>
            </Grid>

            <Grid item xs={12} md={6}>
              <Card raised sx={cardStyle}>
                <CardHeader title="Execution Details" />
                <CardContent sx={cardContentStyle}>
                  <Grid container spacing={2}> {/* Create a Grid container to layout details side by side */}
                  <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Process:</strong></Typography>
                      <Typography variant="body1">{process || 'Not available'}</Typography>
                    </Grid>
                    <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Method:</strong></Typography>
                      <Typography variant="body1">{method || 'Not available'}</Typography>
                    </Grid>
                    <Grid item xs={12}> 
                      <Typography variant="subtitle1" gutterBottom><strong>Status:</strong></Typography>
                      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                        <StatusChip status={taskStatus} />
                      </Box>
                    </Grid>
                </Grid>
                </CardContent>
              </Card>
            </Grid>


            {(status?.toLowerCase() !=="success" && status?.toLowerCase() !=="failure") && <Grid item xs={12}  >
              <Card raised sx={cardStyle}>
                <CardHeader title="Live Logs" />
                <CardContent sx={cardContentStyle}>
                  <LogComponent wsLogs={liveLogs} />
                </CardContent>
              </Card>
            </Grid>}

            {taskStatus?.toLowerCase() === "success" && taskOutput && (
              <Grid item xs={12}  >
                <Card raised sx={cardStyle}>
                  <CardHeader title="Task Results" />
                  <CardContent sx={cardContentStyle}>
                    {
                      taskOutput.map((output, index) => (
                        Object.keys(output).map((key) => (
                          <><Typography variant="subtitle1"><strong>{key}: </strong></Typography>
                            <Typography variant="body1" gutterBottom>
                              { <a download onClick={() => { downloadFile(output[key]); } } style={{ marginLeft: '10px', textAlign: 'center' }}>
                                   {getFileNameFromURL(output[key]) || 'Not available'}
                              </a> }
                            </Typography></>
                        ))
                      )
                    )}
                  </CardContent>
                </Card>
              </Grid>
            )}

            {/* {showErrorLog && (    */}
            {taskStatus?.toLowerCase() === "failure" && (
              <Grid item xs={12}>
                <Card raised sx={cardStyle}>
                  <CardHeader title="Error Feedback" />
                  <CardContent sx={cardContentStyle}>
                    <Typography variant="subtitle1"><strong>User Name:</strong> {uName}</Typography>
                    <Typography variant="subtitle1"><strong>User ID:</strong> {uIat}</Typography>
                    { /*<Typography variant="subtitle1"><strong>Task Result:</strong> {taskResult}</Typography>*/ }
                    <Typography variant="subtitle1"><strong>Job Status:</strong> {taskStatus}</Typography>
                    <Typography variant="subtitle1"><strong>Job ID:</strong> {job_id}</Typography>
                    <Typography variant="subtitle1"><strong>User Comments:</strong></Typography>
                    <TextField
                      fullWidth
                      multiline
                      rows={4}
                      variant="outlined"
                      value={userComment}
                      onChange={(e) => setUserComment(e.target.value)}
                      placeholder="Enter your comments here."
                    />
                    <Button
                      variant="contained"
                      color="primary"
                      onClick={handleSaveComment}
                      disabled={isSaving || isSent} // Disable button if saving
                      sx={{ mt: 2 }}
                    >
                      {isSaving ? 'Sending' : 'Send Feedback'}
                    </Button>
                    {commentSuccessMessage && (
                      <Typography variant="body1" color="success.main" sx={{ mt: 2 }}>
                        {commentSuccessMessage}
                      </Typography>
                    )}
                  </CardContent>
                </Card>
              </Grid>
            )}

          </Grid>
        </Container>

      {loading ? (
        <div className="spinner-container">
          <ScaleLoader color="#36d7b7" loading={loading} />
        </div>
      ) : (
          toolResultsFromMongo && (
          <div align="center"> 
          {
            toolResultsFromMongo.map((result, index) => (
              <React.Fragment key={index}>
                    {(result.umap_plot || result.umap_plot_3d) && (
                      <>
                        <h2>UMAP</h2>
                        <div style={{alignItems: 'center' }}>
                          <FormControl>
                            {/* <FormLabel id="demo-row-radio-buttons-group-label">Dimension</FormLabel> */}
                            <RadioGroup
                              row
                              aria-labelledby="demo-row-radio-buttons-group-label"
                              name="row-radio-buttons-group"
                              value={plotDimension}
                              onChange={(event) => setPlotDimension(event.target.value)}
                            >
                              <FormControlLabel value="2D" control={<Radio color="secondary" />} label="2D" />
                              <FormControlLabel value="3D" control={<Radio color="secondary" />} label="3D" />
                            </RadioGroup>
                          </FormControl>

                            <FormControl sx={{ m: 1, minWidth: 120 }} size="small">
                              <InputLabel id="plot-options-label">Color</InputLabel>
                              <Select
                                labelId="plot-options-label"
                                id="plot-options"
                                value={clusteringPlotType}
                                onChange={(event) => {
                                  const selectedPlotType = event.target.value;
                                  setClusteringPlotType(selectedPlotType);
                                  fetchPlotData(selectedPlotType, result.obs, result.umap, result.umap_3d); // Call the javascript function as soon as the selection changes
                                }}
                              >
                                {Object.keys(result.cell_metadata).map((key) => (
                                  <MenuItem key={key} value={key}>{key}</MenuItem>
                                ))}
                              </Select>
                            </FormControl>
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
                        <h2>Violin</h2>
                        <ReactPlotly plot_data={result.violin_plot} />
                      </>
                    )}
                    {result.scatter_plot && (
                      <>
                        <h2>Scatter</h2>
                        <ReactPlotly plot_data={result.scatter_plot} />
                      </>
                    )}
                    {result.highest_expr_genes_plot && (
                      <>
                        <h2>Highest expression Genes</h2>
                        <ReactPlotly plot_data={result.highest_expr_genes_plot} />
                      </>
                    )}
                  </React.Fragment>
          ))}
          </div>
      )
      )}
      </div>
      <div className="right-rail">
          <RightRail />
      </div>
  </div>
  );
}

export default TaskDetailsComponent;
