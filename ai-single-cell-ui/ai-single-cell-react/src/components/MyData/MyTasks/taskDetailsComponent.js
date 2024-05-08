import React, { useState } from 'react';
import { useLocation } from 'react-router-dom';
import useWebSocket from './useWebSocket'; // Custom hook for WebSocket
import { 
  Container, Typography, Chip, Box, CircularProgress, Paper, Grid, 
  Card, CardContent, Link, CardHeader 
} from '@mui/material';
import { green, red, yellow } from '@mui/material/colors';
import RightRail from '../../RightNavigation/rightRail';
import LogComponent from '../../common_components/liveLogs';
import { SERVER_URL } from '../../../constants/declarations';
import axios from 'axios';
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import { ScaleLoader } from 'react-spinners';
import FormControl from '@mui/material/FormControl';
import Radio from '@mui/material/Radio';
import RadioGroup from '@mui/material/RadioGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import ReactPlotly from '../../publishDatasets/components/reactPlotly';

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

function TaskDetailsComponent() {
  const location = useLocation();
  const { taskId, method, datasetURL, datasetTitle , tool} = location.state || {};   
  const [taskStatus, setTaskStatus] = useState(null); // Set to null initially
  const [liveLogs, setLiveLogs] = useState('');
  const [loading, setLoading] = useState(true);
  const [toolResultsFromMongo, setToolResultsFromMongo] = useState([]);

  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [ isError, setIsError ] = useState(false);
  const [plotDimension, setPlotDimension] = useState('2D');

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
        const response = await axios.post(`${SERVER_URL}/benchmarks/api/getPreProcessResults`, { processIds });
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
      if (data.task_status) {
        setTaskStatus(data.task_status);
        if(data.task_status === "SUCCESS" || data.task_status === "FAILURE"){
          if(data.task_status === "SUCCESS") {
            fetchProcessResults(data.task_result.process_ids);
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

  const handleLogMessage = (event) => {
    setLiveLogs(event.data);
    // Auto-scroll to the bottom of the logs
    const logsElement = document.getElementById("_live_logs");
    if (logsElement) {
      logsElement.scrollTop = logsElement.scrollHeight;
    }
  };

    // Use the WebSocket hook
    useWebSocket(taskId, handleStatusMessage, handleLogMessage);

  return (

    <div className="task-details-container eighty-twenty-grid">

      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}

      <div className="main-content">
        <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
          <Box display="flex" justifyContent="center">
            <Typography variant="h4" gutterBottom component="div">
                Task Details for Task ID: {taskId || 'Loading ...'}
            </Typography>
          </Box>
          
          <Grid container spacing={3}>
            <Grid item xs={12} md={6}>
              <Card raised sx={cardStyle}>
                <CardHeader title="Dataset Information" />
                <CardContent sx={cardContentStyle}>
                  <Typography variant="subtitle1"><strong>Dataset Title:</strong></Typography>
                  <Typography variant="body1" gutterBottom>{datasetTitle || 'Not available'}</Typography>
                  <Typography variant="subtitle1"><strong>Dataset URL:</strong></Typography>
                  <Typography variant="body1" gutterBottom>
                    <Link href={datasetURL} target="_blank" rel="noopener">
                      {datasetURL || 'Not available'}
                    </Link>
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
                      <Typography variant="subtitle1" gutterBottom><strong>Tool:</strong></Typography>
                      <Typography variant="body1">{tool || 'Not available'}</Typography>
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

            <LogComponent wsLogs = {liveLogs}/>
          </Grid>
        </Container>

      {loading ? (
        <div className="spinner-container">
          <ScaleLoader color="#36d7b7" loading={loading} />
        </div>
      ) : (
        (tool === "Quality Control" || tool === "Normalization") && (
          <div>
          {toolResultsFromMongo &&
            toolResultsFromMongo.map((result, index) => (
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
