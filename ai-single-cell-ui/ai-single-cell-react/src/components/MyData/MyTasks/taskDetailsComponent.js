import React, { useState } from 'react';
import { useLocation } from 'react-router-dom';
import useWebSocket from './useWebSocket'; // Custom hook for WebSocket
import { 
  Container, Typography, Chip, Box, CircularProgress, Paper, Grid, 
  Card, CardContent, Link, CardHeader 
} from '@mui/material';
import { green, red, yellow } from '@mui/material/colors';
import { WEB_SOCKET_URL } from "../../../constants/declarations";
import RightRail from '../../RightNavigation/rightRail';

function StatusChip({ status }) {
  const getStatusColor = () => {
    switch (status?.toLowerCase()) { // Ensure status is defined
      case 'success': return green[500];
      case 'failed': return red[500];
      case 'processing': return yellow[700];
      default: return yellow[700];
    }
  };

  return status ? (
    <Chip label={status.toUpperCase()} style={{ backgroundColor: getStatusColor(), color: '#fff' }} />
  ) : (
    <Chip label="LOADING" style={{ backgroundColor: yellow[700], color: '#fff' }} />
  );
}

function TaskDetailsComponent() {
  const location = useLocation();
  const { taskId, method, datasetURL, datasetTitle , tool} = location.state || {};   
  const [taskStatus, setTaskStatus] = useState(null); // Set to null initially
  const [liveLogs, setLiveLogs] = useState('');

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
  

  const handleStatusMessage = (event) => {
    try {
      const data = JSON.parse(event.data);
      if (data[taskId]) {
        setTaskStatus(data[taskId]);
      }
    } catch (error) {
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

  // useWebSocket(`${WEB_SOCKET_URL}/taskStatus/${taskId}`, handleStatusMessage);
  // useWebSocket(`${WEB_SOCKET_URL}/log/${taskId}`, handleLogMessage);

   // Use the WebSocket hook
   useWebSocket(taskId, handleStatusMessage, handleLogMessage);

  return (

    <div className="task-details-container eighty-twenty-grid">
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

            <Grid item xs={12}>
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
                      dangerouslySetInnerHTML={createMarkup(liveLogs || 'Waiting for logs...')}
                    />
                  </Paper>
                </CardContent>
              </Card>
            </Grid>
          </Grid>
        </Container>
      </div>
      <div className="right-rail">
          <RightRail />
      </div>
  </div>
  );
}

export default TaskDetailsComponent;
