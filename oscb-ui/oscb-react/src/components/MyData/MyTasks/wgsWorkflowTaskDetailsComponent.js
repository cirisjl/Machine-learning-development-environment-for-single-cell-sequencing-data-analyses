import React, { useState, useEffect, useCallback } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import useWebSocket from './useWebSocket';
import {
  Container, Typography, Chip, Box, Grid, TextField, Button,
  Card, CardContent, CardHeader, Select, MenuItem, InputLabel,
  FormControl, Radio, RadioGroup, FormControlLabel
} from '@mui/material';
import { green, red, yellow } from '@mui/material/colors';
import RightRail from '../../RightNavigation/rightRail';
import LogComponent from '../../common_components/liveLogs';
import axios from 'axios';
// import { Octokit } from "@octokit/rest";
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import { ScaleLoader } from 'react-spinners';
import ReactPlotly from '../../publishDatasets/components/reactPlotly';
import { isUserAuth, getCookie, plotUmapObs, gunzipDict, ShowVitessce } from '../../../utils/utilFunctions';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { Descriptions } from 'antd';
import { styled } from '@mui/material/styles';
import ArrowForwardIosSharpIcon from '@mui/icons-material/ArrowForwardIosSharp';
import MuiAccordion from '@mui/material/Accordion';
import MuiAccordionSummary from '@mui/material/AccordionSummary';
import MuiAccordionDetails from '@mui/material/AccordionDetails';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import Chatbot from "../../RightNavigation/Chatbot";

// GitImports & Constants
import { CELERY_BACKEND_API, NODE_API_URL } from '../../../constants/declarations';
import TaskImageGallery from './taskImageGallery';
import AnnotationTable from './annotationPanel';
import OutlierTable from './outlierPanel';
import { faFolderOpen, faInfoCircle } from "@fortawesome/free-solid-svg-icons";
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemText from '@mui/material/ListItemText';
import IconButton from '@mui/material/IconButton';
import DeleteIcon from '@mui/icons-material/Delete';
import close_icon from '../../../assets/close_icon_u86.svg';
import close_icon_hover from '../../../assets/close_icon_u86_mouseOver.svg';
import FilePreviewModal from "../filePreviewModal";
import FileManagerModal from "../fileManagerModal";
import Toggle from 'react-toggle';
import 'react-toggle/style.css';

// const octokit = new Octokit({ auth: process.env.REACT_APP_TOKEN });
const jwtToken = getCookie('jwtToken');

// --- Styled Components ---
const Accordion = styled((props) => (
  <MuiAccordion disableGutters elevation={0} square {...props} />
))(({ theme }) => ({
  border: `1px solid ${theme.palette.divider}`,
  '&:not(:last-child)': { borderBottom: 0 },
  '&::before': { display: 'none' },
}));

const panelStyle = {
  maxHeight: '150px',       // Restricts the height
  overflowY: 'auto',        // Adds a scrollbar ONLY if content exceeds maxHeight
  border: '1px solid #ccc',
  padding: '16px',
  borderRadius: '8px',
  backgroundColor: '#f9f9f9'
};

const AccordionSummary = styled((props) => (
  <MuiAccordionSummary expandIcon={<ArrowForwardIosSharpIcon sx={{ fontSize: '0.9rem' }} />} {...props} />
))(({ theme }) => ({
  backgroundColor: 'rgba(0, 0, 0, .03)',
  flexDirection: 'row-reverse',
  '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': { transform: 'rotate(90deg)' },
  '& .MuiAccordionSummary-content': { marginLeft: theme.spacing(1) },
  ...theme.applyStyles('dark', { backgroundColor: 'rgba(255, 255, 255, .05)' }),
}));

const AccordionDetails = styled(MuiAccordionDetails)(({ theme }) => ({
  padding: theme.spacing(2),
  borderTop: '1px solid rgba(0, 0, 0, .125)',
}));

const cardStyle = { height: '100%', display: 'flex', flexDirection: 'column' };
const cardContentStyle = { flexGrow: 1, overflow: 'auto' };

// --- Helper Functions & Small Components ---
const StatusChip = ({ status }) => {
  const getStatusColor = () => {
    switch (status?.toLowerCase()) {
      case 'success': return green[500];
      case 'failure': return red[500];
      case 'started': return yellow[700];
      default: return yellow[700];
    }
  };
  return (
    <Chip
      label={(status || "In Progress").toUpperCase()}
      style={{ backgroundColor: getStatusColor(), color: '#fff' }}
    />
  );
};

const getFileNameFromURL = (fileUrl) => fileUrl ? fileUrl.substring(fileUrl.lastIndexOf('/') + 1) : '';

const downloadFile = (fileUrl) => {
  if (!fileUrl) return;
  const apiUrl = `${NODE_API_URL}/download`;
  const pwd = "jobResults";

  fetch(`${apiUrl}?fileUrl=${fileUrl}&authToken=${jwtToken}&pwd=${pwd}`)
    .then(res => res.blob())
    .then(blob => {
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = getFileNameFromURL(fileUrl);
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    })
    .catch(err => console.error('Error downloading file:', err));
};

// Reusable Download Link Component to prevent repetitive a-tags
const DownloadLink = ({ url, label }) => (
  <button type="button"
    onClick={(e) => { e.preventDefault(); downloadFile(url); }}
    style={{ background: 'none', border: 'none', padding: 0, font: 'inherit', cursor: 'pointer', textDecoration: 'underline', color: 'blue', marginLeft: '10px' }}
  >
    {label || getFileNameFromURL(url) || 'Not available'}
  </button>
);

// --- Main Component ---
export default function WgsWorkflowTaskDetailsComponent() {
  const location = useLocation();
  const navigate = useNavigate();
  const { job_id, methodMap, datasetURL, description, process, output, results, status } = location.state || {};

  const [taskStatus, setTaskStatus] = useState(null);
  const [taskOutput, setTaskOutput] = useState(null);
  const [liveLogs, setLiveLogs] = useState('');
  const [loading, setLoading] = useState(true);

  const [uName, setUName] = useState(null);
  const [uIat, setUIat] = useState(null);
  const [taskResult, setTaskResult] = useState("");
  const [message, setMessage] = useState('');
  const [hasMessage, setHasMessage] = useState(false);
  const [isError, setIsError] = useState(false);

  const [userComment, setUserComment] = useState('');
  const [isSaving, setIsSaving] = useState(false);
  const [isSent, setIsSent] = useState(false);
  const [commentSuccessMessage, setCommentSuccessMessage] = useState('');
  const [showErrorLog, setShowErrorLog] = useState(true);

  const [presetQuestions, setPresetQuestions] = useState(null);

  const isTerminalState = useCallback((s) => ["success", "failure"].includes((s || '').toLowerCase()), []);

  const [isFileManagerOpen, setIsFileManagerOpen] = useState(false);
  const [fileManagerTarget, setFileManagerTarget] = useState('dataset');
  const [pwd, setPwd] = useState('/');
  let jwtToken = getCookie('jwtToken');
  const [enabledCheckboxes, setEnabledCheckboxes] = useState([]);
  const [tempFileList, setTempFileList] = useState([]);
  const [hoveredErrPopup, setHoveredErrPopup] = useState(false);
  const [previewBoxOpen, setPreviewBoxOpen] = useState(false);
  const [fileToPreview, setFileToPreview] = useState(null);
  const [fileNames, setFileNames] = useState([]);
  const [dirNames, setDirNames] = useState([]);
  const [publicDatasetFlag, setPublicDatasetFlag] = useState(false);
  const [isAdminuser, setIsAdminuser] = useState(false);

  const handleDelete = () => {
    if (!window.confirm("Are you sure to delete this job?")) return;
    axios.post(`${CELERY_BACKEND_API}/job/revoke/${job_id}`)
      .then(() => axios.delete(`${NODE_API_URL}/deleteJob?jobID=${job_id}`,
        {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${jwtToken}`,
          },
        }))
      .then(() => navigate('/myJobs'))
      .catch(error => console.error('Error deleting job:', error));
  };

  const handleStatusMessage = useCallback(({ data }) => {
    const currentStatus = status || data?.task_status || taskStatus;
    if (!currentStatus) return;
    setTaskStatus(currentStatus);

    if (isTerminalState(currentStatus.toLowerCase())) {
      setLoading(false);
      if (currentStatus.toLowerCase() === "success") {
        if (results?.process_ids) {
          setTaskResult(results);
        } else if (data.task_result?.process_ids) {
          const resultData = data.task_result;
          if (resultData.output) {
            setTaskOutput(resultData.output);
            setTaskResult(resultData);
          }
        }
      }
    }
  }, [isTerminalState, status, taskStatus, results]); // added results to dependencies

  useEffect(() => {
    async function fetchFiles() {
      if (status?.toLowerCase() === "success" && output) {
        setTempFileList([]);
        setEnabledCheckboxes([]);
        setPwd(output.replace("/usr/src/app/storage/" + uName, ""));
        fetchDirContents();
        setIsFileManagerOpen(true);
        setTaskOutput(output);
        setTaskResult(results);
      } else {
        try {
          const res = await fetch(`${CELERY_BACKEND_API}/task/${job_id}`);
          const data = await res.json();
          setTaskResult(data.task_result);
        } catch (error) { console.error('Error fetching task status:', error); }
      }
    }
    const fetchUserInfo = async () => {
      if (!jwtToken) return;
      try {
        const res = await fetch(NODE_API_URL + "/protected", { headers: { Authorization: `Bearer ${jwtToken}` } });
        const data = await res.json();
        if (data?.authData) { setUName(data.authData.username); setUIat(data.authData.iat); }
      } catch (e) { console.error(e); }
    };
    fetchUserInfo();
    fetchFiles();
  }, [job_id, status, output, results]);

  const handleLogMessage = useCallback((event) => {
    setLiveLogs((prev) => prev + event.data);
    const logsElement = document.getElementById("_live_logs");
    if (logsElement) logsElement.scrollTop = logsElement.scrollHeight;
  }, []);

  const handleSaveComment = async () => {
    setIsSaving(true);
    setIsSent(false);
    try {
      const res = await axios.post(`${NODE_API_URL}/errorlogdata`, {
        name: uName, 
        id: uIat, 
        taskResult, 
        taskStatus, 
        job_id, 
        userComments: userComment
      });

      // Issue creation
      const data = await axios.post(`${NODE_API_URL}/issues`,
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            title: `Issue for Job ID: ${job_id}`,
            body: `User: ${uName}\nJob ID: ${job_id}\nStatus: ${taskStatus}\nComment: ${userComment}`
          }),
        }
      );

      console.log(`Issue created: ${data.html_url}`);
      if (res.status === 200 && data.status === 200) {
        setCommentSuccessMessage('Feedback sent successfully.');
        setIsSent(true);
        setShowErrorLog(false);
      }
    } catch (error) {
      setCommentSuccessMessage('Failed to send feedback.');
    } finally {
      setIsSaving(false);
    }
  };

  useWebSocket(job_id, handleStatusMessage, handleLogMessage, setLoading);

// Custom styled components
  const ScrollableListContainer = styled('div')(({ theme }) => ({
      maxHeight: '600px', // Fixed height of the container
      overflowY: 'auto', // Enable vertical scrolling
      border: `1px solid ${theme.palette.divider}`, // Add border to distinguish the container
      borderRadius: theme.shape.borderRadius, // Use theme's border radius
      marginTop: theme.spacing(2),
  }));

  const CustomListItem = styled(ListItem)(({ theme }) => ({
      '&:hover': {
          backgroundColor: theme.palette.action.hover,
      },
      cursor: 'pointer', // Change cursor on hover to indicate an item is clickable
  }));

  const toggleModal = async() => {
    await setIsFileManagerOpen(!isFileManagerOpen);

    if (!isFileManagerOpen) {
      setEnabledCheckboxes([]);
      setTempFileList([]);
      setPwd('/');
      fetchDirContents();
      }
    }

  // Open the shared file manager for one of the two selectors.
  const openFileManager = (target) => {
    setFileManagerTarget(target);
    // Start a fresh selection session. Existing selections remain unchanged
    // unless the user confirms a new list with the Select Files button.
    setTempFileList([]);
    setEnabledCheckboxes([]);
    setPwd('/');
    fetchDirContents();
    setIsFileManagerOpen(true);
  };

  // Adapter passed to FileManagerModal. It routes the selected files to the
  // selector that opened the modal.
  const setActiveSelectedFiles = (files) => {
      const nextFiles = Array.isArray(files) ? files : [];

      // FileManagerModal sends [] when its Close button is used. Preserve the
      // existing selector value in that case; files can still be removed from
      // the list with the delete buttons in this workflow.
      if (nextFiles.length === 0) {
          return;
      }

      if (fileManagerTarget === 'sra') {
          setSelectedSraFiles(nextFiles);
      } else {
          setSelectedFiles(nextFiles);
    }
  };

  const fetchDirContents = async (subdir) => {
      enabledCheckboxes.forEach((checkboxId) => {
          try {
              const checkbox = document.getElementById(checkboxId);
              checkbox.checked = false;
          }
          catch (e) {
              console.log('Exception: ' + e);
          }
      });

      let newDir = ''
      jwtToken = getCookie('jwtToken');

      if (subdir === "publicDatasets") {
          newDir = subdir;
          setPwd(newDir)
      } else if (pwd === "publicDatasets" && subdir === "..") {
          newDir = ""
          setPwd(newDir)
      }
      else if (subdir === '..') {
          let slashIndex = pwd.lastIndexOf('/');
          newDir = slashIndex !== -1 ? pwd.substring(0, slashIndex) : pwd;
          setPwd(newDir)
      }
      else if (subdir) {
          newDir = pwd + '/' + subdir;
          setPwd(pwd + '/' + subdir);
      }
      else {
          newDir = pwd
      }

      fetch(`${NODE_API_URL}/getDirContents?dirPath=${newDir}&authToken=${jwtToken}`)
        .then(response => {
            if (response.status === 403) {
                throw new Error('Please log in first');
            }
            return response.json();
        })
        .catch(error => {
            setIsFileManagerOpen(false);
            if (error.message === 'Please log in first') {
                navigate('/routing');
            } else {
                console.error(error);
            }

        })
        .then(data => {
            if (!data) return;
            setDirNames(data.Directories || []);
            setFileNames(data.Files || []);
        });
  }

  return (
    <div className="task-details-container eighty-twenty-grid">
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage={setMessage} isError={isError} />}

      <div className="main-content">
        <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
          <Box display="flex" justifyContent="center">
            <Typography variant="h4" gutterBottom>Job Details for Job ID: {job_id || 'Loading ...'}</Typography>
          </Box>

          <Grid container spacing={3}>
            {/* Dataset Information Card */}
            <Grid item xs={12} md={6}>
              <Card raised sx={cardStyle}>
                <CardHeader title="Dataset Information" />
                <CardContent sx={cardContentStyle}>
                  <Typography variant="subtitle1"><strong>Job Description:</strong></Typography>
                  <Typography variant="body1" gutterBottom>{description || 'Not available'}</Typography>
                  <Typography variant="subtitle1"><strong>Dataset:</strong></Typography>
                  <Typography variant="body1" gutterBottom>
                    {Array.isArray(datasetURL)
                      ? datasetURL.map((input, idx) => <DownloadLink key={idx} url={input} />)
                      : <DownloadLink url={datasetURL} />
                    }
                  </Typography>
                </CardContent>
              </Card>
            </Grid>

            {/* Execution Details Card */}
            <Grid item xs={12} md={6}>
              <Card raised sx={cardStyle}>
                <CardHeader title="Execution Details" />
                <CardContent sx={cardContentStyle}>
                  <Grid container spacing={2}>
                    <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Workflow:</strong></Typography>
                      <Typography variant="body1">{process || 'Not available'}</Typography>
                    </Grid>
                    <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Method:</strong></Typography>
                      {methodMap && Object.keys(methodMap).length > 0 ? (
                        Object.entries(methodMap).map(([key, value]) => <Typography variant="body1" key={key}>{key}: {value}</Typography>)
                      ) : <Typography variant="body2" color="textSecondary">Not available</Typography>}
                    </Grid>
                    <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Status:</strong></Typography>
                      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}><StatusChip status={taskStatus} /></Box>
                    </Grid>
                    <Grid item xs={6}>
                      <Typography variant="subtitle1" gutterBottom><strong>Action:</strong></Typography>
                      <Button variant="contained" color="error" onClick={handleDelete}>
                        {(status && ["success", "failure"].includes(status.toLowerCase())) ? 'Delete' : 'Terminate'}
                      </Button>
                    </Grid>
                  </Grid>
                </CardContent>
              </Card>
            </Grid>

            {/* Live Logs Card */}
            {!["success", "failure"].includes(status?.toLowerCase()) && (
              <Grid item xs={12}>
                <Card raised sx={cardStyle}>
                  <CardHeader title="Live Logs" />
                  <CardContent sx={cardContentStyle}><LogComponent wsLogs={liveLogs} /></CardContent>
                </Card>
              </Grid>
            )}

            {/* Task Results Card */}
            {taskStatus?.toLowerCase() === "success" && taskOutput && (
              <Grid item xs={12}>
                <Card raised sx={cardStyle}>
                  <CardHeader title="Task Results" />
                  <CardContent sx={cardContentStyle}>
                    {previewBoxOpen && (
                      <div className="file-preview-layer">
                        <FilePreviewModal
                          selectedFile={fileToPreview}
                          setPreviewBoxOpen={setPreviewBoxOpen}
                          jwtToken={jwtToken}
                          forResultFile={false}
                        />
                      </div>
                    )}
                    {isFileManagerOpen && (
                      <FileManagerModal
                        setFileToPreview={setFileToPreview}
                        tempFileList={tempFileList}
                        setEnabledCheckboxes={setEnabledCheckboxes}
                        fileNames={fileNames}
                        dirNames={dirNames}
                        jwtToken={jwtToken}
                        fetchDirContents={fetchDirContents}
                        pwd={pwd}
                        setPwd={setPwd}
                        setPreviewBoxOpen={setPreviewBoxOpen}
                        selectedFiles={selectedFiles}
                        setSelectedFiles={setActiveSelectedFiles}
                        setFormErrors={setFormErrors}
                        setTempFileList={setTempFileList}
                        enabledCheckboxes={enabledCheckboxes}
                        toggleModal={toggleModal}
                        isAdminuser={isAdminuser}
                        publicDatasetFlag={publicDatasetFlag}
                      />
                    )}
                    {/* {taskOutput.map((outputObj, index) => (
                      Object.keys(outputObj).map((key) => (
                        <div key={`${index}-${key}`}>
                          <Typography variant="subtitle1"><strong>{key}: </strong></Typography>
                          <Typography variant="body1" gutterBottom><DownloadLink url={outputObj[key]} /></Typography>
                        </div>
                      ))
                    ))} */}
                    {taskResult?.wf_results?.figures && (
                      <div>
                        <Typography variant="subtitle1"><strong>Figures: </strong></Typography>
                        <TaskImageGallery figures={taskResult.wf_results.figures} />
                      </div>
                    )}
                  </CardContent>
                </Card>
              </Grid>
            )}

            {/* Error Feedback Card */}
            {taskStatus?.toLowerCase() === "failure" && showErrorLog && (
              <Grid item xs={12}>
                <Card raised sx={cardStyle}>
                  <CardHeader title="Error Feedback" />
                  <CardContent sx={cardContentStyle}>
                    <Typography variant="subtitle1"><strong>User Name:</strong> {uName}</Typography>
                    <Typography variant="subtitle1"><strong>User ID:</strong> {uIat}</Typography>
                    <Typography variant="subtitle1"><strong>Job Status:</strong> {taskStatus}</Typography>
                    <Typography variant="subtitle1"><strong>Job ID:</strong> {job_id}</Typography>
                    <Typography variant="subtitle1" sx={{ mt: 2 }}><strong>User Comments:</strong></Typography>
                    <TextField fullWidth multiline rows={4} variant="outlined" value={userComment} onChange={(e) => setUserComment(e.target.value)} placeholder="Enter your comments here." />
                    <Button variant="contained" color="primary" onClick={handleSaveComment} disabled={isSaving || isSent} sx={{ mt: 2 }}>
                      {isSaving ? 'Sending...' : 'Send Feedback'}
                    </Button>
                    {commentSuccessMessage && <Typography variant="body1" color="success.main" sx={{ mt: 2 }}>{commentSuccessMessage}</Typography>}
                  </CardContent>
                </Card>
              </Grid>
            )}
          </Grid>
        </Container>
        
      </div>

      <div className="right-rail">
        <RightRail />
        <Chatbot presetQuestions={presetQuestions} />
      </div>
    </div>
  );
}