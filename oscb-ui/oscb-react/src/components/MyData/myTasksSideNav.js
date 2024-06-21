import React, { useState, useEffect } from 'react';
import { faArrowUpRightFromSquare, faAngleRight, faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { getCookie } from '../../utils/utilFunctions';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';
import HourglassEmptyIcon from '@material-ui/icons/HourglassEmpty';
import { useNavigate } from 'react-router-dom';
import TextWithEllipsis from '../RightNavigation/textWithEllipsis';
import {
    Accordion,
    AccordionSummary,
    AccordionDetails
  } from '@material-ui/core';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';

const MyTasksSideNav = () => {
    const [expanded, setExpanded] = useState(false);
    const [tasks, setTasks] = useState([]);
    const [changesFound, setChangesFound] = useState(false);
    let jwtToken = getCookie('jwtToken');
    const navigate = useNavigate();
    const NODE_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3001`;
    const WEB_SOCKET_URL = `ws://${process.env.REACT_APP_HOST_URL}:5000/taskStatus`;

    const timestampScheme = {
        year: 'numeric',
        month: 'short',
        day: 'numeric',
        hour: 'numeric',
        minute: 'numeric',
        hour12: true
    };

    useEffect(() => {
        if (jwtToken && expanded) {
            const fetchTasks = async () => {
                const response = await fetch(`${NODE_API_URL}/getTasks?authToken=${jwtToken}`);
                const data = await response.json();
                data.sort((a, b) => b.created_on - a.created_on);
                setTasks(data);

                // Create a list to store incomplete tasks
                const incompleteTasks = [];

                // Iterate over each task and check if its status is null
                data.forEach(task => {
                    if (task.status === null) {
                        incompleteTasks.push(task.job_id);
                    }
                });

                if (incompleteTasks.length > 0) {
                    let webSocketParam = incompleteTasks.join(',');
                    const socket = new WebSocket(`${WEB_SOCKET_URL}/${webSocketParam}`);
                    socket.onopen = () => {
                        console.log('Socket connected');
                    };
                    socket.onclose = () => {
                        console.log('Socket disconnected');
                    };

                    let finishedTasks = [];
                    let failedTasks = [];
                    socket.onmessage = async (event) => {
                        const data = JSON.parse(event.data);
                        Object.keys(data).forEach(jobId => {
                            const status = data[jobId];
                            if (status === 'Success') {
                                finishedTasks.push(jobId);
                            }
                            else if (status === 'Failed') {
                                failedTasks.push(jobId);
                            }
                        });
                        if (finishedTasks.length + failedTasks.length > 0) {
                            await updateTaskStatus(failedTasks, 'Failed');
                            await updateTaskStatus(finishedTasks, 'Success');
                            // Close the WebSocket connection
                            socket.close(1000, 'See you again!');
                            setChangesFound(!changesFound);
                        }
                    };
                }
            };
            fetchTasks();
        }
    }, [changesFound, expanded]);

    const updateTaskStatus = async (jobIds, status) => {
        try {
            const jobIdString = jobIds.join(',');
            const response = await fetch(`${NODE_API_URL}/updateTaskStatus`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    jobIds: jobIdString,
                    status: status
                })
            });
            const data = await response.json();
            console.log(data);
            return data; // return the data from the function
        } catch (error) {
            console.error(error);
            throw error; // throw the error so that the caller can handle it
        }
    };

    const toggleExpand = () => {
        setExpanded(!expanded);
    };

    if (jwtToken)
        return (
            <div className="expandable">
                <div className="header" onClick={toggleExpand}>
                    <span className={`arrow ${expanded ? 'expanded' : ''}`}>
                        <FontAwesomeIcon icon={faAngleRight} />
                    </span>
                    <span className="title">My Jobs</span>
                    &nbsp;
                    <FontAwesomeIcon
                        icon={faArrowUpRightFromSquare}
                        className="hoverable-icon"
                        onClick={() => { navigate('/myTasks') }}
                        style={{ textAlign: 'right' }}
                    />
                </div>
                {expanded && (
                    <div className="content">
                        <div style={{ maxHeight: '360px', overflow: 'auto' }}>

                        {tasks.length === 0 ? (
                            <div>
                                <div>
                                    <div role="alert" aria-live="polite" aria-atomic="true" className="alert m-2 alert-info">
                                        <h4 className="mb-1">
                                            <FontAwesomeIcon icon={faInfoCircle} />
                                            <span>Your task list is empty.</span>
                                        </h4> 
                            </div>
                            </div>
                            </div>
                        ) : (
                            <ul>
                            {tasks.map((task, index) => (
                                <div key={tasks.job_id}>
                                    <Accordion key={tasks.job_id}>
                                      <AccordionSummary expandIcon={<ExpandMoreIcon />} aria-controls="panel-content" id="panel-header">
                                        <div className="panel-summary">
                                            <div className='task-summary'>
                                                <div className='display-flex'>
                                                    {task.status === 'Success' ? (
                                                        <CheckCircleIcon style={{ color: 'green' }} />
                                                    ) : task.status === 'Failed' ? (
                                                        <CancelIcon style={{ color: 'red' }} />
                                                    ) : (
                                                        <HourglassEmptyIcon style={{ color: 'gray' }} />
                                                    )}
                                                    <p><TextWithEllipsis text={task.task_title} maxLength={23} /></p>
                                                </div>
                                                <span className='time-stamp-display'>- {new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(task.created_on))}</span>
                                            </div>
                                        {/* <li style={{
                                    backgroundColor: 'transparent', // Set initial background color
                                    transition: 'background-color 0.3s', // Add transition effect
                                    cursor: 'pointer' // Show pointer cursor on hover
                                }}
                                    onMouseEnter={(e) => { e.target.style.backgroundColor = '#f2f2f2' }} // Change background color on hover
                                    onMouseLeave={(e) => { e.target.style.backgroundColor = 'transparent' }} // Revert back to initial background color on mouse leave 
                                    key={index}> */}
                                      {/* <a
                                        href={`/resultfiles?jobId=${task.job_id}&results_path=${task.results_path}`}
                                        style={{ textDecoration: 'none', color: 'inherit' }}
                                    > 
                                        {task.status === 'Success' ? (
                                        <CheckCircleIcon style={{ color: 'green' }} />
                                    ) : task.status === 'Failed' ? (
                                        <CancelIcon style={{ color: 'red' }} />
                                    ) : (
                                        <HourglassEmptyIcon style={{ color: 'gray' }} />
                                    )}
                                    &nbsp;{task.job_id}
                                    </a> */}
                                {/* </li> */}
                                        </div>
                                      </AccordionSummary>
                                      <AccordionDetails>
                                      <a
                                        href={`/resultfiles?jobId=${task.job_id}&results_path=${task.results_path}&task_title=${task.task_title}`}
                                        style={{ textDecoration: 'none', color: 'inherit' }}
                                     > 
                                        <span className='font-size'><b>Task Id</b> - {task.job_id}</span>
                                      </a>
                                      </AccordionDetails>
                                    </Accordion>
                                </div>
                            ))
                            }
                        </ul>
                        )}
                        </div>
                    </div>
                )}
            </div>
        );
};

export default MyTasksSideNav;