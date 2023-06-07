import React, { useState, useEffect } from 'react';
import { faArrowUpRightFromSquare, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { getCookie } from '../../utils/utilFunctions';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';
import HourglassEmptyIcon from '@material-ui/icons/HourglassEmpty';
import { useNavigate } from 'react-router-dom';

const MyTasksSideNav = () => {
    const [expanded, setExpanded] = useState(false);
    const [tasks, setTasks] = useState([]);
    const [changesFound, setChangesFound] = useState(false);
    let jwtToken = getCookie('jwtToken');
    const navigate = useNavigate();
    const NODE_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3001`;
    const WEB_SOCKET_URL = `ws://${process.env.REACT_APP_HOST_URL}:80/taskStatus`;

    useEffect(() => {
        if (jwtToken && expanded) {
            const fetchTasks = async () => {
                const response = await fetch(`${NODE_API_URL}/getTasks?authToken=${jwtToken}`);
                const data = await response.json();
                data.sort((a, b) => a.created_datetime - b.created_datetime);
                setTasks(data);

                // Create a list to store incomplete tasks
                const incompleteTasks = [];

                // Iterate over each task and check if its status is null
                data.forEach(task => {
                    if (task.status === null) {
                        incompleteTasks.push(task.task_id);
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
                        Object.keys(data).forEach(taskId => {
                            const status = data[taskId];
                            if (status === 'complete') {
                                finishedTasks.push(taskId);
                            }
                            else if (status === 'fail') {
                                failedTasks.push(taskId);
                            }
                        });
                        if (finishedTasks.length + failedTasks.length > 0) {
                            await updateTaskStatus(failedTasks, 'fail');
                            await updateTaskStatus(finishedTasks, 'complete');
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

    const updateTaskStatus = async (taskIds, status) => {
        try {
            const taskIdString = taskIds.join(',');
            const response = await fetch(`${NODE_API_URL}/updateTaskStatus`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    taskIds: taskIdString,
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
                    <span className="title">My Tasks</span>
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
                        <div style={{ maxHeight: '360px', overflow: 'auto' }}><ul>
                            {tasks.map((task, index) => (
                                <li style={{
                                    backgroundColor: 'transparent', // Set initial background color
                                    transition: 'background-color 0.3s', // Add transition effect
                                    cursor: 'pointer' // Show pointer cursor on hover
                                }}
                                    onMouseEnter={(e) => { e.target.style.backgroundColor = '#f2f2f2' }} // Change background color on hover
                                    onMouseLeave={(e) => { e.target.style.backgroundColor = 'transparent' }} // Revert back to initial background color on mouse leave 
                                    onClick={() => navigate('/resultfiles', { state: { taskId: task.task_id } })} key={index}>{task.status === 'complete' ? (
                                        <CheckCircleIcon style={{ color: 'green' }} />
                                    ) : task.status === 'fail' ? (
                                        <CancelIcon style={{ color: 'red' }} />
                                    ) : (
                                        <HourglassEmptyIcon style={{ color: 'gray' }} />
                                    )}
                                    &nbsp;{task.task_id}
                                </li>
                            ))}
                        </ul></div>
                    </div>
                )}
            </div>
        );
};

export default MyTasksSideNav;