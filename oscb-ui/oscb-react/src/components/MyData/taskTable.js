import styled from 'styled-components';
import React, { useEffect, useState } from 'react';
import { getCookie } from '../../utils/utilFunctions';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';
import HourglassEmptyIcon from '@material-ui/icons/HourglassEmpty';
import { Typography } from '@material-ui/core';
import Intl from 'intl';
import 'intl/locale-data/jsonp/en-US';
import { useNavigate } from 'react-router-dom';

const TableWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  margin-top: 50px;
`;

const Table = styled.table`
  border-collapse: collapse;
  width: 95%;
  max-width: 1000px;
  min-width: 500px;
`;

const TableHeader = styled.thead`
  background-color: #f1f1f1;
  font-weight: bold;
  text-align: left;
`;

const TableHeaderCell = styled.th`
  padding: 15px;
`;

const TableRow = styled.tr`
  border-bottom: 1px solid #ddd;
`;

const TableCell = styled.td`
  padding: 15px;
`;


const NODE_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3001`;
const WEB_SOCKET_URL = `ws://${process.env.REACT_APP_HOST_URL}:80/taskStatus`;

const TaskTable = () => {
    const [tasks, setTasks] = useState([]);
    const [changesFound, setChangesFound] = useState(false);
    let jwtToken = getCookie('jwtToken');
    const navigate = useNavigate();
    const timestampScheme = {
        year: 'numeric',
        month: 'short',
        day: 'numeric',
        hour: 'numeric',
        minute: 'numeric',
        hour12: true
    };

    useEffect(() => {
        if (!jwtToken)
            navigate('/routing');
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
    }, [changesFound]);

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


    if (jwtToken)
        return (
            <TableWrapper>
                <h1>My Tasks</h1>
                <Table>
                    <TableHeader>
                        <TableRow>
                            <TableHeaderCell>Task Title</TableHeaderCell>
                            <TableHeaderCell>Task ID</TableHeaderCell>
                            <TableHeaderCell>Tool</TableHeaderCell>
                            <TableHeaderCell>Created</TableHeaderCell>
                            <TableHeaderCell>Finished</TableHeaderCell>
                            <TableHeaderCell>Status</TableHeaderCell>
                            <TableHeaderCell>Result</TableHeaderCell>
                        </TableRow>
                    </TableHeader>
                    <tbody>
                        {tasks.map((task) => (
                            <TableRow key={task.job_id}>
                                <TableCell>{task.task_title}</TableCell>
                                <TableCell>
                                    <Typography variant="body2" style={{ whiteSpace: 'pre-wrap' }}>
                                        {task.job_id}
                                    </Typography>
                                </TableCell>
                                <TableCell>{task.tool}</TableCell>
                                <TableCell>{new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(task.created_datetime))}</TableCell>
                                <TableCell>
                                    {task.finish_datetime ? (
                                        new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(task.finish_datetime))
                                    ) : (
                                        ''
                                    )}
                                </TableCell>
                                <TableCell style={{ textAlign: 'center' }}>
                                    {task.status === 'Success' ? (
                                        <CheckCircleIcon style={{ color: 'green' }} />
                                    ) : task.status === 'Failed' ? (
                                        <CancelIcon style={{ color: 'red' }} />
                                    ) : (
                                        <HourglassEmptyIcon style={{ color: 'gray' }} />
                                    )}
                                </TableCell>
                                {task.status ? (
                                <TableCell>                                    
                                    <a href={`/resultfiles?jobId=${task.job_id}&results_path=${task.results_path}`}
                                    style={{ textDecoration: 'none', color: 'inherit' }}> View</a></TableCell>
                                     
                                ) : (
                                    <TableCell></TableCell>
                                )}
                            </TableRow>
                        ))}
                    </tbody>
                </Table>
            </TableWrapper>
        );
};

export default TaskTable;