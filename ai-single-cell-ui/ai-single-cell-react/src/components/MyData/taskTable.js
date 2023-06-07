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
    }, [changesFound]);

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


    if (jwtToken)
        return (
            <TableWrapper>
                <h1>My Tasks</h1>
                <Table>
                    <TableHeader>
                        <TableRow>
                            <TableHeaderCell>Task ID</TableHeaderCell>
                            <TableHeaderCell>Workflow</TableHeaderCell>
                            <TableHeaderCell>Created</TableHeaderCell>
                            <TableHeaderCell>Finished</TableHeaderCell>
                            <TableHeaderCell>Dataset ID</TableHeaderCell>
                            <TableHeaderCell>Status</TableHeaderCell>
                            <TableHeaderCell>Result</TableHeaderCell>
                        </TableRow>
                    </TableHeader>
                    <tbody>
                        {tasks.map((task) => (
                            <TableRow key={task.task_id}>
                                <TableCell>
                                    <Typography variant="body2" style={{ whiteSpace: 'pre-wrap' }}>
                                        {task.task_id}
                                    </Typography>
                                </TableCell>
                                <TableCell>{task.workflow}</TableCell>
                                <TableCell>{new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(task.created_datetime))}</TableCell>
                                <TableCell>
                                    {task.finish_datetime ? (
                                        new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(task.finish_datetime))
                                    ) : (
                                        ''
                                    )}
                                </TableCell>
                                <TableCell>{task.dataset_id}</TableCell>
                                <TableCell style={{ textAlign: 'center' }}>
                                    {task.status === 'complete' ? (
                                        <CheckCircleIcon style={{ color: 'green' }} />
                                    ) : task.status === 'fail' ? (
                                        <CancelIcon style={{ color: 'red' }} />
                                    ) : (
                                        <HourglassEmptyIcon style={{ color: 'gray' }} />
                                    )}
                                </TableCell>
                                {task.status ? (
                                    <TableCell><a onClick={() => navigate('/resultfiles', { state: { taskId: task.task_id } })}>View</a></TableCell>
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