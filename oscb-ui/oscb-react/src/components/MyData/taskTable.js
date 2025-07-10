import { faEye, faTrash } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import styled from 'styled-components';
import React, { useEffect, useState } from 'react';
import axios from 'axios';
import moment from 'moment';
import { getCookie } from '../../utils/utilFunctions';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';
import HourglassEmptyIcon from '@material-ui/icons/HourglassEmpty';
import { Typography } from '@material-ui/core';
import Intl from 'intl';
import 'intl/locale-data/jsonp/en-US';
import { useNavigate } from 'react-router-dom';
import { NODE_API_URL, WEB_SOCKET_URL, CELERY_BACKEND_API } from '../../constants/declarations'
import Pagination from '../publishDatasets/components/tablePaginationComponent';
import '../publishDatasets/publishDatasets.css';

const TableWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  margin-top: 50px;
`;

const Table = styled.table`
  border-collapse: collapse;
  width: 95%;
  /* max-width: 1000px; */
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

    // Destructure the pagination object for easier access to its properties
    // const { page, pageSize, totalCount } = pagination;

    // Calculate the starting and ending result numbers for the current page
    // const startResult = (page - 1) * pageSize + 1;
    // const endResult = Math.min(page * pageSize, totalCount); // Ensure not to exceed totalCount
    const [pagination, setPagination] = useState({});
    // Function to reset all state variables
    const resetState = () => {
        setPagination({});
    };

    const fetchTasks = async (currentPage) => {
        const response = await fetch(`${NODE_API_URL}/getTasks?authToken=${jwtToken}&page=${currentPage}`);
        const data = await response.json();
        data.results.sort((a, b) => a.created_on - b.created_on);
        setTasks(data.results);
        setPagination(data.pagination);

        // Create a list to store incomplete tasks
        const incompleteTasks = [];

        // Iterate over each task and check if its status is null
        data.results.forEach(task => {
            if (task.status === null) {
                incompleteTasks.push(task.job_id);
            }
        });

        if (incompleteTasks.length > 0) {
            let webSocketParam = incompleteTasks.join(',');
            const socket = new WebSocket(`${WEB_SOCKET_URL}/taskStatus/${webSocketParam}`);
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
                    else if (status === 'Failure') {
                        failedTasks.push(jobId);
                    }
                });
                if (finishedTasks.length + failedTasks.length > 0) {
                    // await updateTaskStatus(failedTasks, 'Failure');
                    // await updateTaskStatus(finishedTasks, 'Success');
                    // Close the WebSocket connection
                    socket.close(1000, 'See you again!');
                    setChangesFound(!changesFound);
                }
            };
        }
    };

    const onPageChange = (newPage) => {
        fetchTasks(newPage);
    };

    const handleDelete = (jobID) => {
        console.log("Delete job: ", jobID);
        const confirmDelete = window.confirm("Are you sure to delete this job?");
        if (!confirmDelete) {
            return; // If user clicks cancel, do nothing
        }
        axios.delete(`${NODE_API_URL}/deleteJob?jobID=${jobID}`)
            .then(response => {
                axios.post(`${CELERY_BACKEND_API}/task/revoke/${jobID}`).then(response => {
                    console.log('Job is deleted successfully');
                    fetchTasks(pagination.page);
                })
                    .catch(error => {
                        console.error('Error deleting job:', error);
                    });
            })
            .catch(error => {
                console.error('Error deleting job:', error);
            });
    };

    useEffect(() => {
        if (!jwtToken)
            navigate('/routing');       
        fetchTasks(1);
    }, [changesFound]);


    // const updateTaskStatus = async (jobIds, status) => {
    //     try {
    //         const jobIdString = jobIds.join(',');
    //         const response = await fetch(`${NODE_API_URL}/updateTaskStatus`, {
    //             method: 'PUT',
    //             headers: {
    //                 'Content-Type': 'application/json'
    //             },
    //             body: JSON.stringify({
    //                 jobIds: jobIdString,
    //                 status: status
    //             })
    //         });
    //         const data = await response.json();
    //         console.log(data);
    //         return data; // return the data from the function
    //     } catch (error) {
    //         console.error(error);
    //         throw error; // throw the error so that the caller can handle it
    //     }
    // };


    if (jwtToken)
        return (
            <><div className='table-results'>
                <TableWrapper>
                    <h1>My Jobs</h1>
                    <Table>
                        <TableHeader>
                            <TableRow>
                                <TableHeaderCell>Description</TableHeaderCell>
                                <TableHeaderCell>Job ID</TableHeaderCell>
                                <TableHeaderCell>Process</TableHeaderCell>
                                <TableHeaderCell>Method</TableHeaderCell>
                                <TableHeaderCell>Created on</TableHeaderCell>
                                <TableHeaderCell>Completed on</TableHeaderCell>
                                <TableHeaderCell>Status</TableHeaderCell>
                                <TableHeaderCell>Actions</TableHeaderCell>
                            </TableRow>
                        </TableHeader>
                        <tbody>
                            {tasks.map((task) => (
                                <TableRow key={task.job_id}>
                                    <TableCell>
                                        <Typography variant="body2" style={{ whiteSpace: 'pre-wrap' }}>
                                            {task.description}
                                        </Typography>
                                    </TableCell>
                                    <TableCell>
                                        <Typography variant="body2" style={{ whiteSpace: 'pre-wrap' }}>
                                            {task.job_id}
                                        </Typography>
                                    </TableCell>
                                    <TableCell>{task.process}</TableCell>
                                    <TableCell>{task.method}</TableCell>
                                    <TableCell>{new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(task.created_on).local()))}</TableCell>
                                    <TableCell>
                                        {new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(task.completed_on).local()))}
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
                                            <a onClick={() => handleDelete(task.job_id)}
                                                style={{ textDecoration: 'none', color: 'inherit' }}> <FontAwesomeIcon icon={faTrash} /></a>
                                            &nbsp;&nbsp;&nbsp;&nbsp;
                                            <a
                                                onClick={() => {
                                                if (task.process && task.category.toLowerCase() === 'workflow') {
                                                    navigate("/mydata/workflowTaskDetails", {
                                                        state: {
                                                        job_id: task.job_id,
                                                        methodMap: task.methodMap,
                                                        datasetURL: task.datasetURL,
                                                        description: task.description,
                                                        process: task.process,
                                                        output: task.output,
                                                        results: task.results,
                                                        status: task.status
                                                        }
                                                    });
                                                    } else {
                                                    navigate("/mydata/taskDetails", {
                                                        state: {
                                                        job_id: task.job_id,
                                                        method: task.method,
                                                        datasetURL: task.datasetURL,
                                                        description: task.description,
                                                        process: task.process,
                                                        output: task.output,
                                                        results: task.results,
                                                        status: task.status
                                                        }
                                                    });
                                                    }
                                                }}
                                                style={{ textDecoration: 'none', color: 'inherit' }}
                                                >
                                                <FontAwesomeIcon icon={faEye} />
                                                </a>
                                        </TableCell>
                                    ) : (
                                        <TableCell></TableCell>
                                    )}
                                </TableRow>
                            ))}
                        </tbody>
                    </Table>
                </TableWrapper>
                <p></p>
                <p></p>
                <p></p>
                <div className='table-pagination' align="center">
                    <Pagination
                        pagination={pagination}
                        onPageChange={onPageChange}
                    />
                </div>
                </div>
            </>
        );
};

export default TaskTable;