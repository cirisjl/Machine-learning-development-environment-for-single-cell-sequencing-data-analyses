import { faEye, faTrash } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React, { useEffect, useState, useMemo } from 'react';
import { Button, Space, Table } from 'antd';
import axios from 'axios';
import moment from 'moment';
import { getCookie } from '../../utils/utilFunctions';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';
import HourglassEmptyIcon from '@material-ui/icons/HourglassEmpty';
import { faQuestionCircle } from '@fortawesome/free-solid-svg-icons';
import Intl from 'intl';
import 'intl/locale-data/jsonp/en-US';
import { useNavigate } from 'react-router-dom';
import { NODE_API_URL, WEB_SOCKET_URL, CELERY_BACKEND_API } from '../../constants/declarations'
import { ScaleLoader } from 'react-spinners';
// import Button from '@material-ui/core/Button';


const TaskTable = () => {
    const [jobs, setJobs] = useState([]);
    const [loading, setLoading] = useState(false);
    const [changesFound, setChangesFound] = useState(false);
    const [filters, setFilters] = useState({});
    const [globalSearchTerm, setGlobalSearchTerm] = useState('');
    const [filteredInfo, setFilteredInfo] = useState({});
    const [sortedInfo, setSortedInfo] = useState({});
    const handleChange = (pagination, filters, sorter) => {
        console.log('Various parameters', pagination, filters, sorter);
        setPagination({
            ...pagination,
        });
        setFilteredInfo(filters);
        setSortedInfo(sorter);
    };
    const clearFilters = () => {
        setFilteredInfo({});
    };
    const clearAll = () => {
        setFilteredInfo({});
        setSortedInfo({});
    };

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

    const [pagination, setPagination] = useState({
        current: 1,
        position: ["topCenter"],
        pageSize: 10, // default number of rows per page
        pageSizeOptions: ['5', '10', '20', '50'], // options for the number of rows per page
        showSizeChanger: true, // show the dropdown to select page size
        });

    const fetchJobs = async (currentPage, searchQuery) => {
        setLoading(true);
        try{
            const response = await fetch(`${NODE_API_URL}/getJobs?q=${searchQuery}&page=${currentPage}`,
                {method: 'POST',
                    headers: {
                    'Content-Type': 'application/json',
                    'Authorization': `Bearer ${jwtToken}`,
                    },
                }
            );
            const data = await response.json();
            data.results.sort((a, b) => a['Created on'] - b['Created on']);
            console.log(data);
            setFilters(data.facets);
            setJobs(data.results);
            setPagination(data.pagination);
            setLoading(false)
            console.log('Fetched jobs:', data);

            // Create a list to store incomplete jobs
            const incompleteTasks = [];

            // Iterate over each task and check if its status is null
            data.results.forEach(task => {
                if (task.Status === null) {
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
        } catch (error) {
            console.error('Error fetching data:', error);
            setLoading(false);
        }
    };

    useEffect(() => {   
        fetchJobs(pagination.page, globalSearchTerm);
    }, []); 

    const handleSearchSubmit = (event) => {
        event.preventDefault();
        fetchJobs(1, globalSearchTerm);
        console.log("Search Handled");
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
                    fetchJobs(pagination.page, globalSearchTerm);
                })
                    .catch(error => {
                        console.error('Error deleting job:', error);
                    });
            })
            .catch(error => {
                console.error('Error deleting job:', error);
            });
    };

    const customSort = (a, b, sortField, sortOrder) => {
        const valueA = a[sortField];
        const valueB = b[sortField];

        // Handle null values
        if (valueA === null && valueB === null) {
            return 0; // Both are null, consider them equal
        }
        if (valueA === null) {
            return 1; // 'a' is null, push it to the end
        }
        if (valueB === null) {
            return -1; // 'b' is null, push it to the end
        }

        // Compare non-null values based on sortOrder
        if (sortOrder === 'ascend') {
            // For numbers
            if (typeof valueA === 'number' && typeof valueB === 'number') {
                return valueA - valueB;
            }
            // For strings
            return valueA.toString().localeCompare(valueB.toString(), 'en', { numeric: true });
        } else { // Descending
            // For numbers
            if (typeof valueA === 'number' && typeof valueB === 'number') {
                return valueB - valueA;
            }
            // For strings
            return valueB.toString().localeCompare(valueA.toString(), 'en', { numeric: true });
        }
    };

    const columns = useMemo(() => {
        if (jobs.length === 0) {
            return [];
        }

        // const baseColumns = Object.keys(jobs[0])
        //     .filter(key => visibleColumns[key])
        //     .map(key => ({
        //         title: key,
        //         dataIndex: key,
        //         key: key,
        //         render: value => {
        //             let res = '';
        //             if (value && typeof value === 'object' && value.label) {
        //                 res = value.label;
        //             } else {
        //                 res = value;
        //             }
        //             return (
        //                 <div 
        //                     data-title={res} 
        //                     className="cell-ellipsis"
        //                     title={res}
        //                 >
        //                     {res}
        //                 </div>
        //             );
        //         }
        //     }));

        const baseColumns = [
            {
                title: 'Description',
                dataIndex: 'Description',
                showSorterTooltip: { target: 'full-header' },
                filters: filters['Description'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo.Description || null,
                // specify the condition of filtering result
                // here is that finding the name started with `value`
                // onFilter: (value, record) => record.Description.indexOf(value) === 0,
                onFilter: (value, record) => record.Description.includes(value),
                sorter: (a, b) => a.Description.length - b.Description.length,
                sortOrder: sortedInfo.columnKey === 'Description' ? sortedInfo.order : null,
                // sortDirections: ['descend'],
                ellipsis: true,
            },
            {
                title: 'Category',
                dataIndex: 'Category',
                showSorterTooltip: { target: 'full-header' },
                filters: filters['Category'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo.Category || null,
                // specify the condition of filtering result
                // here is that finding the name started with `value`
                filterSearch: true,
                // onFilter: (value, record) => record.Category.indexOf(value) === 0,
                onFilter: (value, record) => record.Category.includes(value),
                sorter: (a, b) => customSort(a, b, 'Category', sortedInfo.order),
                sortOrder: sortedInfo.columnKey === 'Category' ? sortedInfo.order : null,
                // sortDirections: ['descend'],
                ellipsis: true,
            },
            {
                title: 'Process',
                dataIndex: 'Process',
                showSorterTooltip: { target: 'full-header' },
                filterSearch: true,
                filters: filters['Process'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo.Process || null,
                // specify the condition of filtering result
                // here is that finding the name started with `value`
                // onFilter: (value, record) => record.Process.indexOf(value) === 0,
                onFilter: (value, record) => record.Process.includes(value),
                sorter: (a, b) => a.Process.length - b.Process.length,
                sortOrder: sortedInfo.columnKey === 'Process' ? sortedInfo.order : null,
                // sortDirections: ['descend'],
                ellipsis: true,
            },
            {
                title: 'Method',
                dataIndex: 'Method',
                showSorterTooltip: { target: 'full-header' },
                filterSearch: true,
                filters: filters['Method'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo.Method || null,
                // specify the condition of filtering result
                // here is that finding the name started with `value`
                // onFilter: (value, record) => record.Method.indexOf(value) === 0,
                onFilter: (value, record) => record.Method.includes(value),
                sorter: (a, b) => a.Method.length - b.Method.length,
                sortOrder: sortedInfo.columnKey === 'Method' ? sortedInfo.order : null,
                // sortDirections: ['descend'],
                ellipsis: true,
            },
            {
                title: 'job_id',
                dataIndex: 'job_id',
                showSorterTooltip: { target: 'full-header' },
                filterSearch: true,
                filters: filters['job_id'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo.job_id || null,
                // specify the condition of filtering result
                // here is that finding the name started with `value`
                // onFilter: (value, record) => record.job_id.indexOf(value) === 0,
                onFilter: (value, record) => record.job_id.includes(value),
                sorter: (a, b) => a.job_id.length - b.job_id.length,
                sortOrder: sortedInfo.columnKey === 'job_id' ? sortedInfo.order : null,
                // sortDirections: ['descend'],
                ellipsis: true,
            },
            {
                title: 'Created on',
                dataIndex: 'Created on',
                defaultSortOrder: 'descend',
                filterSearch: true,
                filters: filters['Created on'].map(filter => ({ text: new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(filter._id).local())) + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo['Created on'] || null,
                // onFilter: (value, record) => record['Created on'].includes(value),
                onFilter: (value, record) => record['Created on'].includes(value),
                sorter: (a, b) => a['Created on'] - b['Created on'],
                sortOrder: sortedInfo.columnKey === 'Created on' ? sortedInfo.order : null,
                ellipsis: true,
                render: value => (
                    <div>
                        {value ? new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(value).local())) : 'N/A'}
                    </div>
                )
            },
            {
                title: 'Completed on',
                dataIndex: 'Completed on',
                defaultSortOrder: 'descend',
                filterSearch: true,
                filters: filters['Completed on'].map(filter => ({ text: new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(filter._id).local())) + "(" + filter.count + ")", value: filter._id })),
                filteredValue: filteredInfo['Completed on'] || null,
                // onFilter: (value, record) => record['Completed on'].includes(value),
                onFilter: (value, record) => record['Completed on'].includes(value),
                sorter: (a, b) => a['Completed on'] - b['Completed on'],
                sortOrder: sortedInfo.columnKey === 'Completed on' ? sortedInfo.order : null,
                ellipsis: true,
                render: value => (
                    <div>
                        {value ? new Intl.DateTimeFormat('en-US', timestampScheme).format(new Date(moment.utc(value).local())) : 'N/A'}
                    </div>
                )
            },   
        ];

        const statusColumn = {
            title: 'Status',
            key: 'Status',
            showSorterTooltip: { target: 'full-header' },
            filterSearch: true,
            filters: filters['Status'].map(filter => ({ text: filter._id + "(" + filter.count + ")", value: filter._id })),
            filteredValue: filteredInfo.Status || null,
            // specify the condition of filtering result
            // here is that finding the name started with `value`
            onFilter: (value, record) => record.Status.includes(value),
            // onFilter: (value, record) => record.Status.indexOf(value) === 0,
            sorter: (a, b) => a.Status.length - b.Status.length,
            sortOrder: sortedInfo.columnKey === 'Status' ? sortedInfo.order : null,
            // sortDirections: ['descend'],
            render: item => {
                return (
                    <div style={{ textAlign: 'center' }}>
                        {item["Status"] === 'Success' ? (
                            <CheckCircleIcon style={{ color: 'green' }} />
                        ) : item["Status"] === 'Failure' ? (
                            <CancelIcon style={{ color: 'red' }} />
                        ) : (
                            <HourglassEmptyIcon style={{ color: 'gray' }} />
                        )}
                    </div>
                );
            }
        };

        const actionColumn = {
            title: 'Actions',
            key: 'actions',
            fixed: 'right',
            width: 120,
            render: item => {
                return (
                    <div className="action-buttons">
                        <Button
                            onClick={() => handleDelete(item["Job ID"])}
                            className="action-button">
                            <FontAwesomeIcon icon={faTrash} />
                        </Button>

                        <Button
                            onClick={() => {
                                if (item["Category"] && item["Category"].toLowerCase() === 'workflow') {
                                    navigate("/mydata/workflowTaskDetails", {
                                        state: {
                                            job_id: item["job_id"],
                                            methodMap: item["Method"],
                                            datasetURL: item["datasetURL"],
                                            description: item["Description"],
                                            process: item["Process"],
                                            output: item["output"],
                                            results: item["results"],
                                            status: item["Status"]
                                        }
                                    });
                                } else {
                                    navigate("/mydata/taskDetails", {
                                        state: {
                                            job_id: item["job_id"],
                                            method: item["Method"],
                                            datasetURL: item["datasetURL"],
                                            description: item["Description"],
                                            process: item["Process"],
                                            output: item["output"],
                                            results: item["results"],
                                            status: item["Status"]
                                        }
                                    });
                                }
                            }}
                            className="action-button">
                            <FontAwesomeIcon icon={faEye} />
                        </Button>
                    </div>
                );
            }
        };

        return [...baseColumns, statusColumn, actionColumn];
    }, [jobs, pagination]);

    useEffect(() => {
        if (!jwtToken)
            navigate('/routing');       
        fetchJobs(1);
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
            <><h1 style={{ textAlign: "left" }}>My Jobs</h1>

            <div className='study-keyword-search'>
                <span className="text-search search-title">Search by text <FontAwesomeIcon icon={faQuestionCircle} /></span>
                <div>
                    <form onSubmit={handleSearchSubmit}>
                    <input
                        type="text"
                        autoComplete="off"
                        className="w-full dark:bg-gray-950 pl-8 form-input-alt h-9 pr-3 focus:shadow-xl"
                        placeholder="Search..."
                        value={globalSearchTerm}
                        onChange={(e) => setGlobalSearchTerm(e.target.value)}
                    />
                    
                    {/* <svg className="absolute left-2.5 text-gray-400 top-1/2 transform -translate-y-1/2" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 32 32">
                        <path d="M30 28.59L22.45 21A11 11 0 1 0 21 22.45L28.59 30zM5 14a9 9 0 1 1 9 9a9 9 0 0 1-9-9z" fill="currentColor"></path>
                    </svg>     */}

                    </form>
                </div>
            </div>

            <div className='table-results'>
                <Space style={{ marginBottom: 16 }}>
                    <Button onClick={clearFilters}>Clear filters</Button>
                    <Button onClick={clearAll}>Clear filters and sorters</Button>
                </Space>

                {loading ? ( 
                    <div className="spinner-container">
                        <ScaleLoader color="#36d7b7" loading={loading} />
                    </div>
                    ) : jobs && jobs.length > 0 ? (                
                    <Table
                        className="table-container"
                        columns={columns}
                        dataSource={jobs}
                        rowKey="job_id"
                        pagination={pagination}
                        onChange={handleChange}
                        showSorterTooltip={{ target: 'sorter-icon' }}
                        onRow={(record,) => {
                            return {
                                onDoubleClick: () => { 
                                    if (record["Category"] && record["Category"].toLowerCase() === 'workflow') {
                                        navigate("/mydata/workflowTaskDetails", {
                                            state: {
                                                job_id: record["job_id"],
                                                methodMap: record["Method"],
                                                datasetURL: record["datasetURL"],
                                                description: record["Description"],
                                                process: record["Process"],
                                                output: record["output"],
                                                results: record["results"],
                                                status: record["Status"]
                                            }
                                        });
                                    } else {
                                        navigate("/mydata/taskDetails", {
                                            state: {
                                                job_id: record["job_id"],
                                                method: record["Method"],
                                                datasetURL: record["datasetURL"],
                                                description: record["Description"],
                                                process: record["Process"],
                                                output: record["output"],
                                                results: record["results"],
                                                status: record["Status"]
                                            }
                                        });
                                    }
                                }
                            };
                        }}
                    />) : (
                        <div>
                            <p>No jobs found.</p>
                        </div>
                    )}
                <div className="pagination-info">
                    <span>* Click <FontAwesomeIcon icon={faTrash} /> to <strong>remove</strong> jobs.</span><br/>
                    <span>
                        * Click <FontAwesomeIcon icon={faEye} /> or <strong>double-click</strong> the row to view job details
                    </span>
                </div>

            </div>
            </>
        );
};

export default TaskTable;