import { faEdit, faEye, faTrash } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React, {useState} from 'react';
import { useTable, useRowSelect } from 'react-table';
import Button from '@material-ui/core/Button';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Checkbox from '@material-ui/core/Checkbox';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import ListItemText from '@material-ui/core/ListItemText';
import { Table } from 'antd';
import axios from 'axios';
import { CELERY_BACKEND_API } from '../../../constants/declarations'


const ResultsTable = ({ data, onSelectDataset, selectedDatasets, multiple, pagination, onSelectSubItem, showEdit=false, showDelete=false }) => {

    const [anchorEl, setAnchorEl] = useState(null);
    const [subItemsData, setSubItemsData] = useState({});

    const [paginationState, setPagination] = useState({
        current: 1,
        position: ["topCenter"],
        pageSize: 10, // default number of rows per page
        pageSizeOptions: ['5', '10', '20', '50'], // options for the number of rows per page
        showSizeChanger: true, // show the dropdown to select page size
      });

    // Handle pagination change (page number and page size)
    const handleTableChange = (pagination) => {
        setPagination({
        ...pagination,
        });
    };
    // Destructure the pagination object for easier access to its properties
    // const { page, pageSize, totalCount } = pagination;

    // Calculate the starting and ending result numbers for the current page
    // const startResult = (page - 1) * pageSize + 1;
    // const endResult = Math.min(page * pageSize, totalCount); // Ensure not to exceed totalCount

    const [visibleColumns, setVisibleColumns] = useState({
        'Benchmarks ID': true,
        'Task': true,
        'Title': true,
        'Category': true,
        'Species': true,
        'Cell Count Estimate': true,
        'Organ Part': true,
        'Development Stage': false, // Optional initially not visible
        'Author': false, // Optional initially not visible
        'Submission Date': false, // Optional initially not visible
        'Source': false, // Optional initially not visible
    });

    const handleMenuClick = (event) => {
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const toggleColumnVisibility = (column) => {
        console.log("Toggle Column visibility");
        setVisibleColumns(prevVisibleColumns => ({
            ...prevVisibleColumns,
            [column]: !prevVisibleColumns[column],
        }));
    };

    const resetColumnVisibility = () => {
        setVisibleColumns({
            'Benchmarks ID': true,
            'Task': true,
            'Title': true,
            'Species': true,
            'Cell Count Estimate': true,
            'Organ Part': true,
            'Development Stage': false,
            'Author': false,
            'Submission Date': false,
            'Source': false,
        });
    };
    
    const handleVisualize = (datasetId) => {
        console.log("Dataset Id: ", datasetId);
        window.open(`/mydata/view-dataset-info?datasetId=${datasetId}`, '_blank');
    };

    const handleEdit = (datasetId) => {
        console.log("Edit dataset Id: ", datasetId);
        window.open(`/mydata/edit-dataset-info?datasetId=${datasetId}`, '_blank');
    };

    const handleDelete = (datasetId) => {
        console.log("Delete dataset Id: ", datasetId);
        // Implement your visualization logic here
    };

    const columns = React.useMemo(() => {
        if (data.length === 0) {
            return [];
        }

        const baseColumns = Object.keys(data[0])
        .filter(key => visibleColumns[key])
        .map(key => ({
            title: key,
            dataIndex: key,
            key: key,
            render: value => {
                let res = '';
                if (value && typeof value === 'object' && value.label) {
                    res = value.label;
                } else {
                    res = value;
                }
                return(
                    <div 
                        data-title={res} 
                        className="cell-ellipsis"
                        title={res}
                        >
                        {res}
                    </div>
                )
                
            }
        }));

        const actionColumn = {
            title: 'Actions',
            key: 'actions',
            render: item => {
                return(
                <div className="action-buttons">
                    <input
                        type="checkbox"
                        style={{ cursor:'pointer' }}
                        onChange={() => onSelectDataset(item)}
                        checked={!!selectedDatasets[item["Id"]]}
                        // disabled={isDisabled() && !isSelected(item["Id"])} // Disable if multiple is false and a dataset is already selecte
                    />
                    <button
                        onClick={() => handleEdit(item["Id"])}
                        className="action-button">
                        <FontAwesomeIcon icon={faEdit} />
                    </button>

                    <button
                        onClick={() => handleDelete(item["Id"], item)}
                        className="action-button">
                        <FontAwesomeIcon icon={faTrash} />
                    </button>

                    <button
                        onClick={() => handleVisualize(item["Id"])}
                        className="action-button">
                        <FontAwesomeIcon icon={faEye} />
                    </button>
                </div>
                );
            }
        };
          
        return [actionColumn, ...baseColumns];
    }, [data, selectedDatasets, visibleColumns]);

    const fetchSubItems = async (process_ids) => {
        const process_ids_key = process_ids.join(',');
        
        // Check if data already exists for this process_ids_key
        if (subItemsData[process_ids_key]) {
            return; // Data already exists, no need to fetch again
        }
    
        try {
            const response = await axios.post(CELERY_BACKEND_API + '/getPreProcessResults', { process_ids: process_ids, record_type: "table" });
            setSubItemsData(prevData => ({
                ...prevData,
                [process_ids_key]: response.data // Store the result with concatenated process_ids as key
            }));
        } catch (error) {
            console.error("Error fetching sub-items:", error);
        }
    };
    

    const expandedRowRender = (record) => {

        const process_ids_key = record.process_ids.join(',');

        // Fetch sub-items only if not already fetched
        if (!subItemsData[process_ids_key]) {
            fetchSubItems(record.process_ids);
        }
        const subColumns = [
            {
                title: 'Description',
                dataIndex: 'description',
                key: 'description',
            },
            {
                title: 'Stage',
                dataIndex: 'stage',
                key: 'stage',
            },
            {
                title: 'Process',
                dataIndex: 'process',
                key: 'process',
            },
            {
                title: 'Method',
                dataIndex: 'method',
                key: 'method',
            },
            {
                title: 'nCells',
                dataIndex: 'nCells',
                key: 'nCells',
            },
            {
                title: 'Action',
                key: 'operation',
                render: (text, subRecord) => (
                    <Checkbox 
                    onChange={() => onSelectSubItem(record, subRecord)} 
                    checked={selectedDatasets[record.Id]?.selectedSubItem?.process_id === subRecord.process_id}
                    />
                ),
            },
        ];

        const subData = subItemsData[process_ids_key] || [];
        return <Table columns={subColumns} dataSource={subData} pagination={false} />;
    };

    return (
        <div>
            {/* Dropdown for editing columns */}
            <div className="dropdown">
                <div className='total-results-count'>
                    {/* Do not remove this div. It will remove styles. If you want to remove this div, Add alternate styles to the edit columns without breaking it. */}
                </div>
                <Button variant="contained" aria-controls="edit-columns-menu" aria-haspopup="true" onClick={handleMenuClick} >
                    Edit Columns
                </Button>
                <Menu
                    id="edit-columns-menu"
                    anchorEl={anchorEl}
                    keepMounted
                    open={Boolean(anchorEl)}
                    onClose={handleMenuClose}
                    style={{height: "400px" }} // Adjusts the vertical position
                    anchorOrigin={{
                        vertical: 'bottom', 
                        horizontal: 'right', 
                      }}
                      transformOrigin={{
                        vertical: 'top', 
                        horizontal: 'right',
                      }}
                      getContentAnchorEl={null} // This will make anchorOrigin work as expected
                >
                    <FormGroup>
                        {Object.keys(visibleColumns).map((column) => (
                           <MenuItem key={column} onClick={(event) => event.stopPropagation()}>
                           <FormControlLabel
                               control={
                                   <Checkbox
                                       checked={visibleColumns[column]}
                                       onChange={() => toggleColumnVisibility(column)}
                                       onClick={(event) => event.stopPropagation()} // Prevent triggering the menu item's onClick
                                   />
                               }
                               label={column}
                               // Remove the onClick here to avoid overriding Checkbox's behavior
                           />
                       </MenuItem>
                        ))}
                        {/* Reset Menu Item */}
                        <MenuItem onClick={resetColumnVisibility}>
                            <ListItemText primary="Reset" />
                        </MenuItem>
                    </FormGroup>
                </Menu>
            </div>

            <Table
                className="table-container"
                columns={columns}
                dataSource={data}
                rowKey="Id"
                pagination={paginationState}
                onChange={handleTableChange}
                expandable={{
                    expandedRowRender: (record) => {
                        fetchSubItems(record.process_ids); // Pass process_ids array
                        return expandedRowRender(record);
                    },
                    rowExpandable: (record) => Array.isArray(record.process_ids) && record.process_ids.length > 0,
                }}
            />

        </div>
    );
};

export default ResultsTable;
