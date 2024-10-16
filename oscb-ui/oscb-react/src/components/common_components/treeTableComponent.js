import { faEdit, faEye } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React, { useState, useMemo } from 'react';
import { Table } from 'antd';
import Checkbox from '@material-ui/core/Checkbox';
import Menu from '@material-ui/core/Menu';
import Button from '@material-ui/core/Button';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import MenuItem from '@material-ui/core/MenuItem';
import ListItemText from '@material-ui/core/ListItemText';

const TreeTable = ({ data, onSelectDataset, selectedDatasets, multiple, pagination }) => {
    const [anchorEl, setAnchorEl] = useState(null);

    // Destructure the pagination object for easier access to its properties
    // const { page, pageSize, totalCount } = pagination;

    // // Calculate the starting and ending result numbers for the current page
    // const startResult = (page - 1) * pageSize + 1;
    // const endResult = Math.min(page * pageSize, totalCount); // Ensure not to exceed totalCount

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
        console.log(visibleColumns[column]);
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

    const handleVisualize = (benchmarksId) => {
        console.log("Visualizing benchmarksId: ", benchmarksId);
        window.open(`/benchmarks/viewDetails?benchmarksId=${benchmarksId}`, '_blank');
    };

    const columns = useMemo(() => {
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
                    return (
                        <div 
                            data-title={res} 
                            className="cell-ellipsis"
                            title={res}
                        >
                            {res}
                        </div>
                    );
                }
            }));

        const actionColumn = {
            title: 'Actions',
            key: 'actions',
            render: item => {
                return (
                    <div className="action-buttons">
                        <Checkbox
                            style={{ cursor: 'pointer' }}
                            onChange={() => onSelectDataset(item)}
                            checked={!!selectedDatasets[item["Benchmarks ID"]]}
                        />
                        <Button
                            onClick={() => handleVisualize(item["Benchmarks ID"])}
                            className="action-button"
                        >
                            <FontAwesomeIcon icon={faEye} />
                        </Button>
                    </div>
                );
            }
        };

        return [actionColumn, ...baseColumns];
    }, [data, selectedDatasets, visibleColumns]);

    return (
        <div>
            {/* Dropdown for editing columns */}
            <div className="dropdown">
                <div className='total-results-count'>
                    {/* Do not remove this div. It will remove styles. If you want to remove this div, Add alternate styles to the edit columns without breaking it. */}
                </div>
                <Button variant="contained" onClick={handleMenuClick}>
                    Edit Columns
                </Button>
                <Menu
                    id="edit-columns-menu"
                    anchorEl={anchorEl}
                    keepMounted
                    open={Boolean(anchorEl)}
                    onClose={handleMenuClose}
                    style={{ height: "400px" }} // Adjusts the vertical position
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
            />
        </div>
    );
};

export default TreeTable;
