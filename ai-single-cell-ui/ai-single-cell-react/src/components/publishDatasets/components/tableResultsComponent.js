import { faEdit, faEye } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React from 'react';
import { useTable, useRowSelect } from 'react-table';


const ResultsTable = ({ data, onSelectDataset, selectedDatasets }) => {

    const handleEdit = (dataset) => {
        console.log("Editing dataset: ", dataset);
        // Implement your edit logic here
    };
    
    const handleVisualize = (dataset) => {
        console.log("Visualizing dataset: ", dataset);
        // Implement your visualization logic here
    };

    const columns = React.useMemo(() => {
        if (data.length === 0) {
            return [];
        }

        const baseColumns = Object.keys(data[0]).map(key => ({
            Header: key,
            accessor: item => {
                const value = item[key];
                if (value && typeof value === 'object' && value.label) {
                    return value.label;
                }
                return value;
            }
        }));

        const actionColumn = {
            id: 'actions',
            Header: 'Actions',
            Cell: ({ row }) => {
                console.log("Rendering action cell for row:", row);
                <div className="action-buttons">
                    <input
                        type="checkbox"
                        // {...row.getToggleRowSelectedProps()}
                        checked={selectedDatasets[row.id] === true}
                        onChange={() => onSelectDataset(row.id)}
                        style={{ marginRight: '5px' }}
                    />
                    <button
                        onClick={() => handleEdit(row.original)}
                        className="action-button"
                    >Edit
                        {/* <FontAwesomeIcon icon={faEdit} /> */}
                    </button>
                    <button
                        onClick={() => handleVisualize(row.original)}
                        className="action-button"
                    >
                        <FontAwesomeIcon icon={faEye} />
                    </button>
                </div>
            }
        };

        return [actionColumn, ...baseColumns];
    }, [data, selectedDatasets]);


    const {
        getTableProps,
        getTableBodyProps,
        headerGroups,
        rows,
        prepareRow,
    } = useTable({ columns, data },useRowSelect);

    return (
        <table {...getTableProps()} className="table-container">
        <thead>
                {headerGroups.map((headerGroup, index) => (
                    <tr {...headerGroup.getHeaderGroupProps()} key={index}>
                        {headerGroup.headers.map((column, colIndex) => (
                            <th {...column.getHeaderProps()} key={colIndex}>{column.render('Header')}</th>
                        ))}
                    </tr>
                ))}
            </thead>
            <tbody {...getTableBodyProps()}>
                {rows.map((row, rowIndex) => {
                    prepareRow(row);
                    return (
                        <tr {...row.getRowProps()} key={rowIndex}>
                            {row.cells.map((cell, cellIndex) => {
                                return <td {...cell.getCellProps()} key={cellIndex}>{cell.value}</td>;
                            })}
                        </tr>
                    );
                })}
            </tbody>
        </table>
    );
};

export default ResultsTable;
