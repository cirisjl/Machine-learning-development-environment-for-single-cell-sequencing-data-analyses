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
            accessor: item => {
                console.log("Rendering action cell for row:", item);
                return(
                <div className="action-buttons">
                    <input
                        type="checkbox"
                        style={{ cursor:'pointer' }}
                        onChange={() => onSelectDataset(item["Id"])}
                        checked={selectedDatasets[item["Id"]] === true} // Checked if this datasetId is true in selectedDatasets

                    />
                    <button
                        onClick={() => handleEdit(item["Id"])}
                        className="action-button"
                    >
                        <FontAwesomeIcon icon={faEdit} />
                    </button>
                    <button
                        onClick={() => handleVisualize(item["Id"])}
                        className="action-button"
                    >
                        <FontAwesomeIcon icon={faEye} />
                    </button>
                </div>
                );
            }
        };
          
        return [actionColumn, ...baseColumns];
    }, [data, selectedDatasets]);

    const handleAction1 = (rowData) => {
        console.log('Action 1 clicked for: ', rowData);
        // Implement your Action 1 logic here
      };
    
      const handleAction2 = (rowData) => {
        console.log('Action 2 clicked for: ', rowData);
        // Implement your Action 2 logic here
      };

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
