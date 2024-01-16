import React from 'react';
import { useTable } from 'react-table';


const ResultsTable = ({ data }) => {
    const columns = React.useMemo(() => {
        if (data.length === 0) {
            return [];
        }

        // Extracting the first item's keys to determine columns
        const sampleItem = data[0];

        return Object.keys(sampleItem).map(key => {
            return {
                Header: key,
                // Customize the accessor to handle nested objects
                accessor: item => {
                    const value = item[key];
                    if (value && typeof value === 'object' && value.label) {
                        return value.label;
                    }
                    return value;
                }
            };
        });
    }, [data]);

    const {
        getTableProps,
        getTableBodyProps,
        headerGroups,
        rows,
        prepareRow,
    } = useTable({ columns, data });

    return (
        <table {...getTableProps()}>
            <thead>
                {headerGroups.map(headerGroup => (
                    <tr {...headerGroup.getHeaderGroupProps()}>
                        {headerGroup.headers.map(column => (
                            <th {...column.getHeaderProps()}>{column.render('Header')}</th>
                        ))}
                    </tr>
                ))}
            </thead>
            <tbody {...getTableBodyProps()}>
                {rows.map(row => {
                    prepareRow(row);
                    return (
                        <tr {...row.getRowProps()}>
                            {row.cells.map(cell => {
                                return <td {...cell.getCellProps()}>{cell.value}</td>;
                            })}
                        </tr>
                    );
                })}
            </tbody>
        </table>
    );
};

export default ResultsTable;
