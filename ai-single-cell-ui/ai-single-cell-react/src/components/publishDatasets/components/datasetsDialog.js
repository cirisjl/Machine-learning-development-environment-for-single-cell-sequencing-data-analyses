import React, { useState, useMemo } from 'react';

// This could be replaced with the actual list of datasets you have
const mockDatasets = ['Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'];

const DatasetSelectionDialog = ({ datasets, onSelect, multiple, onClose , isVisible }) => {
    const [selectedDatasets, setSelectedDatasets] = useState([]);

    const dialogStyle = {
        display: isVisible ? 'block' : 'none',
        // ... other styles
    };
  
    const handleSelectClick = () => {
      onSelect(selectedDatasets);
      onClose(); // Close the dialog after selection
    };
  
    const handleDatasetClick = (dataset) => {
      if (multiple) {
        // For multiple selection, toggle datasets in the array
        setSelectedDatasets(selectedDatasets.includes(dataset)
          ? selectedDatasets.filter((d) => d !== dataset)
          : [...selectedDatasets, dataset]);
      } else {
        // For single selection, set the dataset directly
        setSelectedDatasets([dataset]);
      }
    };

    const [data, setData] = useState({
      initial_data: datasets,
    });
    const [sortColumn, setSortColumn] = useState(null);
    const [sortDirection, setSortDirection] = useState(1);
    const [visibleColumns, setVisibleColumns] = useState([
      'Species',
      'Library Construction Method',
      'Anatomical Entity',
      'Disease Status (Donor)',
      'Cell Count Estimate',
    ]);
  
    const [filterValues, setFilterValues] = useState({});
    const [searchTerm, setSearchTerm] = useState('');
    const [currentPage, setCurrentPage] = useState(1);
    const [itemsPerPage, setItemsPerPage] = useState(25);
    const [totalPages, setTotalPages] = useState(1);
    const [globalSearchTerm, setGlobalSearchTerm] = useState('');
    const [exportedData, setExportedData] = useState([]);
    const [error, setError] = useState(null);
  
    const requiredColumns = [
      'Species',
      'Library Construction Method',
      'Anatomical Entity',
      'Disease Status (Donor)',
      'Cell Count Estimate',
    ];
  
    const columnColors = {
      Species: 'btn-info',
      'Library Construction Method': 'btn-info',
      'Anatomical Entity': 'btn-info',
      'Disease Status (Donor)': 'btn-info',
      'Cell Count Estimate': 'btn-info',
    };
  
    const handleSort = (columnName) => {
      try {
        const newDirection = columnName === sortColumn ? -sortDirection : 1;
  
        const sortedData = [...data.initial_data].sort((a, b) => {
          if (a[columnName] < b[columnName]) return -newDirection;
          if (a[columnName] > b[columnName]) return newDirection;
          return 0;
        });
  
        setData({
          initial_data: sortedData,
        });
  
        setSortColumn(columnName);
        setSortDirection(newDirection);
        setError(null);
      } catch (error) {
        setError(`Error sorting data: ${error.message}`);
      }
    };
  
    const handleColumnToggle = (columnName) => {
      setVisibleColumns((prevVisibleColumns) => {
        if (prevVisibleColumns.includes(columnName)) {
          return prevVisibleColumns.filter((col) => col !== columnName);
        } else {
          return [...prevVisibleColumns, columnName];
        }
      });
    };
  
    const handleFilterChange = (column, value) => {
      setFilterValues((prevFilterValues) => ({
        ...prevFilterValues,
        [column]: value,
      }));
      // Clear global search term when a column filter is changed
      setGlobalSearchTerm('');
    };
  
    const handleGlobalSearchChange = (value) => {
      setGlobalSearchTerm(value);
    };
  
    const handleResetFilters = () => {
      setGlobalSearchTerm('');
      setFilterValues({});
    };
  
    const applySearch = (row) => {
      return (
        (globalSearchTerm === '' ||
          requiredColumns.some(
            (column) =>
              row[column]
                .toString()
                .toLowerCase()
                .includes(globalSearchTerm.toLowerCase())
          )) &&
        requiredColumns.every((column) => {
          const value = row[column];
          const filterValue = filterValues[column];
  
          if (
            filterValue === undefined ||
            filterValue === null ||
            filterValue === ''
          ) {
            return true; // No filter for this column
          }
  
          if (typeof value === 'string' && typeof filterValue === 'string') {
            // Case-insensitive search for strings
            return value.toLowerCase().includes(filterValue.toLowerCase());
          } else if (
            typeof value === 'number' &&
            typeof filterValue === 'string'
          ) {
            // Convert numbers to string and then search
            return value.toString().includes(filterValue);
          }
  
          // Non-string and non-number types are not searchable
          return false;
        })
      );
    };
  
    const renderCell = (column, value, index) => {
      const uniqueKey = `cell_${index}_${column}`;
  
      if (
        column === 'Species' ||
        column === 'Anatomical Entity' ||
        column === 'Disease Status (Donor)'
      ) {
        const values = value.label.split(',');
        const count = values.length;
        var k = '';
        if (column === 'Anatomical Entity') {
          k = 'Anatomic';
        } else if (column === 'Disease Status (Donor)') {
          k = 'Disease';
        } else {
          k = 'Species';
        }
  
        if (count === 1) {
          return <td key={uniqueKey}>{value}</td>;
        } else {
          const highlightClass =
            column === 'Species' ||
            column === 'Anatomical Entity' ||
            column === 'Disease Status (Donor)'
              ? 'highlight'
              : '';
  
          return (
            <td
              key={uniqueKey}
              data-toggle="tooltip"
              data-placement="top"
              title={value}
              className={highlightClass}
            >
              <span className='abc'>{`${k} (${count})`}</span>
            </td>
          );
        }
      } else {
        return (
          <td
            key={uniqueKey}
            data-toggle="tooltip"
            data-placement="top"
            title={value}
          >
            {value}
          </td>
        );
      }
    };
  
    const handleItemsPerPageChange = (value) => {
      setItemsPerPage(value);
      setCurrentPage(1);
    };
  
    const downloadCSV = (csvContent, fileName) => {
      const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
      const link = document.createElement('a');
  
      if (link.download !== undefined) {
        const url = URL.createObjectURL(blob);
        link.setAttribute('href', url);
        link.setAttribute('download', fileName);
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      } else {
        alert(
          'Your browser does not support the download functionality. Please try a different browser.'
        );
      }
    };
  
    const formatCount = (count) => {
      if (count >= 1e6) {
        return `${(count / 1e6).toFixed(1)}M`;
      } else if (count >= 1e3) {
        return `${(count / 1e3).toFixed(1)}k`;
      } else {
        return count.toString();
      }
    };
  
    const handleExport = () => {
      try {
        const exportedData = data.initial_data
          .filter(applySearch)
          .map((result) =>
            visibleColumns.map((column) => result[column]).join(',')
          );
  
        setExportedData(exportedData);
  
        downloadCSV(exportedData.join('\n'), 'exported_data.csv');
        setError(null);
      } catch (error) {
        setError(`Error exporting data: ${error.message}`);
      }
    };
  
    const paginatedData = useMemo(() => {
      const startIndex = (currentPage - 1) * itemsPerPage;
      const endIndex = startIndex + itemsPerPage;
      setTotalPages(Math.ceil(data.initial_data.length / itemsPerPage));
  
      return data.initial_data.slice(startIndex, endIndex);
    }, [data.initial_data, currentPage, itemsPerPage]);
  
  
    return (
        <div style={dialogStyle} className="dialog-container">
            <div className="dialog-backdrop" onClick={onClose} />
            <div className="dialog">
                <h3>Select Datasets</h3>
                <header>
                  <div className="d-flex justify-content-between align-items-center">
                    <div className="text-center flex-grow-1">
                      <h1>Explore Data</h1>
                    </div>
                    <div className="d-flex align-items-center ">
                      <input
                        type="text"
                        placeholder="🔍 Search by text"
                        value={globalSearchTerm}
                        onChange={(e) => handleGlobalSearchChange(e.target.value)}
                        className="ml-5"
                        style={{
                          marginLeft: '550px',
                          padding: '5px',
                        }}
                      />
                      <button
                        className="btn btn-secondary ml-2"
                        onClick={handleResetFilters}
                      >
                        Reset
                      </button>
                    </div>
                  </div>
                </header>

                <div className="border p-3 mb-3">
                  <div className="d-flex justify-content-between align-items-center">
                    {requiredColumns.map((column) => (
                      <div key={column} className="text-center">
                        <span className="mr-2">
                          {column !== 'CellCountEstimate' ? (
                            <>
                              {formatCount(
                                new Set(
                                  data.initial_data.map((row) => row[column])
                                ).size
                              )}{' '}
                              {column}
                            </>
                          ) : (
                            `Total ${column}: ${formatCount(
                              data.initial_data.reduce(
                                (sum, row) =>
                                  sum + (row[column] ? parseInt(row[column], 10) : 0),
                                0
                              )
                            )}`
                          )}
                        </span>
                      </div>
                    ))}
                    <button className="btn btn-success ml-2" onClick={handleExport}>
                      Export
                    </button>
                  </div>
                </div>

                <div className="mb-3">
                  <label className="form-label">Select Columns:</label>
                  <div className="btn-group" role="group">
                    {requiredColumns.map((column) => (
                      <button
                        key={column}
                        type="button"
                        className={`btn ${columnColors[column]} ${
                          visibleColumns.includes(column) ? 'active' : ''
                        }`}
                        style={{
                          marginRight: '8px',
                          backgroundColor: visibleColumns.includes(column)
                            ? '#007bff'
                            : 'transparent',
                          borderColor: visibleColumns.includes(column)
                            ? '#007bff'
                            : '',
                        }}
                        onClick={() => handleColumnToggle(column)}
                      >
                        {column}
                      </button>
                    ))}
                    <div className="d-flex align-items-center ">
                      <input
                        type="text"
                        placeholder="🔍 Search by text"
                        value={globalSearchTerm}
                        onChange={(e) => handleGlobalSearchChange(e.target.value)}
                        className="ml-5"
                        style={{
                          marginLeft: '550px',
                          padding: '5px',
                        }}
                      />
                      <button
                        className="btn btn-secondary ml-2"
                        onClick={handleResetFilters}
                      >
                        Reset
                      </button>
                    </div>
                  </div>
                </div>

                <div className="mb-3">
                  <table className="table table-striped">
                    <thead>
                      <tr>
                        {visibleColumns.map((column) => (
                          <th
                            key={column}
                            style={{ cursor: 'pointer' }}
                            onClick={() => handleSort(column)}
                          >
                            {column}
                            {sortColumn === column && (
                              <span>{sortDirection === 1 ? '↑' : '↓'}</span>
                            )}
                          </th>
                        ))}
                      </tr>
                      <tr>
                        {visibleColumns.map((column) => (
                          <th key={column}>
                            <input
                              type="text"
                              placeholder={`Search ${column}`}
                              value={filterValues[column] || ''}
                              onChange={(e) => handleFilterChange(column, e.target.value)}
                            />
                          </th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {paginatedData
                        .filter(applySearch)
                        .map((result, index) => (
                          <tr key={`row_${index}`}>
                            {visibleColumns.map((column) =>
                              renderCell(column, result[column], index)
                            )}
                          </tr>
                        ))}
                    </tbody>
                  </table>
                </div>

                <div className="pagination">
                  <label className="form-label mr-2">Items Per Page:</label>
                  <select
                    className="form-select"
                    onChange={(e) => handleItemsPerPageChange(Number(e.target.value))}
                    value={itemsPerPage}
                    style={{ width: '40px' }}
                  >
                    <option value={10}>10</option>
                    <option value={20}>20</option>
                    <option value={30}>30</option>
                    <option value={40}>40</option>
                  </select>
                  <button
                    className="btn btn-secondary ml-2"
                    onClick={() =>
                      setCurrentPage((prevPage) => Math.max(prevPage - 1, 1))
                    }
                    disabled={currentPage === 1}
                  >
                    <span>&#8592;</span>
                  </button>
                  <span className="ml-2">
                    Page {currentPage} of {totalPages}
                  </span>
                  <button
                    className="btn btn-secondary ml-2"
                    onClick={() =>
                      setCurrentPage((prevPage) =>
                        Math.min(prevPage + 1, totalPages)
                      )
                    }
                    disabled={currentPage === totalPages}
                  >
                    <span>&#8594;</span>
                  </button>
                </div>
            </div>
        </div>
    );
  };

export default DatasetSelectionDialog;
