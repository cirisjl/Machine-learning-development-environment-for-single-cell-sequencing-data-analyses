import React, { useState } from 'react';

// This could be replaced with the actual list of datasets you have
const mockDatasets = ['Dataset 1', 'Dataset 2', 'Dataset 3', 'Dataset 4'];

const DatasetSelectionDialog = ({ onSelect, multiple, onClose }) => {
    const [selectedDatasets, setSelectedDatasets] = useState([]);
  
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
  
    return (
      <div className="dialog-backdrop">
        <div className="dialog">
          <h3>Select Datasets</h3>
          <ul>
            {mockDatasets.map((dataset) => (
              <li key={dataset} onClick={() => handleDatasetClick(dataset)}>
                {multiple ? (
                  <input
                    type="checkbox"
                    checked={selectedDatasets.includes(dataset)}
                    readOnly
                  />
                ) : (
                  <input
                    type="radio"
                    name="dataset-selection"
                    checked={selectedDatasets.includes(dataset)}
                    readOnly
                  />
                )}
                {dataset}
              </li>
            ))}
          </ul>
          <button onClick={handleSelectClick}>Select</button>
          <button onClick={onClose}>Cancel</button>
        </div>
      </div>
    );
  };

export default DatasetSelectionDialog;
