// DatasetDetailsTable.js
import React from 'react';
import { Descriptions } from 'antd';

const DatasetDetailsTable = ({ datasetDetails, downloadFile, getFileNameFromURL }) => {
  return (
    <Descriptions title="Dataset Details" bordered column={1}>
      <Descriptions.Item label="Dataset ID">{datasetDetails.Id}</Descriptions.Item>
      <Descriptions.Item label="Title">{datasetDetails.Title}</Descriptions.Item>
      <Descriptions.Item label="Author">{datasetDetails.Author}</Descriptions.Item>
      <Descriptions.Item label="Species">{datasetDetails.Species?.label || datasetDetails.Species}</Descriptions.Item>
      <Descriptions.Item label="Cell Count Estimate">{datasetDetails["Cell Count Estimate"]}</Descriptions.Item>
      <Descriptions.Item label="Organ Part">{datasetDetails["Organ Part"]?.label || datasetDetails["Organ Part"]}</Descriptions.Item>
      <Descriptions.Item label="Selected Cell Types">
        {Array.isArray(datasetDetails["Selected Cell Types"]) 
          ? datasetDetails["Selected Cell Types"].map((cellType, index) => (
              <span key={index}>{cellType.label}{index < datasetDetails["Selected Cell Types"].length - 1 ? ', ' : ''}</span>
            ))
          : datasetDetails["Selected Cell Types"]?.label || datasetDetails["Selected Cell Types"]
        }
      </Descriptions.Item>
      <Descriptions.Item label="Disease Status (Specimen)">
        {Array.isArray(datasetDetails["Disease Status (Specimen)"]) 
          ? datasetDetails["Disease Status (Specimen)"].map((status, index) => (
              <span key={index}>{status.label}{index < datasetDetails["Disease Status (Specimen)"].length - 1 ? ', ' : ''}</span>
            ))
          : datasetDetails["Disease Status (Specimen)"]?.label || datasetDetails["Disease Status (Specimen)"]
        }
      </Descriptions.Item>
      <Descriptions.Item label="Submission Date">{datasetDetails["Submission Date"]}</Descriptions.Item>
      <Descriptions.Item label="AnnData File">
        <a
          download
          onClick={() => { downloadFile(datasetDetails["adata_path"]) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(datasetDetails["adata_path"]) || 'Not available'}
        </a>
      </Descriptions.Item>
    </Descriptions>
  );
};

export default DatasetDetailsTable;
