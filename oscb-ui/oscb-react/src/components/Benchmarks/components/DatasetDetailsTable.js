// DatasetDetailsTable.js
import React from 'react';
import { Descriptions } from 'antd';

const DatasetDetailsTable = ({ detail, downloadFile, getFileNameFromURL }) => {
  return (
    <Descriptions title="Dataset Details" bordered column={1}>
      <Descriptions.Item label="Dataset ID">{detail.datasetDetails.Id}</Descriptions.Item>
      <Descriptions.Item label="Title">{detail.datasetDetails.Title}</Descriptions.Item>
      <Descriptions.Item label="Author">{detail.datasetDetails.Author}</Descriptions.Item>
      <Descriptions.Item label="Species">{detail.datasetDetails.Species?.label || detail.datasetDetails.Species}</Descriptions.Item>
      <Descriptions.Item label="Cell Count Estimate">{detail.datasetDetails["Cell Count Estimate"]}</Descriptions.Item>
      <Descriptions.Item label="Organ Part">{detail.datasetDetails["Organ Part"]?.label || detail.datasetDetails["Organ Part"]}</Descriptions.Item>
      <Descriptions.Item label="Selected Cell Types">{detail.datasetDetails["Selected Cell Types"]?.label || detail.datasetDetails["Selected Cell Types"]}</Descriptions.Item>
      <Descriptions.Item label="Disease Status (Specimen)">{detail.datasetDetails["Disease Status (Specimen)"]?.label || detail.datasetDetails["Disease Status (Specimen)"]}</Descriptions.Item>
      <Descriptions.Item label="Submission Date">{detail.datasetDetails["Submission Date"]}</Descriptions.Item>
      <Descriptions.Item label="AnnData File">
        <a
          download
          onClick={() => { downloadFile(detail.datasetDetails["adata_path"]) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.datasetDetails["adata_path"]) || 'Not available'}
        </a>
      </Descriptions.Item>
    </Descriptions>
  );
};

export default DatasetDetailsTable;
