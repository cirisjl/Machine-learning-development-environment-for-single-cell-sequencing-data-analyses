// DatasetDetailsTable.js
import React from 'react';
import { Descriptions } from 'antd';

const DatasetDetailsTable = ({ datasetDetails, downloadFile, getFileNameFromURL }) => {
  return (
    <Descriptions title="Dataset Details" bordered column={1}>
      <Descriptions.Item label="Dataset ID"><strong>{datasetDetails.Id}</strong></Descriptions.Item>
      <Descriptions.Item label="Title">{datasetDetails.Title}</Descriptions.Item>
      <Descriptions.Item label="Author">{datasetDetails.Author}</Descriptions.Item>
      <Descriptions.Item label="Reference (paper)">{datasetDetails["Reference (paper)"]}</Descriptions.Item>
      <Descriptions.Item label="DOI"><a href={datasetDetails.DOI}>{datasetDetails.DOI}</a></Descriptions.Item>
      <Descriptions.Item label="Abstract">{datasetDetails.Abstract}</Descriptions.Item>
      <Descriptions.Item label="Original Downloads"><a href={datasetDetails.Downloads}>{datasetDetails.Downloads}</a></Descriptions.Item>
      <Descriptions.Item label="Species">{datasetDetails.Species?.label || datasetDetails.Species}</Descriptions.Item>
      <Descriptions.Item label="Sample Type">{datasetDetails["Sample Type"]?.label}</Descriptions.Item>
      <Descriptions.Item label="Cell Count Estimate">{datasetDetails["Cell Count Estimate"]}</Descriptions.Item>
      <Descriptions.Item label="Organ Part">{datasetDetails["Organ Part"]?.label || datasetDetails["Organ Part"]}</Descriptions.Item>
      <Descriptions.Item label="Selected Cell Types">
        {datasetDetails["Selected Cell Types"]?.value.join(', ')}
      </Descriptions.Item>
      <Descriptions.Item label="Library Construction Method">{datasetDetails["Library Construction Method"]?.label}</Descriptions.Item>
      <Descriptions.Item label="Nucleic Acid Source">{datasetDetails["Nucleic Acid Source"]?.label}</Descriptions.Item>
      <Descriptions.Item label="Analysis Protocol">{datasetDetails["Analysis Protocol"]}</Descriptions.Item>
      <Descriptions.Item label="Disease Status (Specimen)">
        {Array.isArray(datasetDetails["Disease Status (Specimen)"]) 
          ? datasetDetails["Disease Status (Specimen)"].map((status, index) => (
              <span key={index}>{status.label}{index < datasetDetails["Disease Status (Specimen)"].length - 1 ? ', ' : ''}</span>
            ))
          : datasetDetails["Disease Status (Specimen)"]?.label || datasetDetails["Disease Status (Specimen)"]
        }
      </Descriptions.Item>
      <Descriptions.Item label="Development Stage">
        {Array.isArray(datasetDetails["Development Stage"])
          ? datasetDetails["Development Stage"].map((status, index) => (
            <span key={index}>{status.label}{index < datasetDetails["Development Stage"].length - 1 ? ', ' : ''}</span>
          ))
          : datasetDetails["Development Stage"]?.label || datasetDetails["Development Stage"]
        }
      </Descriptions.Item>
      <Descriptions.Item label="Source"><a href={datasetDetails["Source"]?.label}>{datasetDetails["Source"]?.label}</a></Descriptions.Item>
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
