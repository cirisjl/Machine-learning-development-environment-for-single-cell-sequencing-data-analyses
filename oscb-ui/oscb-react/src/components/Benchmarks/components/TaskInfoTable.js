// TaskInfoTable.js
import React from 'react';
import { Descriptions } from 'antd';

const TaskInfoTable = ({ detail, downloadFile, getFileNameFromURL }) => {
  return (
    <Descriptions title="Task Info" bordered column={1}>
      <Descriptions.Item label="Dataset ID">{detail.datasetId}</Descriptions.Item>
      <Descriptions.Item label="Task Type">{detail.task_type}</Descriptions.Item>
      <Descriptions.Item label="AnnData Path">
        <a
          download
          onClick={() => { downloadFile(detail.adata_path) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.adata_path) || 'Not available'}
        </a>
      </Descriptions.Item>
    </Descriptions>
  );
};

export default TaskInfoTable;
