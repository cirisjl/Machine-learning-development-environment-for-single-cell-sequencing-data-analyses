// TaskInfoTable.js
import React from 'react';
import { Descriptions } from 'antd';

const TaskInfoTable = ({ detail, downloadFile, getFileNameFromURL }) => {
  return (
    <Descriptions title="Task Info" bordered column={1}>
      <Descriptions.Item label="Dataset ID">{detail.datasetId}</Descriptions.Item>
      <Descriptions.Item label="Task Type">{detail.task_type}</Descriptions.Item>
      <Descriptions.Item label="Archive Path">
        <a
          download
          onClick={() => { downloadFile(detail.archive_path) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.archive_path) || 'Not available'}
        </a>
      </Descriptions.Item>
      <Descriptions.Item label="Train Path">
        <a
          download
          onClick={() => { downloadFile(detail.train_path) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.train_path) || 'Not available'}
        </a>
      </Descriptions.Item>
      <Descriptions.Item label="Test Path">
        <a
          download
          onClick={() => { downloadFile(detail.test_path) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.test_path) || 'Not available'}
        </a>
      </Descriptions.Item>
      <Descriptions.Item label="Validation Path">
        <a
          download
          onClick={() => { downloadFile(detail.validation_path) }}
          style={{
            textAlign: 'center',
            cursor: 'pointer',
            textDecoration: 'underline',
            color: 'blue'
          }}>
          {getFileNameFromURL(detail.validation_path) || 'Not available'}
        </a>
      </Descriptions.Item>
    </Descriptions>
  );
};

export default TaskInfoTable;
