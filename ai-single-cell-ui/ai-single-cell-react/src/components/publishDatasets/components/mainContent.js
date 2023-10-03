import React from 'react';
import UploadDataTaskComponent from './uploadDataTask';
import ValidationTaskComponent from './validationTask';
import GetMetaDataComponent from './getMetadataTask';
import ReviewTaskComponent from './reviewTask';
import QualityControlTaskComponent from './qualityControlTask';
import TaskBuilderTaskComponent from './taskBuilderTask';
import BenchmarksTaskComponent from './benchmarksTask';

function MiddleContent({ activeTask, setTaskStatus }) {
  const taskComponents = {
    1: <UploadDataTaskComponent setTaskStatus={setTaskStatus} />,
    2: <ValidationTaskComponent setTaskStatus={setTaskStatus} />,
    3: <QualityControlTaskComponent setTaskStatus={setTaskStatus} />,
    4: <GetMetaDataComponent setTaskStatus={setTaskStatus} />,
    5: <TaskBuilderTaskComponent setTaskStatus={setTaskStatus} />,
    6: <BenchmarksTaskComponent setTaskStatus={setTaskStatus} />,
    7: <ReviewTaskComponent setTaskStatus={setTaskStatus} />,
    // Add other task components here
  };

  return (
    <main>
      {taskComponents[activeTask]}
    </main>
  );
}

export default MiddleContent;
