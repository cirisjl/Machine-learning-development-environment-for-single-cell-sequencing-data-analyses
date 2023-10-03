import React from 'react';
import UploadDataTaskComponent from './uploadDataTask';
import ValidationTaskComponent from './validationTask';
import GetMetaDataComponent from './getMetadataTask';
import ReviewTaskComponent from './reviewTask';

function MiddleContent({ activeTask, setTaskStatus }) {
  const taskComponents = {
    1: <UploadDataTaskComponent setTaskStatus={setTaskStatus} />,
    2: <ValidationTaskComponent setTaskStatus={setTaskStatus} />,
    3: <GetMetaDataComponent setTaskStatus={setTaskStatus} />,
    4: <ReviewTaskComponent setTaskStatus={setTaskStatus} />,
    // Add other task components here
  };

  return (
    <main>
      {taskComponents[activeTask]}
    </main>
  );
}

export default MiddleContent;
