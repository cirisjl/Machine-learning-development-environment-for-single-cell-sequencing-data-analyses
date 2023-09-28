import React from 'react';
import UploadDataTaskComponent from './uploadDataTask';
import ValidationTaskComponent from './validationTask';
import GetMetaDataComponent from './getMetadataTask';
import ReviewTaskComponent from './reviewTask';

function MiddleContent({ activeTask }) {
  const taskComponents = {
    1: <UploadDataTaskComponent />,
    2: <ValidationTaskComponent />,
    3: <GetMetaDataComponent />,
    4: <ReviewTaskComponent />,
    // Add other task components here
  };

  return (
    <main>
      {taskComponents[activeTask]}
    </main>
  );
}

export default MiddleContent;
