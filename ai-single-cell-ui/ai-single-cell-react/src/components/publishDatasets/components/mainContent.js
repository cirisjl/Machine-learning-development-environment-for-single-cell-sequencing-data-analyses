import React from 'react';
import UploadDataTaskComponent from './uploadDataTask';
import ValidationTaskComponent from './validationTask';
import GetMetaDataComponent from './getMetadataTask';
import ReviewTaskComponent from './reviewTask';
import QualityControlTaskComponent from './qualityControlTask';
import TaskBuilderTaskComponent from './taskBuilderTask';
import BenchmarksTaskComponent from './benchmarksTask';

function MiddleContent({ activeTask, setActiveTask, setTaskStatus, taskData, setTaskData, taskStatus}) {
  const taskComponents = {
    1: <UploadDataTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask} />,
    2: <ValidationTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask}/>,
    3: <QualityControlTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask}/>,
    4: <GetMetaDataComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask}/>,
    5: <TaskBuilderTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask} />,
    6: <BenchmarksTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask}/>,
    7: <ReviewTaskComponent setTaskStatus={setTaskStatus} taskData={taskData} setTaskData={setTaskData} setActiveTask = {setActiveTask} activeTask={activeTask}/>,
    // Add other task components here
  };

  return (
    <main>
      {taskComponents[activeTask]}
      
    </main>
  );
}

export default MiddleContent;
