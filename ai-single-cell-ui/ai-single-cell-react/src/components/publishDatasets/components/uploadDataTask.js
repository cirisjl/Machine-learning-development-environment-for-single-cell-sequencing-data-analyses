import React from 'react';
import UppyUploader from '../../MyData/uppy';
import { getCookie, isUserAuth, createUniqueFolderName, moveFilesToNewDirectory } from '../../../utils/utilFunctions';
import { useState, useEffect } from 'react';

function UploadDataTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask , activeTask}) {

  let jwtToken = getCookie('jwtToken');

  let pwd = "tempStorage/";

  // State to manage error messages
  const [fileError, setFileError] = useState('');
  const [titleError, setTitleError] = useState('');

  // Handle the title input change
  const handleTitleChange = (e) => {
    const newTitle = e.target.value;
    setTitleError('');

    // Update the title in the taskData state
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      upload: {
        ...prevTaskData.upload,
        title: newTitle,
      },
    }));
  };

  const handleTask1Completion = () => {
    // Validate file upload and title input
    if (taskData.upload.files === undefined || taskData.upload.files.length === 0) {
      setFileError('Please upload a file.');
    } else {
      setFileError('');
    }

    if (!taskData.upload.title) {
      setTitleError('Please enter a title.');
    } else {
      setTitleError('');
    }

    // If both file and title are provided, continue to the next step
    if ((taskData.upload.files !== undefined && taskData.upload.files.length !== 0 ) && taskData.upload.title !== undefined) {

        //Move the uploaded files from tempStorage to the admin user project folder
       
        // Create a new directory with the title name
        const newDirectoryName = createUniqueFolderName(taskData.upload.title);
        const newDirectoryPath = `projects/${newDirectoryName}`;

        // Move the uploaded files from tempStorage to the new directory
        moveFilesToNewDirectory(newDirectoryPath); 

        setTaskStatus((prevTaskStatus) => ({
          ...prevTaskStatus,
          1: true, // Mark Task 1 as completed
        }));
      
        // Update the state of the task in the taskData state
      setTaskData((prevTaskData) => ({
        ...prevTaskData,
        upload: {
          ...prevTaskData.upload,
          status: 'completed',
        },
      }));

      //The current task is finished, so make the next task active
      setActiveTask(2);      
    }
    console.log(taskData);
  };

  return (
    <div className='upload-task'>
      <div class="separator heading">
          <div class="stripe"></div>
          <h2 class="h-sm font-weight-bold">Input</h2>
          <div class="stripe"></div>
      </div>
      <div className='uppy-uploader-component'>
        <span>Choose your file*</span>
        <UppyUploader toPublishDataset={true} isUppyModalOpen={true} pwd={pwd} authToken={jwtToken} publicDatasetFlag= {true} setFileError={setFileError} setTaskData = {setTaskData}/>
        {taskData.upload.files && 
          <div className="uploaded-files">
            <h3 className="file-list-heading">Uploaded Files:</h3>
            <ul className="file-list">
              {taskData.upload.files.map((filename, index) => (
                <li key={index} className="file-list-item">{filename}</li>
              ))}
            </ul>
          </div>
        }
        {fileError && <div className="error-message">{fileError}</div>}
      </div>
      <div class="separator heading">
          <div class="stripe"></div>
          <h2 class="h-sm font-weight-bold">Parameters</h2>
          <div class="stripe"></div>
      </div>
      <div class="form-group field field-string">
        <label class="control-label" for="root_title">Title<span class="required">*</span></label>
        <input class="form-control" id="root_title" label="Title" required="" placeholder="" type="text" value={taskData.upload.title} fdprocessedid="jwyrb9" onChange={handleTitleChange}></input>
        {titleError && <div className="error-message">{titleError}</div>}
      </div>
      <div className='next-upon-success'>
        <button type="submit" class="btn btn-info" onClick={handleTask1Completion}>Next</button>
      </div>
    </div>
  );
}

export default UploadDataTaskComponent;
