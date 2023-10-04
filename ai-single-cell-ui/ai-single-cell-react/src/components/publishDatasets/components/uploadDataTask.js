import React from 'react';
import UppyUploader from '../../MyData/uppy';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';

function UploadDataTaskComponent({ setTaskStatus }) {

  let jwtToken = getCookie('jwtToken');

  let pwd = "tempStorage/";

  const handleTask1Completion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      1: true, // Mark Task 1 as completed
    }));
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
        <UppyUploader toPublishDataset={true} isUppyModalOpen={true} pwd={pwd} authToken={jwtToken} publicDatasetFlag= {false}/>
      </div>
      <div class="separator heading">
          <div class="stripe"></div>
          <h2 class="h-sm font-weight-bold">Parameters</h2>
          <div class="stripe"></div>
      </div>
      <div class="form-group field field-string">
        <label class="control-label" for="root_title">Title<span class="required">*</span></label>
        <input class="form-control" id="root_title" label="Title" required="" placeholder="" type="text" value="" fdprocessedid="jwyrb9"></input>
      </div>
      <div className='next-upon-success'>
        <button type="submit" class="btn btn-info" onClick={handleTask1Completion}>Next</button>
      </div>
    </div>
  );
}

export default UploadDataTaskComponent;
