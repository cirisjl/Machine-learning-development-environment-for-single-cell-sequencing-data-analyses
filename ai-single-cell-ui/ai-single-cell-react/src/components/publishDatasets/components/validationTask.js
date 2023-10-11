import React, { useState, useEffect } from 'react';
import { CELERY_BACKEND_API, STORAGE } from '../../../constants/declarations';
import axios from 'axios';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';


function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask , activeTask}) {
  const [assayNames, setAssayNames] = useState([]);
  const [selectedAssays, setSelectedAssays] = useState([]);
  // const [adataPaths, setAdataPaths] = useState({});

  let jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();


   
  useEffect(() => {
    let username = ''

    isUserAuth(jwtToken)
    .then((authData) => {
      if(authData.isAdmin) {
          username = authData.username;
      } else {
        console.warn("Unauthorized - you must be an admin to access this page");
        navigate("/accessDenied");
      }

    })
    .catch((error) => {
      console.error(error);
    } 
    );
  
    let path = STORAGE + "/" + username + "/" + taskData.upload.newDirectoryPath + "/" + taskData.upload.files[0];
    const fetchAssayNames = async () => {
      try {
        const response = await axios.post(`${CELERY_BACKEND_API}/tools/convert_to_anndata`, {path}); // Replace with your endpoint and request data.
        setAssayNames(response.data.assay_names.map((name) => ({ label: name, value: name })));
      } catch (error) {
        console.error('Error fetching assay names:', error);
      }
    };
  
    fetchAssayNames();
  }, []);

  const handleSelectionChange = (selectedOptions) => {
    setSelectedAssays(selectedOptions.map((option) => option.value));
  };
  

  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      2: true, // Mark Task 1 as completed
    }));
    //The current task is finished, so make the next task active
    setActiveTask(3);
  };

  return (
    <div>
      <div>
        <h1>Choose Assay Names</h1>
        <Select
          isMulti
          options={assayNames}
          value={assayNames.filter((option) => selectedAssays.includes(option.value))}
          onChange={handleSelectionChange}
        />
        <button>Fetch adata_paths</button>
        {/* {Object.keys(adataPaths).map((name) => (
          <div key={name}>
            <p>{name}</p>
            <p>{adataPaths[name]}</p>
          </div>
        ))} */}
      </div>
      <div className='previous'>
        <button type="submit" class="btn btn-info" onClick={() => setActiveTask(activeTask - 1)} >Previous</button>
      </div>
      <div className='next-upon-success'>
        <button type="submit" class="btn btn-info" onClick={handleTaskCompletion} >Next</button>
      </div>
    </div>
  );
}

export default ValidationTaskComponent;
