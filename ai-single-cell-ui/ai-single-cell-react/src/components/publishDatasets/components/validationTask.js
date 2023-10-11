import React, { useState, useEffect } from 'react';
import { CELERY_BACKEND_API, STORAGE } from '../../../constants/declarations';
import axios from 'axios';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import debounce from 'lodash/debounce';


function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask , activeTask}) {
  const [assayNames, setAssayNames] = useState([]);
  const [selectedAssays, setSelectedAssays] = useState([]);
  // const [adataPaths, setAdataPaths] = useState({});
  const [isFetchingData, setIsFetchingData] = useState(false);

  let jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();


  const fetchAssayNames = debounce(async (path) => {
    setIsFetchingData(true);
    try {
      const response = await fetch(`${CELERY_BACKEND_API}/tools/convert_to_anndata`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ path }),
      });

      if (response.ok) {
        const data = await response.json();
        setAssayNames(data.assay_names.map((name) => ({ label: name, value: name })));
      } else {
        console.error('Error fetching assay names:', response.status);
      }
    } catch (error) {
      console.error('Error fetching assay names:', error);
    } finally {
      setIsFetchingData(false);
    }
  }, 300); // Debounce for 300 milliseconds

  useEffect(() => {
    isUserAuth(jwtToken)
      .then((authData) => {
        if (authData.isAdmin) {
          let username = authData.username;
          let path = STORAGE + "/" + username + "/" + taskData.upload.newDirectoryPath + "/" + taskData.upload.files[0];

          fetchAssayNames(path);
        } else {
          console.warn("Unauthorized - you must be an admin to access this page");
          navigate("/accessDenied");
        }
      })
      .catch((error) => {
        console.error(error);
      });
  }, []);

  // useEffect(() => {
  //   isUserAuth(jwtToken)
  //     .then((authData) => {
  //       if (authData.isAdmin) {
  //         let username = authData.username;
  //         let path = STORAGE + "/" + username + "/" + taskData.upload.newDirectoryPath + "/" + taskData.upload.files[0];
  //         const fetchAssayNames = () => {
  //           try {
  //             fetch(`${CELERY_BACKEND_API}/tools/convert_to_anndata`, {
  //               method: 'POST',
  //               headers: {
  //                 'Content-Type': 'application/json',
  //               },
  //               body: JSON.stringify({ path }),
  //             })
  //               .then((response) => {
  //                 if (response.ok) {
  //                   response.json()
  //                     .then((data) => {
  //                       setAssayNames(data.assay_names.map((name) => ({ label: name, value: name })));
  //                     })
  //                     .catch((error) => {
  //                       console.error('Error parsing JSON response:', error);
  //                     });
  //                 } else {
  //                   console.error('Error fetching assay names:', response.status);
  //                 }
  //               })
  //               .catch((error) => {
  //                 console.error('Error making the request:', error);
  //               });
  //           } catch (error) {
  //             console.error('Error fetching assay names:', error);
  //           }
  //         };
  
  //         fetchAssayNames();
  //       } else {
  //         console.warn("Unauthorized - you must be an admin to access this page");
  //         navigate("/accessDenied");
  //       }
  //     })
  //     .catch((error) => {
  //       console.error(error);
  //     });
  // }, []);
  
  
  

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
