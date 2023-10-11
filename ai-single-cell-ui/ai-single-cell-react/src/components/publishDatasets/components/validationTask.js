import React, { useState, useEffect } from 'react';
import { FLASK_BACKEND_API, STORAGE } from '../../../constants/declarations';
import axios from 'axios';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import debounce from 'lodash/debounce';
import { PropagateLoader } from 'react-spinners';


function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask , activeTask}) {
  const [assayNames, setAssayNames] = useState([]);
  const [selectedAssays, setSelectedAssays] = useState([]);
  const [seuratFiles, setSeuratFiles] = useState([]); // List to store Seurat files
  const [selectedSeuratFile, setSelectedSeuratFile] = useState(null); // Selected Seurat file

  // const [adataPaths, setAdataPaths] = useState({});
  const [isFetchingData, setIsFetchingData] = useState(false);
  const [loading, setLoading] = useState(false); // Initialize loading to true

  let jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();


  const fetchAssayNames = debounce(async (path) => {
    setIsFetchingData(true);
    setLoading(true);
    try {
      const assayNamesPromises = seuratFiles.map((fileInfo) => {
        return fetch(`${FLASK_BACKEND_API}/api/convert_to_anndata`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: fileInfo.value }),
        })
          .then((response) => response.json())
          .then((data) => {
            fileInfo.assayNames = data.assay_names.map((name) => ({ label: name, value: name }));
            fileInfo.selectedAssays = [];
            return fileInfo;
          });
      });
      

      const filesWithAssayNames = await Promise.all(assayNamesPromises);
      setSeuratFiles(filesWithAssayNames);
      setLoading(false);
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
          let newDirectoryPath = taskData.upload.newDirectoryPath;
          let files = taskData.upload.files;
          for (let file of files) {
            let path = STORAGE + "/" + username + "/" + newDirectoryPath + "/" + file;

            if (file.endsWith('.h5Seurat') || file.endsWith('.h5seurat') || file.endsWith('.rds')) {
              setSeuratFiles((seuratFiles) => [
                ...seuratFiles,
                { label: file, value: path, assayNames: [], selectedAssays: [] },
              ]);
            }
            
          }

          // Fetch assay names for Seurat files
          fetchAssayNames();
        } else {
          console.warn("Unauthorized - you must be an admin to access this page");
          navigate("/accessDenied");
        }
      })
      .catch((error) => {
        console.error(error);
      });
  }, []);


    // // When the user selects a Seurat file, fetch assay names for that file
    // useEffect(() => {
    //   if (selectedSeuratFile) {
    //     fetchAssayNames(selectedSeuratFile.value);
    //   }
    // }, [selectedSeuratFile]);
  

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

  // return (
  //   <div>
  //           {loading ? (
  //             // Render the loading spinner while loading is true
  //             <div className="spinner-container">
  //             <PropagateLoader color={'#36D7B7'} loading={loading} size={15} />
  //           </div>             
  //           ) : (
  //             <div>
  //               <div>
  //                 <h1>Select a Seurat File</h1>
  //                 <Select
  //                   options={seuratFiles}
  //                   value={selectedSeuratFile}
  //                   onChange={(selectedOption) => setSelectedSeuratFile(selectedOption)}
  //                 />
  //                 {selectedSeuratFile && (
  //                   <div>
  //                     <h1>Choose Assay Names</h1>
  //                     <Select
  //                       isMulti
  //                       options={assayNames}
  //                       value={assayNames.filter((option) => selectedAssays.includes(option.value))}
  //                       onChange={handleSelectionChange}
  //                     />
  //                   </div>
  //                 )}
  //               </div>
  //               <div className='previous'>
  //                 <button type="submit" class="btn btn-info" onClick={() => setActiveTask(activeTask - 1)} >Previous</button>
  //               </div>
  //               <div className='next-upon-success'>
  //                 <button type="submit" class="btn btn-info" onClick={handleTaskCompletion} >Next</button>
  //               </div>
  //             </div>
  //         )}
  //   </div>
  // );

  return (
    <div>
      {loading ? (
        <div className="spinner-container">
          <PropagateLoader color={'#36D7B7'} loading={loading} size={15} />
        </div>
      ) : (
        <div>
          <div>
            <h1>Select a Seurat File</h1>
            <Select
              options={seuratFiles}
              value={selectedSeuratFile}
              onChange={(selectedOption) => setSelectedSeuratFile(selectedOption)}
            />
            {selectedSeuratFile && (
              <>
                <h1>Choose Assay Names</h1>
                <Select
                  isMulti
                  options={selectedSeuratFile.assayNames}
                  value={selectedSeuratFile.selectedAssays}
                  onChange={(selectedOptions) => {
                    setSelectedSeuratFile((prevSeuratFile) => ({
                      ...prevSeuratFile,
                      selectedAssays: selectedOptions,
                    }));
                  }}
                />
              </>
            )}
          </div>
          <div className='previous'>
            <button type="submit" class="btn btn-info" onClick={() => setActiveTask(activeTask - 1)} >Previous</button>
          </div>
          <div className='next-upon-success'>
            <button type="submit" class="btn btn-info" onClick={handleTaskCompletion} >Next</button>
          </div>
        </div>
      )}
    </div>
  );
}

export default ValidationTaskComponent;
