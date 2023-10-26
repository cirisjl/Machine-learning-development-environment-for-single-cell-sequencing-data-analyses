import React, { useState, useEffect } from 'react';
import { FLASK_BACKEND_API,CELERY_BACKEND_API ,STORAGE } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import { PropagateLoader, RingLoader } from 'react-spinners';

function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask }) {
   const [loading, setLoading] = useState(false);
   const [validationLoading, setValidationLoading] = useState(false);
   const [errorMessage, setErrorMessage] = useState('');


  const jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();

  const fetchAssayNames = async (path, file) => {
    if (!taskData.validation.assayNamesMap[file]) {
      setLoading(true);

      try {
        const response = await fetch(`${CELERY_BACKEND_API}/convert/api/convert_to_anndata`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path }),
        });

        if (response.ok) {
          const data = await response.json();
          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            validation: {
              ...prevTaskData.validation,
              assayNamesMap: {
                ...prevTaskData.validation.assayNamesMap,
                [file]: data.assay_names.map((name) => ({ label: name, value: name })),
              },
            },
          }));
        } else {
          console.error('Error fetching assay names:', response.status);
        }
      } catch (error) {
        console.error('Error fetching assay names:', error);
      } finally {
        setLoading(false);
      }
    }
  };

  useEffect(() => {
    isUserAuth(jwtToken)
      .then((authData) => {
        if (authData.isAdmin) {
          let username = authData.username;
          let newDirectoryPath = taskData.upload.newDirectoryPath;
          let files = taskData.upload.files;

          for (let file of files) {
            let path = STORAGE + "/" + username + "/" + newDirectoryPath + "/" + file;
            if (
              (file.endsWith('.h5Seurat') || file.endsWith('.h5seurat') || file.endsWith('.rds')) &&
              !taskData.validation.seuratFiles.some((fileInfo) => fileInfo.label === file)
            ) {
              setTaskData((prevTaskData) => ({
                ...prevTaskData,
                validation: {
                  ...prevTaskData.validation,
                  seuratFiles: [
                    ...prevTaskData.validation.seuratFiles,
                    { label: file, value: path, assayNames: [], selectedAssays: [] },
                  ],
                },
              }));
            }
          }
        } else {
          console.warn("Unauthorized - you must be an admin to access this page");
          navigate("/accessDenied");
        }
      })
      .catch((error) => {
        console.error(error);
      });
  }, [jwtToken]);


  useEffect(() => {
    console.log(taskData);
  }, [taskData]);

  const handleTaskCompletion = async () => {
    try {
        if(taskData.validation.status !== 'completed') {

          // Prepare the data to send to the backend
          const dataToSend = [];

          const hasSelectedAssays = taskData.validation.seuratFiles.every((file) => file.selectedAssays && file.selectedAssays.length > 0);

          if (!hasSelectedAssays) {
            setErrorMessage("Please select at least one assay for each Seurat file within available assays.");
            return;
          }
          
          
          console.log(" hasSelectedAssays");
          console.log(hasSelectedAssays);
          console.log(errorMessage);

          // set Validation loading to true.
          setValidationLoading(true);

          console.log("Data to send")
          console.log(dataToSend);
        // Send the data to the backend API
        const response = await fetch(`${CELERY_BACKEND_API}/convert/api/convert_sce_to_annData`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify(dataToSend),
        });
    
        if (response.ok) {
        // Handle the response from the API
        const responseData = await response.json();

        // Initialize a list to store file mappings
        const newFileMappings = [];

        if (responseData && responseData.data && responseData.data.length > 0) {

          responseData.data.forEach((entry) => {
            const fileDetails = entry.path; // Use the appropriate property

            // Create an object to represent the file mapping
            const fileMapping = {
              fileDetails: fileDetails,
              adata_path: entry.adata_path,
              assay: entry.assay,
            };

            // Add the file mapping to the list
            newFileMappings.push(fileMapping);
          });

          // Update the fileMappings state with the new list
          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            validation: {
              ...prevTaskData.validation,
              fileMappings: newFileMappings,
              status: 'completed'
            },
          }));
          // After the API call is complete, you can update the task status
          setTaskStatus((prevTaskStatus) => ({
            ...prevTaskStatus,
            2: true, // Mark Task 3 as completed (or update it to the appropriate task)
          }));
    
          // The current task is finished, so make the next task active
          setActiveTask(3); // Move to the next task (or update it to the appropriate task)
          setValidationLoading(false);
        } else {
          setValidationLoading(false);
          console.error('Error making API call:', response.status);
        }
      }
    } else {
      setActiveTask(3); // Move to the next task (or update it to the appropriate task)
    }
  } catch (error) {
      console.error('Error making API call:', error);
    }
  };
  

  const handleSeuratFileChange = (selectedOption) => {
    
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      validation: {
        ...prevTaskData.validation,
        selectedSeuratFile: selectedOption,
      },
    }));

    if (selectedOption) {
      fetchAssayNames(selectedOption.value, selectedOption.label);
    }
  };

  const handleAssayNamesChange = (selectedOptions) => {
    if (taskData.validation.selectedSeuratFile) {
      const updatedSeuratFiles = [...taskData.validation.seuratFiles];
      const index = updatedSeuratFiles.findIndex((file) => file.value === taskData.validation.selectedSeuratFile.value);
      updatedSeuratFiles[index].selectedAssays = selectedOptions;
      setTaskData((prevTaskData) => ({
        ...prevTaskData,
        validation: {
          ...prevTaskData.validation,
          seuratFiles: updatedSeuratFiles,
        },
      }));    
    setErrorMessage('');
    } 
  };

  return (
    <div className='validation-task'>
      { validationLoading ? (
        <div className="spinner-container">
          <RingLoader color={'#36D7B7'} loading={validationLoading} size={150} />
        </div>
      ) : (
        <div className='container'>
          <div>
            <h1 className="header">Select a Seurat File</h1>
            <Select
            options={taskData.validation.seuratFiles.map((fileInfo) => ({ label: fileInfo.label, value: fileInfo.value }))} value={taskData.validation.selectedSeuratFile}
            onChange={handleSeuratFileChange}
            />
              {loading ? (
                <div className="spinner-container">
                  <PropagateLoader color={'#36D7B7'} loading={loading} size={15} />
                </div>
              ) : (
                  <div>
                    {taskData.validation.selectedSeuratFile && (
                      <>
                        <h1 className="header">Choose Assay Names</h1>
                        <Select
                          isMulti
                          options={taskData.validation.assayNamesMap[taskData.validation.selectedSeuratFile.label] || []}
                          value={taskData.validation.seuratFiles[taskData.validation.seuratFiles.findIndex((file) => file.value === taskData.validation.selectedSeuratFile.value)].selectedAssays}
                          onChange={handleAssayNamesChange}
                        />
                      </>
                    )}

              </div>
              )}
          </div>

          {errorMessage && <div className="error-message">{errorMessage}</div>}

          <div className='navigation-buttons'>
            <div className="previous">
              <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                Previous
              </button>
            </div>
            <div className="next-upon-success">
              <button type="submit" className="btn btn-info button" onClick={handleTaskCompletion}>
                Next
              </button>
            </div>
          </div>
        </div>
        )}
    </div>
  );
}

export default ValidationTaskComponent;
