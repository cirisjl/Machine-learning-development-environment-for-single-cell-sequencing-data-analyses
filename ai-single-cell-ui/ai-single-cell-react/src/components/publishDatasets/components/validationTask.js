import React, { useState, useEffect } from 'react';
import { FLASK_BACKEND_API, STORAGE } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import { PropagateLoader } from 'react-spinners';

function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask }) {
   const [loading, setLoading] = useState(false);

  const jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();

  const fetchAssayNames = async (path, file) => {
    if (!taskData.validation.assayNamesMap[file]) {
      setLoading(true);

      try {
        const response = await fetch(`${FLASK_BACKEND_API}/api/convert_to_anndata`, {
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
    // Prepare the data to send to the backend
    const dataToSend = [];
  
    taskData.validation.seuratFiles.forEach((file) => {
      // Check if any assays are selected for this file
      if (file.selectedAssays && file.selectedAssays.length > 0) {
        file.selectedAssays.forEach((assay) => {
          // Create an entry with the complete file details and assay name
          dataToSend.push({
            fileDetails: file.value,
            assayName: assay.value,
          });
        });
      }
    });
  
    try {
      // Send the data to the backend API
      const response = await fetch(`${FLASK_BACKEND_API}/api/convert_sce_to_annData`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(dataToSend),
      });
  
      if (response.ok) {
        // Handle the response from the API
        const responseData = await response.json();
        // Optionally, you can perform actions based on the response
        console.log('API response:', responseData);
  
        // After the API call is complete, you can update the task status
        setTaskStatus((prevTaskStatus) => ({
          ...prevTaskStatus,
          2: true, // Mark Task 3 as completed (or update it to the appropriate task)
        }));
  
        // The current task is finished, so make the next task active
        setActiveTask(3); // Move to the next task (or update it to the appropriate task)
      } else {
        console.error('Error making API call:', response.status);
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
      }));    } 
  };

  return (
    <div className='validation-task'>
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
    </div>
  );
}

export default ValidationTaskComponent;
