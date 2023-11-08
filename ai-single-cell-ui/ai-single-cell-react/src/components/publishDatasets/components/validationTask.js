import React, { useState, useEffect } from 'react';
import { FLASK_BACKEND_API,CELERY_BACKEND_API ,STORAGE } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import { PropagateLoader, RingLoader } from 'react-spinners';
import axios from 'axios';


function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask }) {
   const [loading, setLoading] = useState(false);
   const [validationLoading, setValidationLoading] = useState(false);
   const [errorMessage, setErrorMessage] = useState('');


  const jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();

  function extractFilename(inputPath) {
    // Use the `split` method to split the path by the forward slash ("/") and get the last part, which is the filename
    const parts = inputPath.split('/');
    const filename = parts[parts.length - 1];
    return filename;
  }

  // const fetchAssayNames = async (path, file) => {
  //   if (!taskData.validation.assayNamesMap[file]) {
  //     setLoading(true);

  //     try {
  //       const response = await fetch(`${CELERY_BACKEND_API}/convert/api/convert_to_anndata`, {
  //         method: 'POST',
  //         headers: {
  //           'Content-Type': 'application/json',
  //         },
  //         body: JSON.stringify({ path }),
  //       });

  //       if (response.ok) {
  //         const data = await response.json();


  //         if (data.assay_names.length === 0) {
  //           // Remove the file from seuratFiles
  //           const updatedSeuratFiles = taskData.validation.seuratFiles.filter((fileInfo) => fileInfo.label !== file);
  //           setTaskData((prevTaskData) => ({
  //             ...prevTaskData,
  //             validation: {
  //               ...prevTaskData.validation,
  //               fileMappings: [
  //                 ...prevTaskData.validation.fileMappings,
  //                 { fileDetails: path },
  //               ],
  //               seuratFiles: updatedSeuratFiles,
  //             },
  //           }));
  //         } else {   
  //           setTaskData((prevTaskData) => ({
  //             ...prevTaskData,
  //             validation: {
  //               ...prevTaskData.validation,
  //               assayNamesMap: {
  //                 ...prevTaskData.validation.assayNamesMap,
  //                 [file]: data.assay_names.map((name) => ({ label: name, value: name })),
  //               },
  //             },
  //           }));
  //         }
  //       } else {
  //         console.error('Error fetching assay names:', response.status);
  //       }
  //     } catch (error) {
  //       console.error('Error fetching assay names:', error);
  //     } finally {
  //       setLoading(false);
  //     }
  //   }
  // };


  useEffect(() => {
    isUserAuth(jwtToken)
    .then((authData) => {
      if(authData.isAdmin) {
          let username = authData.username;
          let newDirectoryPath = taskData.upload.newDirectoryPath;
          let files = taskData.upload.files;

          let inputFiles = [];
          for(let file of files) {
            let path = STORAGE + "/" + username + "/" + newDirectoryPath + "/" + file;
            inputFiles.push(path);
          }
          let hasAddedSeuratFiles = false; // Flag to track if seuratFiles have been added

          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            validation: {
              ...prevTaskData.validation,
              inputFiles: inputFiles
            },
          }));

          // Assuming you have already prepared the inputFiles array
          const requestData = {
            inputFiles: inputFiles.map(file => ({
              fileDetails: file
            })),
          };
          // Make the API call
          axios.post(`${CELERY_BACKEND_API}/convert/publishDatasets/validation`, requestData)
          .then(response => {

            // Handle the response from the server
            const results = response.data;
            // Handle the response from the server
            console.log('API Response:', response.data);

                // Iterate over the results array and process the data
            results.forEach(result => {
              console.log("Inside results");
              if (result.inputfile && (result.inputfile.endsWith('.h5Seurat') || result.inputfile.endsWith('.h5seurat') || result.inputfile.endsWith('.rds') || result.inputfile.endsWith('.Robj')) && !taskData.validation.seuratFiles.some((fileInfo) => fileInfo.value === result.inputfile)) {
                if(result.default_assay !== 'RNA') {
                  hasAddedSeuratFiles = true; // Update the flag
                  setTaskData((prevTaskData) => ({
                    ...prevTaskData,
                    validation: {
                      ...prevTaskData.validation,
                      seuratFiles: [
                        ...prevTaskData.validation.seuratFiles,
                        { label: extractFilename(result.inputfile), value: result.inputfile, assayNames: result.assay_names, selectedAssays: [] },
                      ],
                    },
                  }));
                  console.log("Inside results if");
                } else if(result.default_assay && result.default_assay === 'RNA') {
                  // Add the result directly to the qc_results array
                  setTaskData(prevTaskData => ({
                    ...prevTaskData,
                    quality_control: {
                      ...prevTaskData.quality_control,
                      qc_results: [
                        ...prevTaskData.quality_control.qc_results,
                        result, // Adding the entire result
                      ],
                    },
                  }));
                  console.log("Inside results else if");

                }
              } else {
                let fileDetails = {
                  fileDetails: result.inputfile,
                  format: result.format,
                  adata_path: result.adata_path
                };
                // Add the fileDetails directly to the fileMappings
                setTaskData((prevTaskData) => ({
                  ...prevTaskData,
                  validation: {
                    ...prevTaskData.validation,
                    fileMappings: [
                      ...prevTaskData.validation.fileMappings,
                      fileDetails,
                    ],
                  },
                }));
                console.log("Inside results else");

              }
            });

            console.log("done iterating");
            console.log(hasAddedSeuratFiles)

            if (!hasAddedSeuratFiles || taskData.validation.seuratFiles.length === 0) {
              console.log("done iterating if block");

                // validation step is successful, move to next task as there are no seurat or rds datasets
                setTaskData((prevTaskData) => ({
                  ...prevTaskData,
                  validation: {
                    ...prevTaskData.validation,
                    status: 'completed'
                  },
                }));

                setTaskStatus((prevTaskStatus) => ({
                  ...prevTaskStatus,
                  2: true, // Mark Task 2 as completed 
                }));
          
                // The current task is finished, so make the next task active
                setActiveTask(3); // Move to the next task (or update it to the appropriate task)
            }

            console.log("end of validation");

          })
          .catch(error => {
            console.error('API Error:', error.response);
            // console.error('Error Detail:', error.response.data.detail);
          });
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

          const hasSelectedAssays = taskData.validation.seuratFiles.every((file) => file.selectedAssays && file.selectedAssays.length > 0);

          if (!hasSelectedAssays) {
            setErrorMessage("Please select at least one assay for each Seurat file within available assays.");
            return;
          }

          // set Validation loading to true.
          setValidationLoading(true);

          console.log("Data to send")
          console.log(dataToSend);
        // Send the data to the backend API
        const response = await fetch(`${CELERY_BACKEND_API}/convert/publishDatasets/validation`, {
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
        // const newFileMappings = [];

        if (responseData && responseData.data && responseData.data.length > 0) {

          responseData.data.forEach((entry) => {
            // const fileDetails = entry.path; // Use the appropriate property

            // // Create an object to represent the file mapping
            // const fileMapping = {
            //   fileDetails: fileDetails,
            //   adata_path: entry.adata_path,
            //   assay: entry.assay,
            // };

            // // Add the file mapping to the list
            // newFileMappings.push(fileMapping);
            // Add the result directly to the qc_results array
            setTaskData(prevTaskData => ({
              ...prevTaskData,
              quality_control: {
                ...prevTaskData.quality_control,
                qc_results: [
                  ...prevTaskData.quality_control.qc_results,
                  entry, // Adding the entire result
                ],
              },
            }));

          });

          // Update the fileMappings state with the new list
          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            validation: {
              ...prevTaskData.validation,
              // fileMappings: newFileMappings,
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

    setErrorMessage('');
    
    setTaskData((prevTaskData) => ({
      ...prevTaskData,
      validation: {
        ...prevTaskData.validation,
        selectedSeuratFile: selectedOption,
      },
    }));

    // if (selectedOption) {
    //   fetchAssayNames(selectedOption.value, selectedOption.label);
    // }
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
                          options={taskData.validation.seuratFiles[taskData.validation.seuratFiles.findIndex((file) => file.value === taskData.validation.selectedSeuratFile.value)].assayNames}
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
