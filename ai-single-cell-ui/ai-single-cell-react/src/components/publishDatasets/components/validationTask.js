import React, { useState, useEffect } from 'react';
import { FLASK_BACKEND_API, STORAGE } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Select from 'react-select';
import { PropagateLoader } from 'react-spinners';

function ValidationTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask }) {
  const [seuratFiles, setSeuratFiles] = useState([]);
  const [selectedSeuratFile, setSelectedSeuratFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [assayNamesMap, setAssayNamesMap] = useState({}); // Store fetched assay names
  const [addedFiles, setAddedFiles] = useState([]); // Track files that have been added


  const jwtToken = getCookie('jwtToken');
  const navigate = useNavigate();

  const fetchAssayNames = async (path, file) => {
    if (!assayNamesMap[file]) {
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
          setAssayNamesMap((prevMap) => ({
            ...prevMap,
            [file]: data.assay_names.map((name) => ({ label: name, value: name })),
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
              !seuratFiles.some((fileInfo) => fileInfo.label === file) && !addedFiles.includes(file)
            ) {
              setSeuratFiles((seuratFiles) => [
                ...seuratFiles,
                { label: file, value: path, assayNames: [], selectedAssays: [] },
              ]);

              // Add the file to the list of added files
            setAddedFiles((addedFiles) => [...addedFiles, { label: file, value: path }]);
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
  }, [jwtToken, taskData, navigate]);

  const handleTaskCompletion = () => {
    // Perform the necessary actions for completing Task 1
    // For example, submit a form, validate input, etc.

    // After Task 1 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      2: true, // Mark Task 1 as completed
    }));
    // The current task is finished, so make the next task active
    setActiveTask(3);
  };

  const handleSeuratFileChange = (target) => {
    
    const selectedValue = target.value;
    const selectedLabel = target.options[target.selectedIndex].text;

    const newOption = {"label": selectedLabel, "value": selectedValue};

    setSelectedSeuratFile(newOption);

    // Fetch assay names when a Seurat file is selected
    if (newOption) {
      fetchAssayNames(newOption.value, newOption.label);
    }
  };

  const handleAssayNamesChange = (selectedOptions) => {
    if (selectedSeuratFile) {
      const updatedSeuratFiles = [...seuratFiles];
      const index = updatedSeuratFiles.findIndex((file) => file.value === selectedSeuratFile.value);
      updatedSeuratFiles[index].selectedAssays = selectedOptions;
      setSeuratFiles(updatedSeuratFiles);
    }
  };

  return (
    <div>
        <div>
          <div>
            <h1>Select a Seurat File</h1>
            <select value={selectedSeuratFile.value} onChange={(e) => handleSeuratFileChange(e.target)}>
                <option value="">Select a Seurat File</option>
                {addedFiles.map((file, index) => (
                  <option key={index} value={file.value}>
                    {file.label}
                  </option>
                ))}
              </select>
              {loading ? (
                <div className="spinner-container">
                  <PropagateLoader color={'#36D7B7'} loading={loading} size={15} />
                </div>
              ) : (
                  <div>
                    {selectedSeuratFile && (
                      <>
                        <h1>Choose Assay Names</h1>
                        <Select
                          isMulti
                          options={assayNamesMap[selectedSeuratFile.label] || []}
                          value={selectedSeuratFile.selectedAssays}
                          onChange={handleAssayNamesChange}
                        />
                      </>
                    )}
              </div>
              )}
          </div>
          <div className="previous">
            <button type="submit" className="btn btn-info" onClick={() => setActiveTask(activeTask - 1)}>
              Previous
            </button>
          </div>
          <div className="next-upon-success">
            <button type="submit" className="btn btn-info" onClick={handleTaskCompletion}>
              Next
            </button>
          </div>
        </div>
    </div>
  );
}

export default ValidationTaskComponent;
