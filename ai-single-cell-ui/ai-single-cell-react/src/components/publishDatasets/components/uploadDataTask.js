import React from 'react';
import UppyUploader from '../../MyData/uppy';
import { getCookie, isUserAuth, createUniqueFolderName, moveFilesToNewDirectory } from '../../../utils/utilFunctions';
import { SERVER_URL} from '../../../constants/declarations';
import axios from 'axios';
import { useState, useEffect } from 'react';
import close_icon from '../../../assets/close_icon_u86.svg';
import close_icon_hover from '../../../assets/close_icon_u86_mouseOver.svg';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { useNavigate } from 'react-router-dom';
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemText from '@mui/material/ListItemText';
import IconButton from '@mui/material/IconButton';
import DeleteIcon from '@mui/icons-material/Delete';
import { styled } from '@mui/material/styles';
import { Select, MenuItem } from '@mui/material';
import { FormControl, InputLabel } from '@mui/material';

function UploadDataTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask , activeTask}) {

  let pwd = "tempStorage/";
  const navigate = useNavigate();


  // State to manage error messages
  const [fileError, setFileError] = useState('');
  const [titleError, setTitleError] = useState('');
  const [errorMessage, setErrorMessage] = useState('');
  const [hoveredErrPopup, setHoveredErrPopup] = useState(false);

  const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
  const [selectedFiles, setSelectedFiles] = useState(taskData.upload.files);
  let [selectedAliases, setSelectedAliases] = useState(taskData.upload.files);
  const acceptedMultiFileNames = ['molecules.txt', 'annotation.txt', 'barcodes.tsv', 'genes.tsv', 'matrix.mtx', 'barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz', 'features.tsv', 'count_matrix.mtx', 'features.tsv.gz', 'count_matrix.mtx.gz'];
  const acceptedMultiFileSets = [
      ['molecules.txt', 'annotation.txt'],
      ['barcodes.tsv', 'genes.tsv', 'matrix.mtx'],
      ['barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz'],
      ['barcodes.tsv', 'features.tsv', 'count_matrix.mtx'],
      ['barcodes.tsv.gz', 'features.tsv.gz', 'count_matrix.mtx.gz']
  ];

  // Custom styled components
const ScrollableListContainer = styled('div')(({ theme }) => ({
  maxHeight: '400px', // Fixed height of the container
  overflowY: 'auto', // Enable vertical scrolling
  border: `1px solid ${theme.palette.divider}`, // Add border to distinguish the container
  borderRadius: theme.shape.borderRadius, // Use theme's border radius
  marginTop: theme.spacing(2),
}));

const CustomListItem = styled(ListItem)(({ theme }) => ({
  '&:hover': {
    backgroundColor: theme.palette.action.hover,
  },
  cursor: 'pointer', // Change cursor on hover to indicate an item is clickable
}));

  function getAliasOptions(fileName) {
    if (fileName.endsWith('.txt')) {
        return ['molecules', 'annotation'];
    } else if (fileName.endsWith('.tsv')) {
        return ['genes', 'cells', 'features'];
    } else if (fileName.endsWith('.tsv.gz')) {
        return ['genes', 'cells', 'features'];
    } else if (fileName.endsWith('.mtx')) {
        return ['matrix', 'count_matrix'];
    } else if (fileName.endsWith('.mtx.gz')) {
        return ['matrix', 'count_matrix'];
    }
    else {
        return [];
    }
};

function getStandardFileName(fileName, fileType) {
    const acceptedFileTypes = ["molecules", "annotation", "cells", "genes", "matrix", "features", "count_matrix"];
    if (!acceptedFileTypes.includes(fileType)) {
        return fileName;
    }
    const txt = { "molecules": "molecules.txt", "annotation": "annotation.txt" }
    const tsv = { "cells": "barcodes.tsv", "genes": "genes.tsv", "features": "features.tsv" }
    const tsv_gz = { "cells": "barcodes.tsv.gz", "genes": "genes.tsv.gz", "features": "features.tsv.gz" }
    const mtx = {"matrix": "matrix.mtx", "count_matrix": "count_matrix.mtx"}
    const mtx_gz = {"matrix": "matrix.mtx.gz", "count_matrix": "count_matrix.mtx.gz"}

    if (fileName.endsWith('.txt')) {
        return txt[fileType];
    } else if (fileName.endsWith('.tsv')) {
        return tsv[fileType];
    } else if (fileName.endsWith('.tsv.gz')) {
        return tsv_gz[fileType];
    } else if (fileName.endsWith('.mtx')) {
        return mtx[fileType];
    } else if (fileName.endsWith('.mtx.gz')) {
        return mtx_gz[fileType];
    }
};

  const handleMouseOver = () => {
    setHoveredErrPopup(true);
  };

  const handleMouseOut = () => {
    setHoveredErrPopup(false);
  };

  const handleCrossButtonClick = () => {
    setErrorMessage('');
  }
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

  useEffect(() => {
    setSelectedFiles(taskData.upload.files);
    setSelectedAliases(taskData.upload.files)
  }, [taskData.upload.files]); // Dependency array

  const removeFile = async (item, indexToRemove) => {
    try {
      // Send request to backend to delete the file
      await axios.delete(`${SERVER_URL}/api/storage/delete-file?fileName=${item}&authToken=${getCookie('jwtToken')}&newDirectoryPath=tempStorage`);

      // If successful, update the state to remove the file from the list
      setSelectedFiles(selectedFiles.filter((_, index) => index !== indexToRemove));
      setSelectedAliases(selectedAliases.filter((_, index) => index !== indexToRemove));
      const newTaskDataFiles = taskData.upload.files.filter((_, index) => index !== indexToRemove);
        // Update the taskData state
        setTaskData(previousTaskData => ({
            ...previousTaskData,
            upload: {
                ...previousTaskData.upload,
                files: newTaskDataFiles
            }
        }));


    } catch (error) {
      console.error("Error deleting file:", error);
      // Handle error (e.g., show error message to the user)
    }
};
  const handleTask1Completion = async () => {
    // Validate file upload and title input
    if (taskData.upload.files === undefined || taskData.upload.files.length === 0) {
      setFileError('Please upload a file.');
    } else {
      setFileError('');
    }
    if (selectedFiles.length > 1) {
        let isFileSelectionValid = false;
        acceptedMultiFileSets.forEach(function (multiFileSet) {
            for (let i = 0; i < multiFileSet.length; i++) {
                if (!selectedAliases.includes(multiFileSet[i])) {
                    break;
                }
                else if (i == multiFileSet.length - 1)
                    isFileSelectionValid = true;
            }
        });
        console.log('Selected Aliases: ' + selectedAliases);
        if (!isFileSelectionValid) {
            setErrorMessage("The set of selected files do not comply with the standard multi-file dataset requirements.");
            return;
        }
        for (let i = 0; i < selectedFiles.length; i++) {
            const fileName = selectedFiles[i];
            if (!acceptedMultiFileNames.includes(fileName)) {
                    selectedFiles[i] = selectedAliases[i];
            
                fetch(`${SERVER_URL}/api/storage/renameFile?oldName=tempStorage/${fileName}&newName=tempStorage/${selectedFiles[i]}`, {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json'
                    },
                })
                console.log("Make an API Call to rename the files that are stored");
            }
        }
        // Update the state of the task in the taskData state
        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          upload: {
            ...prevTaskData.upload,
            files: selectedAliases
          },
        }));
    }
    else if(selectedFiles.length === 1) {
        const acceptedFormats = [".tsv", ".csv", ".txt.gz", ".txt", ".h5ad", "rds", "h5seurat", "tsv.gz", "mtx.gz", "h5", "xlsx", "hdf5", "gz", "Robj", "zip", "rar", "tar", "tar.bz2", "tar.xz"];
        if (!acceptedFormats.some(format => selectedFiles[0].endsWith(format))) {
            setErrorMessage("The selected file is not of an accepted standard format.");
            return;
        }
    }
    if (!taskData.upload.title) {
      setTitleError('Please enter a title.');
    } else {
      setTitleError('');
    }

    // If both file and title are provided, continue to the next step
    if ((taskData.upload.files !== undefined && taskData.upload.files.length !== 0 ) && taskData.upload.title !== undefined) {
        
        try {
            if(taskData.upload.status !== 'completed') {  
          // Create a new directory with the title name
          const newDirectoryName = createUniqueFolderName(taskData.upload.title);

          let isMultiFileDataset = taskData.upload.files.length > 1 ? true : false;

          let newDirectoryPath = `projects/${newDirectoryName}`;
          // Move the uploaded files from tempStorage to the new directory
          await moveFilesToNewDirectory(newDirectoryPath); 

          // Update the state of the task in the taskData state
          setTaskData((prevTaskData) => ({
            ...prevTaskData,
            upload: {
              ...prevTaskData.upload,
              status: 'completed',
              newDirectoryPath: newDirectoryPath,
              isMultiFileDataset: isMultiFileDataset
            },
          }));

          setTaskStatus((prevTaskStatus) => ({
            ...prevTaskStatus,
            1: true, // Mark Task 1 as completed
          }));
        
      }

      //The current task is finished, so make the next task active
      setActiveTask(2);     
      
      } catch (error) {
        console.error('Error moving files:', error);
        // Handle the error as needed (e.g., display an error message to the user)
      }
    }
    console.log(taskData);
  };

  useEffect(() => {
    isUserAuth(getCookie('jwtToken'))
    .then((authData) => {
      if (authData.isAdmin) {
        console.log("User is admin and has access to this page");

      }  else {
        console.warn("Unauthorized - you must be an admin to access this page");
        navigate("/accessDenied");
      }
    })
    .catch((error) => {
      console.error(error);
    });
  }, []);


  return (
    <div className='upload-task'>
      {(errorMessage !== '') && (
        <div className='message-box' style={{ backgroundColor: 'lightpink', zIndex: 9999 }}>
            <div style={{ textAlign: 'center' }}>
                <p>{errorMessage}</p>
                <div style={{ position: "absolute", right: "12px", top: "20px", cursor: "pointer" }}>
                    <img src={hoveredErrPopup ? close_icon_hover : close_icon} alt="close-icon" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut} onClick={handleCrossButtonClick} />
                </div>
            </div>
      </div>)}
      <div className="separator heading">
          <div className="stripe"></div>
          <h2 className="h-sm font-weight-bold">Input</h2>
          <div className="stripe"></div>
      </div>
      <div className='uppy-uploader-component'>
        <div className="info-icon" onClick={() => { setIsInfoModalOpen(true); }}>
            <FontAwesomeIcon icon={faInfoCircle} size="1.2x" />
        </div>
        <span>Choose your file*</span>
        {isInfoModalOpen && <div className="modal" style={{ zIndex: 9999, width: "30%", height: "40%" }}>
            <div className='clear-icon'>
                <img src={hoveredErrPopup ? close_icon_hover : close_icon} alt="close-icon" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut} onClick={() => setIsInfoModalOpen(false)} />
            </div>
            <div className="modal-content">
                <div>
                    <p>
                        Accepted Formats for Single-file Datasets: csv, tsv, txt, txt.gz, h5ad, rds, h5, hdf5. h5seurat, Robj
                    </p>
                    <p>
                        Standard File Structure for Multi-file Datasets:
                    </p>
                    <ul>
                        <li>Molecules(txt)&nbsp;+&nbsp;Annotation(txt)</li>
                        <li>Cells(tsv)&nbsp;+&nbsp;Genes(tsv)&nbsp;+&nbsp;Matrix(mtx)</li>
                        <li>Cells(tsv.gz)&nbsp;+&nbsp;Genes(tsv.gz)&nbsp;+&nbsp;Matrix(mtx.gz)</li>
                        <li>Cells(tsv)&nbsp;+&nbsp;Features(tsv)&nbsp;+&nbsp;Count Matrix(mtx)</li>
                        <li>Cells(tsv.gz)&nbsp;+&nbsp;Features(tsv.gz)&nbsp;+&nbsp;Count Matrix(mtx.gz)</li>
                    </ul>
                </div>
            </div>
        </div>}
        <UppyUploader toPublishDataset={true} isUppyModalOpen={true} pwd={pwd} authToken={getCookie('jwtToken')} publicDatasetFlag= {true} setFileError={setFileError} setTaskData = {setTaskData}/>

        {selectedFiles.length > 0 && 
          <div id="files-selected">
            <ScrollableListContainer>
              <List dense>
                {selectedFiles.length > 1 ? (
                    <>
                      {selectedFiles.map((item, index) => {
                        const showDropdown = !acceptedMultiFileNames.includes(item);
                        return (
                          <div key={index} className="file-selections">
                            <CustomListItem key={index}>
                              <IconButton edge="start" aria-label="delete" onClick={() => removeFile(item, index)}>
                                  <DeleteIcon />
                              </IconButton>
                              {showDropdown && (
                                <FormControl sx={{ m: 1, minWidth: 120 }} size="small">
                                  <Select
                                    displayEmpty
                                    value={selectedAliases[index]}
                                    onChange={(e) => {
                                      const updatedAliases = [...selectedAliases];
                                      updatedAliases[index] = getStandardFileName(item, e.target.value);
                                      setSelectedAliases(updatedAliases);
                                    }}
                                    renderValue={(selected) => {
                                      if (selected && selected.length === 0) {
                                        return <em>Set a standard file type</em>;
                                      }
                                      return selected;
                                    }}
                                  >
                                    {getAliasOptions(item).map((alias, aliasIndex) => (
                                      <MenuItem key={aliasIndex} value={alias}>{alias}</MenuItem>
                                    ))}
                                  </Select>
                                </FormControl>
                              )}
                              <ListItemText primary={item} />
                            </CustomListItem>
                          </div>
                        );
                      })}
                    </>
                  ) : (
                    <CustomListItem>
                      <IconButton edge="start" aria-label="delete" onClick={() => removeFile(selectedFiles[0], 0)}>
                          <DeleteIcon />
                      </IconButton>                    
                      <ListItemText primary={selectedFiles[0]} />
                    </CustomListItem>
                )}
              </List>
            </ScrollableListContainer>
            {selectedFiles.length > 1 && 
              <div style={{ color: 'red' }}>
                Notice: Files will be renamed to standard names of their corresponding type.
              </div>
            }
          </div>
        }
        {fileError && <div className="error-message">{fileError}</div>}
      </div>
      <div className="separator heading">
          <div className="stripe"></div>
          <h2 className="h-sm font-weight-bold">Parameters</h2>
          <div className="stripe"></div>
      </div>
      <div className="form-group field field-string">
        <label className="control-label" for="root_title">Title<span className="required">*</span></label>
        <input className="form-control" id="root_title" label="Title" required="" placeholder="" type="text" value={taskData.upload.title} fdprocessedid="jwyrb9" onChange={handleTitleChange}></input>
        {titleError && <div className="error-message">{titleError}</div>}
      </div>
      <div className='next-upon-success'>
        <button type="submit" className="btn btn-info button" onClick={handleTask1Completion}>Next</button>
      </div>
    </div>
  );
}

export default UploadDataTaskComponent;
