import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { faFolderOpen } from "@fortawesome/free-solid-svg-icons";
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import './ModalWindow.css';
import { getCookie } from '../../utils/utilFunctions';
import Form from "@rjsf/core";
import close_icon from '../../assets/close_icon_u86.svg';
import close_icon_hover from '../../assets/close_icon_u86_mouseOver.svg';

import schema from "../../schema/react-json-schema/uploadDataSchema.json";
import RightRail from "../RightNavigation/rightRail";
import { useLocation } from 'react-router-dom';
import updateSchema from "./../updateDataSchema.json";
import FilePreviewModal from "./filePreviewModal";
import FileManagerModal from "./fileManagerModal";

export default function UploadData() {
    const [isFileManagerOpen, setIsFileManagerOpen] = useState(false);
    const [selectedFiles, setSelectedFiles] = useState([]);
    const [pwd, setPwd] = useState('/');
    const SERVER_URL = "http://" + process.env.REACT_APP_HOST_URL + ":3001";
    let jwtToken = getCookie('jwtToken');
    const [formData, setFormData] = useState({});
    const [isButtonDisabled, setIsButtonDisabled] = useState(false);
    const [enabledCheckboxes, setEnabledCheckboxes] = useState([]);
    const [tempFileList, setTempFileList] = useState([]);
    const [errorMessage, setErrorMessage] = useState('');
    const navigate = useNavigate();
    const [hoveredErrPopup, setHoveredErrPopup] = useState(false);
    const location = useLocation();
    let mode = location.state?.mode || '';
    let formInfo = location.state?.formInfo || '';
    const [currentFileList, setCurrentFileList] = useState(formInfo?.files?.map(file => file.file_loc) || []);
    const [previewBoxOpen, setPreviewBoxOpen] = useState(false);
    const [fileToPreview, setFileToPreview] = useState(null);
    const [fileNames, setFileNames] = useState([]);
    const [dirNames, setDirNames] = useState([]);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    let [selectedAliases, setSelectedAliases] = useState([]);
    const acceptedMultiFileNames = ['molecules.txt', 'annotation.txt', 'barcodes.tsv', 'genes.tsv', 'matrix.mtx', 'barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz', 'features.tsv', 'count_matrix.mtx', 'features.tsv.gz', 'count_matrix.mtx.gz'];
    const acceptedMultiFileSets = [
        ['molecules.txt', 'annotation.txt'],
        ['barcodes.tsv', 'genes.tsv', 'matrix.mtx'],
        ['barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz'],
        ['barcodes.tsv', 'features.tsv', 'count_matrix.mtx'],
        ['barcodes.tsv.gz', 'features.tsv.gz', 'count_matrix.mtx.gz']
    ];
    const handleMouseOver = () => {
        setHoveredErrPopup(true);
    };

    const handleMouseOut = () => {
        setHoveredErrPopup(false);
    };

    const handleCrossButtonClick = () => {
        setErrorMessage('');
    }

    useEffect(() => {
        const hookForUpdate = async () => {
            if (mode === 'update') {
                console.log('Mode: ' + mode);
                formInfo.reference = (formInfo.reference !== null) ? formInfo.reference : '';
                formInfo.summary = (formInfo.summary !== null) ? formInfo.summary : '';
                setSelectedFiles(currentFileList);
                setFormData(formInfo);
            }
        };
        hookForUpdate();
    }, []);


    useEffect(() => {
        const timeoutId = setTimeout(() => {
            setErrorMessage('');
        }, 30000);
        // Return a cleanup function to cancel the timeout when the component unmounts
        return () => clearTimeout(timeoutId);
    }, [errorMessage]);

    useEffect(() => {
        setSelectedAliases(selectedFiles.map(file => {
            const parts = file.split('/');
            return parts[parts.length - 1];
        }));
    }, [selectedFiles]);

    const SubmitButton = ({ disabled, onClick }) => {
        const handleClick = () => {
            setIsButtonDisabled(true);
            onClick();
        };

        return (
            <button
                type="submit"
                onClick={handleClick}
                disabled={isButtonDisabled || disabled}
                style={{ cursor: "pointer" }}
            >
                Submit
            </button>
        );
    };

    // toggle modal window visibility
    const toggleModal = async () => {
        await setIsFileManagerOpen(!isFileManagerOpen);

        if (!isFileManagerOpen) {
            setSelectedFiles([]);
            setTempFileList([]);
            fetchDirContents();
        }
    }

    function getAliasOptions(fileName) {
        if (fileName.endsWith('.txt')) {
            return ['molecules', 'annotation'];
        } else if (fileName.endsWith('.tsv')) {
            return ['genes', 'cells', 'features'];
        } else if (fileName.endsWith('.tsv.gz')) {
            return ['genes', 'cells', 'features'];
        } else if (fileName.endsWith('.mtx')) {
            return ['matrix', 'count matrix'];
        } else if (fileName.endsWith('.mtx.gz')) {
            return ['matrix', 'count matrix'];
        }
        else {
            return [];
        }
    };

    function getStandardFileName(fileName, fileType) {
        const acceptedFileTypes = ["molecules", "annotation", "cells", "genes", "matrix", "features", "count matrix"];
        if (!acceptedFileTypes.includes(fileType)) {
            return fileName;
        }
        const txt = { "molecules": "molecules.txt", "annotation": "annotation.txt" }
        const tsv = { "cells": "barcodes.tsv", "genes": "genes.tsv", "features": "features.tsv" }
        const tsv_gz = { "cells": "barcodes.tsv.gz", "genes": "genes.tsv.gz", "features": "features.tsv.gz" }
        const mtx = {"matrix": "matrix.mtx", "count matrix": "count_matrix.mtx"}
        const mtx_gz = {"matrix": "matrix.mtx.gz", "count matrix": "count_matrix.mtx.gz"}

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


    const fetchDirContents = async (subdir) => {

        enabledCheckboxes.forEach((checkboxId) => {
            try {

                const checkbox = document.getElementById(checkboxId);
                checkbox.checked = false;
            }
            catch (e) {
                console.log('Exception: ' + e);
            }
        });


        let newDir = ''
        jwtToken = getCookie('jwtToken');

        if (subdir == '..') {
            let slashIndex = pwd.lastIndexOf('/');
            newDir = slashIndex !== -1 ? pwd.substring(0, slashIndex) : pwd;
            setPwd(newDir)
        }
        else if (subdir) {
            newDir = pwd + '/' + subdir;
            setPwd(pwd + '/' + subdir);
        }
        else
            newDir = pwd

        fetch(`${SERVER_URL}/getDirContents?dirPath=${newDir}&authToken=${jwtToken}`)
            .then(response => {
                if (response.status === 403) {
                    throw new Error('Please log in first');
                }
                return response.json();
            })
            .catch(error => {
                setIsFileManagerOpen(false);
                if (error.message === 'Please log in first') {
                    navigate('/routing');
                } else {
                    console.error(error);
                }

            })
            .then(data => {
                setDirNames(data.Directories);
                setFileNames(data.Files);
            });

    }

    const handleSubmit = (event) => {
        if (selectedFiles === undefined || selectedFiles.length === 0) {
            setErrorMessage('Select at least one file.');
            return;
        }
        else if (jwtToken === undefined || jwtToken.length === 0) {
            setErrorMessage('Please log in first.');
            return;
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
                const file = selectedFiles[i];
                const lastSlashIndex = file.lastIndexOf('/');
                const fileName = file.substring(lastSlashIndex + 1, file.length);
                if (!acceptedMultiFileNames.includes(fileName)) {
                    if (lastSlashIndex !== -1) {
                        const prefix = file.substring(0, lastSlashIndex + 1);
                        selectedFiles[i] = prefix + selectedAliases[i];
                    } else {
                        selectedFiles[i] = selectedAliases[i]; // No slash found, push the original string
                    }
                    fetch(`${SERVER_URL}/renameFile?oldName=${file}&newName=${selectedFiles[i]}&authToken=${jwtToken}`, {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json'
                        },
                    })
                }
            }
        }
        else {
            const acceptedFormats = [".tsv", ".csv", ".txt.gz", ".txt", ".h5ad", "rds", "h5seurat", "tsv.gz", "mtx.gz", "h5", "xlsx", "hdf5", "gz"];
            if (!acceptedFormats.some(format => selectedFiles[0].endsWith(format))) {
                setErrorMessage("The selected file is not of an accepted standard format.");
                return;
            }
        }

        formData['authToken'] = jwtToken;
        formData['files'] = selectedFiles;

        if (mode === 'update') {
            formData.currentFileList = currentFileList;
            fetch(`${SERVER_URL}/updateDataset`, {
                method: 'PUT',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify(formData),
            })
                .then(response => {
                    if (response.status === 200) {
                        navigate('/dashboard', { state: { message: 'Dataset updated successfully.' } });
                    }
                    else {
                        console.log('Error updating dataset:', response);
                    }
                })
                .catch(error => {
                    setErrorMessage('Error updating dataset.');
                    console.error('Error updating dataset:', error);
                });
            return;
        }
        fetch(`${SERVER_URL}/createDataset`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(formData),
        })
            .then(response => {
                if (response.status === 201) {
                    navigate('/dashboard', { state: { message: 'Dataset created successfully.' } });
                }
                else if (response.status === 400) {
                    setErrorMessage(`Dataset '${formData.title}' already exists. Choose a different name.`);
                }
                else {
                    console.log('Error creating dataset:', response);
                }
            })
            .catch(error => {
                setErrorMessage('Error creating dataset.');
                console.error('Error creating dataset:', error);
            });
        setIsButtonDisabled(true);
        setTimeout(() => {
            setIsButtonDisabled(false);
        }, 5000);
    };

    if (jwtToken === '' || jwtToken === undefined) {
        navigate("/routing");
    }

    else return (
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                {(errorMessage !== '') && (
                    <div className='message-box' style={{ backgroundColor: 'lightpink', zIndex: 9999 }}>
                        <div style={{ textAlign: 'center' }}>
                            <p>{errorMessage}</p>
                            <div style={{ position: "absolute", right: "12px", top: "20px", cursor: "pointer" }}>
                                <img src={hoveredErrPopup ? close_icon_hover : close_icon} alt="close-icon" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut} onClick={handleCrossButtonClick} />
                            </div>
                        </div>
                    </div>)}
                <div>
                    <h2 style={{ textAlign: "left" }}><span>{(mode === 'update' ? 'Update Dataset' : "Create Dataset")}</span></h2>
                    {isInfoModalOpen && <div className="modal" style={{ zIndex: 9999, width: "30%", height: "40%" }}>
                        <div className='clear-icon'>
                            <img src={hoveredErrPopup ? close_icon_hover : close_icon} alt="close-icon" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut} onClick={() => setIsInfoModalOpen(false)} />
                        </div>
                        <div className="modal-content">
                            <div>
                                <p>
                                    Accepted Formats for Single-file Datasets: csv, tsv, txt, txt.gz, h5ad, rds, h5, hdf5
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
                    {previewBoxOpen && <FilePreviewModal selectedFile={fileToPreview} setPreviewBoxOpen={setPreviewBoxOpen} jwtToken={jwtToken} forResultFile={false} />}
                    {isFileManagerOpen && <FileManagerModal setFileToPreview={setFileToPreview} tempFileList={tempFileList} setEnabledCheckboxes={setEnabledCheckboxes} fileNames={fileNames} dirNames={dirNames} jwtToken={jwtToken} fetchDirContents={fetchDirContents} pwd={pwd} setPwd={setPwd} setPreviewBoxOpen={setPreviewBoxOpen} selectedFiles={selectedFiles} setSelectedFiles={setSelectedFiles} setErrorMessage={setErrorMessage} setTempFileList={setTempFileList} enabledCheckboxes={enabledCheckboxes} toggleModal={toggleModal} />}
                    <div>        <div>
                        <div id="upload-data-div">
                            <div className="info-icon" onClick={() => { setIsInfoModalOpen(true); }}>
                                <FontAwesomeIcon icon={faInfoCircle} size="1.2x" />
                            </div>
                            <b>Choose your files *</b> <br />
                            <div id="files-selected">
                                {selectedFiles.length > 1 ? (
                                    <>
                                        {selectedFiles.map((item, index) => {
                                            const itemName = item.substring(item.lastIndexOf('/') + 1);
                                            const showDropdown = !acceptedMultiFileNames.includes(itemName);

                                            return (
                                                <div key={index} onClick={() => console.log(item)}>
                                                    {item}&nbsp;&nbsp;
                                                    {showDropdown && (
                                                        <select
                                                            onChange={(e) => {
                                                                const updatedAliases = [...selectedAliases];
                                                                updatedAliases[index] = getStandardFileName(item, e.target.value);
                                                                setSelectedAliases(updatedAliases);
                                                            }}
                                                        >
                                                            <option selected>Set a standard file type</option>
                                                            {getAliasOptions(itemName).map((alias, aliasIndex) => (
                                                                <option key={aliasIndex}>{alias}</option>
                                                            ))}
                                                        </select>
                                                    )}
                                                </div>
                                            );
                                        })}


                                        <div style={{ color: 'red' }}>
                                            Notice: Files will be renamed to standard names of their corresponding type.
                                        </div>
                                    </>
                                ) : (
                                    <div>
                                        {selectedFiles[0]}
                                    </div>
                                )}
                            </div>

                            <button type="button" onClick={toggleModal} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
                                <FontAwesomeIcon icon={faFolderOpen} />
                            </button>

                        </div>
                        <br />
                        <h2 style={{ textAlign: "left" }}><span>Parameters</span></h2>


                        <Form
                            schema={(mode === 'update') ? updateSchema : schema}
                            formData={formData}
                            onChange={({ formData }) => setFormData(formData)}
                            onSubmit={handleSubmit}
                            SubmitButton={SubmitButton}
                        />

                    </div></div>
                </div >
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}