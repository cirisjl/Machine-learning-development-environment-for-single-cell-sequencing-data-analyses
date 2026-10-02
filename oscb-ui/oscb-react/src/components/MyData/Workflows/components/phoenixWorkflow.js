import React, { useState, useEffect } from 'react';
import FilePreviewModal from "../../filePreviewModal";
import FileManagerModal from "../../fileManagerModal";
import Schema from '../../../../schema/react-json-schema/Workflows/phoenixWorkflow.json'
import { uiSchema } from '../../../../schema/UI-schema/Workflows/phoenixWorkflow';
import Form from '@rjsf/core';
import validator from '@rjsf/validator-ajv8';
import SelectComponent from '../../Tools/components/selectComponent';
import GeneRangeSlider from '../../Tools/components/geneRangeSlider';
import MultiSelectComponent from '../../Tools/components/multiselectComponent';
import RangeSlider from '../../Tools/components/sliderComponent';
import SwitchComponent from '../../Tools/components/switchComponent';
import UseDefaultSwitch from '../../Tools/components/useDefaultSwitch';
import ClusterLabelInput from '../../Tools/components/customInputComponent';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import AlertMessageComponent from '../../../publishDatasets/components/alertMessageComponent';
import { CELERY_BACKEND_API, NODE_API_URL } from '../../../../constants/declarations';
import { useNavigate } from 'react-router-dom';
import { isUserAuth, getCookie, copyFilesToPrivateStorage } from '../../../../utils/utilFunctions';
import Chatbot from "../../../RightNavigation/Chatbot";
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { faFolderOpen, faInfoCircle } from "@fortawesome/free-solid-svg-icons";
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemText from '@mui/material/ListItemText';
import IconButton from '@mui/material/IconButton';
import DeleteIcon from '@mui/icons-material/Delete';
import close_icon from '../../../../assets/close_icon_u86.svg';
import close_icon_hover from '../../../../assets/close_icon_u86_mouseOver.svg';
import { styled } from '@mui/material/styles';


export function PhoenixWorkFlowComponent(props) {
  const [selectedFiles, setSelectedFiles] = useState([]);
  const [selectedSraFiles, setSelectedSraFiles] = useState([]);
  const [isFileManagerOpen, setIsFileManagerOpen] = useState(false);
  const [fileManagerTarget, setFileManagerTarget] = useState('dataset');
  const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
  const [formErrors, setFormErrors] = useState("");
  const [formData, setFormData] = useState({});
  const [pwd, setPwd] = useState('/');
  let jwtToken = getCookie('jwtToken');
  const [enabledCheckboxes, setEnabledCheckboxes] = useState([]);
  const [tempFileList, setTempFileList] = useState([]);
  const [hoveredErrPopup, setHoveredErrPopup] = useState(false);
  const [previewBoxOpen, setPreviewBoxOpen] = useState(false);
  const [fileToPreview, setFileToPreview] = useState(null);
  const [fileNames, setFileNames] = useState([]);
  const [dirNames, setDirNames] = useState([]);
  const [publicDatasetFlag, setPublicDatasetFlag] = useState(false);
  const [isAdminuser, setIsAdminuser] = useState(false);

  const [message, setMessage] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
  const [isError, setIsError] = useState(false);
  const [presetQuestions, setPresetQuestions] = useState(null);

  let [selectedAliases, setSelectedAliases] = useState([]);
  let [selectedSraAliases, setSelectedSraAliases] = useState([]);

  const navigate = useNavigate();

   useEffect(() => {
      const containsPublicDatasets = pwd.startsWith('publicDatasets') || pwd.startsWith('/publicDatasets');
      setPublicDatasetFlag(containsPublicDatasets);
    }, [pwd]);

    useEffect(() => {
      isUserAuth(jwtToken)
        .then(async (authData) => {
            console.log(authData)
            if (authData.isAuth) {
                setIsAdminuser(authData.isAdmin);
                setPresetQuestions(getPresetQuestions());
            } else {
                console.warn("Unauthorized - pLease login first to continue");
                navigate("/routing");
            }
          })
          .catch((error) => console.error(error));
      // eslint-disable-next-line react-hooks/exhaustive-deps
        }, []);

  const extractDir = (inputFile) => {
    const fileLocParts = inputFile.split('/');
    fileLocParts.pop(); // Remove the file name from the array
    const output = fileLocParts.join('/'); // Join the remaining parts with '/'
    return output;
  };

  const filterCategoryMap = {
    phoenix: '/workflows/phoenix',
    // Add more filter categories and their corresponding URL paths as needed
  };

  // Custom styled components
  const ScrollableListContainer = styled('div')(({ theme }) => ({
      maxHeight: '600px', // Fixed height of the container
      overflowY: 'auto', // Enable vertical scrolling
      border: `1px solid ${theme.palette.divider}`, // Add border to distinguish the container
      borderRadius: theme.shape.borderRadius, // Use theme's border radius
      marginTop: theme.spacing(2),
  }));

  const removeFile = async (item, indexToRemove) => {
      try {

          // If successful, update the state to remove the file from the list
          setSelectedFiles(selectedFiles.filter((_, index) => index !== indexToRemove));
          setSelectedAliases(selectedAliases.filter((_, index) => index !== indexToRemove));

      } catch (error) {
          console.error("Error deleting file:", error);
          // Handle error (e.g., show error message to the user)
      }
  };

  const removeSraFile = async (item, indexToRemove) => {
    try {

      // If successful, update the state to remove the file from the list
      setSelectedSraFiles(selectedSraFiles.filter((_, index) => index !== indexToRemove));
      setSelectedSraAliases(selectedSraAliases.filter((_, index) => index !== indexToRemove));

    } catch (error) {
      console.error("Error deleting file:", error);
      // Handle error (e.g., show error message to the user)
    }
  };

  const CustomListItem = styled(ListItem)(({ theme }) => ({
      '&:hover': {
          backgroundColor: theme.palette.action.hover,
      },
      cursor: 'pointer', // Change cursor on hover to indicate an item is clickable
  }));

  const handleMouseOver = () => {
    setHoveredErrPopup(true);
  };

  const handleMouseOut = () => {
    setHoveredErrPopup(false);
  };

  // const handleCrossButtonClick = () => {
  //   setHasMessage(true);
  //   setMessage("");
  //   setIsError(true);
  // }

  useEffect(() => {
      setSelectedAliases((selectedFiles || []).map(file => {
          const parts = file.split('/');
          return parts[parts.length - 1];
      }));
      // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedFiles]);

  useEffect(() => {
    setSelectedSraAliases((selectedSraFiles || []).map(file => {
      const parts = file.split('/');
      return parts[parts.length - 1];
    }));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedSraFiles]);

  // FileManagerModal uses this callback to close itself after Select/Close.
  const toggleModal = async() => {
    await setIsFileManagerOpen(!isFileManagerOpen);

    if (!isFileManagerOpen) {
      setEnabledCheckboxes([]);
      setTempFileList([]);
      setPwd('/');
      fetchDirContents();
      }
    }

  // Open the shared file manager for one of the two selectors.
  const openFileManager = (target) => {
    setFileManagerTarget(target);
    // Start a fresh selection session. Existing selections remain unchanged
    // unless the user confirms a new list with the Select Files button.
    setTempFileList([]);
    setEnabledCheckboxes([]);
    setPwd('/');
    fetchDirContents();
    setIsFileManagerOpen(true);
  };

  // Adapter passed to FileManagerModal. It routes the selected files to the
  // selector that opened the modal.
  const setActiveSelectedFiles = (files) => {
      const nextFiles = Array.isArray(files) ? files : [];

      // FileManagerModal sends [] when its Close button is used. Preserve the
      // existing selector value in that case; files can still be removed from
      // the list with the delete buttons in this workflow.
      if (nextFiles.length === 0) {
          return;
      }

      if (fileManagerTarget === 'sra') {
          setSelectedSraFiles(nextFiles);
      } else {
          setSelectedFiles(nextFiles);
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

      if (subdir === "publicDatasets") {
          newDir = subdir;
          setPwd(newDir)
      } else if (pwd === "publicDatasets" && subdir === "..") {
          newDir = ""
          setPwd(newDir)
      }
      else if (subdir === '..') {
          let slashIndex = pwd.lastIndexOf('/');
          newDir = slashIndex !== -1 ? pwd.substring(0, slashIndex) : pwd;
          setPwd(newDir)
      }
      else if (subdir) {
          newDir = pwd + '/' + subdir;
          setPwd(pwd + '/' + subdir);
      }
      else {
          newDir = pwd
      }

      fetch(`${NODE_API_URL}/getDirContents?dirPath=${newDir}&authToken=${jwtToken}`)
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
            if (!data) return;
            setDirNames(data.Directories || []);
            setFileNames(data.Files || []);
        });
  }

  const widgets = {
    SelectComponent: SelectComponent,
    geneRangeSlider: GeneRangeSlider,
    MultiSelectComponent: MultiSelectComponent,
    toggle: (props) => (
      <Toggle
        checked={props.value}
        onChange={(e) => props.onChange(e.target.checked)}
      />
    ),
    GeneRangeSlider: GeneRangeSlider,
    RangeSlider: RangeSlider,
    SwitchComponent: SwitchComponent,
    ClusterLabelInput: ClusterLabelInput,
    UseDefaultSwitch: UseDefaultSwitch
  };

  const getPresetQuestions = () => {
    let presetQuestionsList = [];
    presetQuestionsList.push({
      title: "What is **PHoeNIx** (Portable Healthcare Nextgen Informatics pipeline)?", "prompt": "What is **PHoeNIx** (Portable Healthcare Nextgen Informatics pipeline)?" });
    return presetQuestionsList;
  };

  const isCSVFile = (path) => {
    return typeof path === 'string' &&
      path.toLowerCase().endsWith('.csv');
  };

  const isFolderPath = (path) => {
    if (typeof path !== 'string' || path.trim() === '') {
      return false;
    }

    // Explicit folder ending
    if (path.endsWith('/') || path.endsWith('\\')) {
      return true;
    }

    // Get the last part of the path
    const name = path.split(/[\\/]/).pop();

    // Consider it a folder if the last part has no extension
    return !name.includes('.');
  };

  const handleSubmit = ({ formData }) => {
    console.log('formData', formData);
    formData = formData.parameters;
    formData.method = formData['mode'];
    if (formData['description'] == null || formData['description'].trim() === '') {
      setHasMessage(true);
      setMessage("Please provide a job description before submitting the form.");
      setIsError(true);
      return;
    }
    if (!formData['mode']) {
      formData['mode'] = 'PHOENIX'; // Set a default mode if not provided
    }
    if (!formData['kraken2db']) {
      formData['kraken2db'] = 'Standard-8'; // Set a default Kraken2 database if not provided
    }
    
    // Verify the authenticity of the user
    isUserAuth(getCookie('jwtToken'))
      .then((authData) => {
        if (authData.isAuth) {
          // Prepare the data to send. For workflow we only use a single dataset for phoenix.
          //TODO - Change this logic accordingly based on other workflows in the future

          if (selectedFiles.length < 1) {
            setHasMessage(true);
            setMessage("Please select a Samplesheet Input or Directory Input of the project before submitting the form.");
            setIsError(true);
            setFormErrors("Please select a Samplesheet Input or Directory Input of the project before submitting the form.");
            return;
          } else if (selectedFiles.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE Samplesheet Input or Directory Input.");
            setIsError(true);
            return;
          } else if (!isFolderPath(selectedFiles[0]) && !isCSVFile(selectedFiles[0])) {
            setHasMessage(true);
            setMessage("Input must be either a CSV file or a directory.");
            setIsError(true);
            setFormErrors("Input must be either a CSV file or a directory.");
            return;
          } else {
            let filesFromPublic = false;
            
            // Logic to Copy files from public storage to user private storage if it is a public Dataset.
            for (const file of selectedFiles) {

                if (file.startsWith("publicDataset") || file.startsWith("/publicDatasets")) {
                    filesFromPublic = true;
                    break;
                }
            }

            if (filesFromPublic) {
              copyFilesToPrivateStorage(selectedFiles, authData.username)
                .then(response => {
                    if (response.success) {
                        console.log(response.message);
                    } else {
                        console.error(response.message);
                    }
                })
                .catch(error => {
                    console.error("Error:", error.message);
                });

              for (let i = 0; i < selectedFiles.length; i++) {
                selectedFiles[i] = selectedFiles[i].replace(/^\/?publicDatasets\//, '/');
              }
            }

            let updatedFiles = selectedFiles.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData.input = updatedFiles[0];
            if (isCSVFile(updatedFiles[0])) {
              const directory = extractDir(formData.input);
              formData.outdir = directory + "/phoenix";
            } else {
              formData['outdir'] = updatedFiles[0] + "/phoenix";
            }
          }
          formData['use_sra'] = formData['sra_params'].use_sra ? formData['sra_params'].use_sra : false;
          formData['scaffold_ext'] = formData['scaffolds_params'].scaffold_ext ? formData['scaffolds_params'].scaffold_ext : '.scaffolds.fa.gz';
          formData['centar'] = formData['scaffolds_params'].centar ? formData['scaffolds_params'].centar : false;
          formData['griphin_out'] = formData['scaffolds_params'].griphin_out ? formData['scaffolds_params'].griphin_out : '${launchDir}';
          formData['coverage'] = formData['optional_params'].coverage ? formData['optional_params'].coverage : 30;
          formData['create_ncbi_sheet'] = formData['optional_params'].create_ncbi_sheet ? formData['optional_params'].create_ncbi_sheet : false;
          formData['email'] = formData['optional_params'].email;
          formData['multiqc_title'] = formData['optional_params'].multiqc_title ? formData['optional_params'].multiqc_title : '';

          if (formData['sra_params'].use_sra === true && selectedSraFiles.length < 1) {
            setHasMessage(true);
            setMessage("Please select a csv file that has one SRR number per line.");
            setIsError(true);
            setFormErrors("Please select a csv file that has one SRR number per line.");
            return;
          } else if (formData['sra_params'].use_sra === true && selectedSraFiles.length === 1 && isCSVFile(selectedSraFiles[0]) ) {
            formData['use_sra'] = formData['sra_params'].use_sra;
            let updatedSraFiles = selectedSraFiles.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['input_sra'] = updatedSraFiles[0];
          } else if (formData['sra_params'].use_sra === true && selectedSraFiles.length > 1){
            setHasMessage(true);
            setMessage("Please select ONLY ONE CSV file that has one SRR number per line.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE CSV file that has one SRR number per line.");
            return;
          }
          // delete formData['sra_params'];
          // delete formData['scaffolds_params'];
          // delete formData['species_params'];
          // delete formData['optional_params'];
          setFormErrors("");

          formData.userID = authData.username;
          console.log("Adding user data");
          console.log(props.selectedWorkflow);
          console.log("In submit function");
          console.log("Submitted data:", formData);
          console.log(selectedFiles);

          const RELATIVE_PATH = filterCategoryMap[props.selectedWorkflow];
          console.log(CELERY_BACKEND_API + RELATIVE_PATH);
          fetch(CELERY_BACKEND_API + RELATIVE_PATH, {
            method: 'POST',
            headers: {
              'Content-Type': 'application/json',
            },
            body: JSON.stringify(formData),
          })
          .then(response => {
            // Check the status code
            if (response.ok) {
              return response.json();
            } else {
              throw new Error('Error while making a call to the celery Tools API');
            }
          })
          .then(response => {

            const jobId = response.job_id;
            console.log(jobId);

            setHasMessage(true);
            setMessage(response.status ? response.status : "Job Successfully Submitted.");
            setIsError(false);
            navigate("/mydata/wgsWorkflowTaskDetails", { state: { job_id: jobId, methodMap: formData.method, datasetURL: formData.input, description: formData.description, process: props.selectedWorkflow }});
          })
          .catch(error => {
            // Handle any errors that occur during the API call
            console.error("Form submission error:", error);
            setHasMessage(true);
            setMessage("Form submission error:", error);
            setIsError(true);
          });
        } else {
          console.warn("Unauthorized - please login first to continue");
          navigate("/routing");
        }
      })
      .catch((error) => {
        setHasMessage(true);
        setMessage("An error occurred while submitting the form.");
        setIsError(true);
      });
  };

  return (
    <div className='tools-container workflow-container common-class-tools-and-workflows'>

      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage={setMessage} isError={isError} />}
      {formErrors && <div style={{ color: 'red', marginBottom: '10px' }}>{formErrors}</div>}

      <div className="separator heading">
        <div className="stripe"></div>
        <h2 className="h-sm font-weight-bold">
          Input
        </h2>
        <div className="stripe"></div>
      </div>

      <div className={isFileManagerOpen ? "file-modal-overlay" : ""}>
          {isInfoModalOpen && <div className="modal" style={{ zIndex: 9999, width: "50%", height: "40%" }}>
              <div className='clear-icon'>
                  <img src={hoveredErrPopup ? close_icon_hover : close_icon} alt="close-icon" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut} onClick={() => setIsInfoModalOpen(false)} />
              </div>
              <div className="modal-content">
                  <div>
                      <p>
                        The full PHoeNIx pipeline (--mode PHOENIX or --mode CDC_PHOENIX) only runs on Illumina paired-end reads. Multiple samples can be run using a samplesheet.csv file.
                      </p>
                      <p>
                        The Samplesheet Input must be a comma-separated file (csv) with at least 3 columns and a header row, as shown in the example below. DO NOT HAVE ANY SPACES IN THIS FILE. An example samplesheet could be downloaded from the following link: <a href="https://github.com/CDCgov/mycosnp-nf/blob/b4216fec5d1ec2e7d8b136539748c3623a02fb45/assets/samplesheet.csv" target="_blank" rel="noopener noreferrer">Download Example Samplesheet</a>
                      </p>
                        <pre><code>
                          sample,fastq_1,fastq_2 <br/>
                          SAMPLE_1,$PATH/AEG588A1_S1_L002_R1_001.fastq.gz,$PATH/AEG588A1_S1_L002_R2_001.fastq.gz <br />
                          SAMPLE_2,$PATH/AEG588A2_S2_L002_R1_001.fastq.gz,$PATH/AEG588A2_S2_L002_R2_001.fastq.gz <br />
                          SAMPLE_3,$PATH/AEG588A3_S3_L002_R1_001.fastq.gz,$PATH/AEG588A3_S3_L002_R2_001.fastq.gz <br />
                        </code></pre>
                      <p>If you choose a Directory, the system will create a samplesheet from the directory of fastq files. The script will search 1 directory deep and attempt to determine sample id names and pairing/multilane information and will automatically create a samplesheet.</p>
                  </div>
              </div>
          </div>}
        {previewBoxOpen && (
          <div className="file-preview-layer">
            <FilePreviewModal
              selectedFile={fileToPreview}
              setPreviewBoxOpen={setPreviewBoxOpen}
              jwtToken={jwtToken}
              forResultFile={false}
            />
          </div>
        )}
        {isFileManagerOpen && (
          <FileManagerModal
            setFileToPreview={setFileToPreview}
            tempFileList={tempFileList}
            setEnabledCheckboxes={setEnabledCheckboxes}
            fileNames={fileNames}
            dirNames={dirNames}
            jwtToken={jwtToken}
            fetchDirContents={fetchDirContents}
            pwd={pwd}
            setPwd={setPwd}
            setPreviewBoxOpen={setPreviewBoxOpen}
            selectedFiles={fileManagerTarget === 'sra' ? selectedSraFiles : selectedFiles}
            setSelectedFiles={setActiveSelectedFiles}
            setFormErrors={setFormErrors}
            setTempFileList={setTempFileList}
            enabledCheckboxes={enabledCheckboxes}
            toggleModal={toggleModal}
            isAdminuser={isAdminuser}
            publicDatasetFlag={publicDatasetFlag}
          />
        )}
      </div>

      <div>
        <div id="upload-data-div">
          <div className="info-icon" onClick={() => { setIsInfoModalOpen(true); }}>
              <FontAwesomeIcon icon={faInfoCircle} size="1.2x" />
          </div>
          <b>Choose a Samplesheet CSV file or input Directory *</b> <br />
          {selectedFiles && selectedFiles.length > 0 &&
            <div id="files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedFiles.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile(item, index)}>
                                <DeleteIcon />
                            </IconButton>
                            <ListItemText primary={item} />
                          </CustomListItem>
                        </div>
                      );
                    })}
                  </>
                </List>
              </ScrollableListContainer>
              {selectedFiles && selectedFiles.length > 1 &&
                  <div style={{ color: 'red' }}>
                      Notice: Please select ONLY ONE CSV file or ONLY ONE Directory.
                  </div>
              }
            </div>
          }
          <div className="folder-button">
              <button type="button" onClick={() => openFileManager('input')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
                  <FontAwesomeIcon icon={faFolderOpen} />
              </button>
          </div>
        </div>
        <div>
          <p>
            The full PHoeNIx pipeline (--mode PHOENIX or --mode CDC_PHOENIX) only runs on Illumina paired-end reads. Multiple samples can be run using a samplesheet.csv file.
          </p>
          <p>
            The Samplesheet Input must be a comma-separated file (csv) with at least 3 columns and a header row, as shown in the example below. DO NOT HAVE ANY SPACES IN THIS FILE. An example samplesheet could be downloaded from the following link: <a href="https://github.com/CDCgov/mycosnp-nf/blob/b4216fec5d1ec2e7d8b136539748c3623a02fb45/assets/samplesheet.csv" target="_blank" rel="noopener noreferrer">Download Example Samplesheet</a>
          </p>
            <pre><code>
              sample,fastq_1,fastq_2 <br />
              SAMPLE_1,$PATH/AEG588A1_S1_L002_R1_001.fastq.gz,$PATH/AEG588A1_S1_L002_R2_001.fastq.gz <br />
              SAMPLE_2,$PATH/AEG588A2_S2_L002_R1_001.fastq.gz,$PATH/AEG588A2_S2_L002_R2_001.fastq.gz <br />
              SAMPLE_3,$PATH/AEG588A3_S3_L002_R1_001.fastq.gz,$PATH/AEG588A3_S3_L002_R2_001.fastq.gz <br />
            </code></pre>
          <p>If you choose a Directory, the system will create a samplesheet from the directory of fastq files. The script will search 1 directory deep and attempt to determine sample id names and pairing/multilane information and will automatically create a samplesheet.</p>
        </div>
      </div>

      <div className="separator heading">
        <div className="stripe"></div>
        <h2 className="h-sm font-weight-bold">
          Input SRA for --mode SRA and CDC_SRA (Optional)
        </h2>
        <div className="stripe"></div>
      </div>

      <div>
        <div id="upload-sra-data-div">
          <b>Choose a CSV file that has one SRR number per line</b> <br />
          {selectedSraFiles.length > 0 && (
            <div id="sra-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedSraFiles.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeSraFile(item, index)}>
                              <DeleteIcon />
                            </IconButton>
                            <ListItemText primary={item} />
                          </CustomListItem>
                        </div>
                      );
                    })}
                  </>
                </List>
              </ScrollableListContainer>
              {selectedSraFiles && selectedSraFiles.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE CSV file.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('sra')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          The samplesheet has strict requirement for there to be one column with each SRR number on each line. No header is necessary. DO NOT HAVE ANY SPACES IN THIS FILE. For example:
          <pre><code>
            SRR23709130<br />
            SRR23709131<br />
          </code></pre>
        </div>
      </div>

      <div className="form-component">
        <Form
          schema={Schema}
          uiSchema={uiSchema()}
          widgets={widgets}
          formData={formData}
          onChange={({ formData }) => setFormData(formData)}
          onSubmit={handleSubmit}
          onError={(errors) => console.log("Form errors:", errors)}
          validator={validator}
        />
      </div>
      <Chatbot presetQuestions={presetQuestions} />
    </div>
  )
};