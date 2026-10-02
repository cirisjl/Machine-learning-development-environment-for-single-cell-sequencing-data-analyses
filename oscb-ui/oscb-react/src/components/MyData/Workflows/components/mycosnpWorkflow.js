import React, { useState, useEffect } from 'react';
import FilePreviewModal from "../../filePreviewModal";
import FileManagerModal from "../../fileManagerModal";
import Schema from '../../../../schema/react-json-schema/Workflows/mycosnpWorkflow.json'
import { uiSchema } from '../../../../schema/UI-schema/Workflows/mycosnpWorkflow';
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


export function MycoSNPWorkFlowComponent(props) {
  const [selectedFiles, setSelectedFiles] = useState([]);
  const [selectedSraFiles, setSelectedSraFiles] = useState([]);
  const [selectedVcfFiles, setSelectedVcfFiles] = useState([]);
  // Pre - MycoSNP workflow options
  const [selectedGambitDb, setSelectedGambitDb] = useState([]);
  const [selectedGambitH5Dir, setSelectedGambitH5Dir] = useState([]);
  const [selectedSubtypeDb, setSelectedSubtypeDb] = useState([]);
  // Main MycoSNP workflow
  const [selectedFastaFiles, setSelectedFastaFiles] = useState([]);
  const [selectedRefDir, setSelectedRefDir] = useState([]);
  const [selectedRefMaskedFasta, setSelectedRefMaskedFasta] = useState([]);
  const [selectedRefFai, setSelectedRefFai] = useState([]);
  const [selectedRefDict, setSelectedRefDict] = useState([]);
  const [selectedSnpeffconfig, setSelectedSnpeffconfig] = useState([]);
  const [selectedRefBwa, setSelectedRefBwa] = useState([]);
  const [selectedSkipSamplesFile, setSelectedSkipSamplesFile] = useState([]);
  const [activeOpenFiles, setActiveOpenFiles] = useState([]);

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
  let [selectedVcfAliases, setSelectedVcfAliases] = useState([]);
  // Pre - MycoSNP workflow options
  let [selectedGambitDbAliases, setSelectedGambitDbAliases] = useState([]);
  let [selectedGambitH5DirAliases, setSelectedGambitH5DirAliases] = useState([]);
  let [selectedSubtypeDbAliases, setSelectedSubtypeDbAliases] = useState([]);
  // Main MycoSNP workflow
  let [selectedFastaAliases, setSelectedFastaAliases] = useState([]);
  let [selectedRefDirAliases, setSelectedRefDirAliases] = useState([]);
  let [selectedRefMaskedFastaAliases, setSelectedRefMaskedFastaAliases] = useState([]);
  let [selectedRefFaiAliases, setSelectedRefFaiAliases] = useState([]);
  let [selectedRefDictAliases, setSelectedRefDictAliases] = useState([]);
  let [selectedSnpeffconfigAliases, setSelectedSnpeffconfigAliases] = useState([]);
  let [selectedRefBwaAliases, setSelectedRefBwaAliases] = useState([]);
  let [selectedSkipSamplesFileAliases, setSelectedSkipSamplesFileAliases] = useState([]);

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

  // const removeFile = async (item, indexToRemove) => {
  //     try {

  //         // If successful, update the state to remove the file from the list
  //         setSelectedFiles(selectedFiles.filter((_, index) => index !== indexToRemove));
  //         setSelectedAliases(selectedAliases.filter((_, index) => index !== indexToRemove));

  //     } catch (error) {
  //         console.error("Error deleting file:", error);
  //         // Handle error (e.g., show error message to the user)
  //     }
  // };

  const removeFile = async (target, item, indexToRemove) => {
    try {
      if (target === 'add_sra_file') {
        // If successful, update the state to remove the file from the list
        setSelectedSraFiles(selectedSraFiles.filter((_, index) => index !== indexToRemove));
        setSelectedSraAliases(selectedSraAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'add_vcf_file') {
        // If successful, update the state to remove the file from the list
        setSelectedVcfFiles(selectedVcfFiles.filter((_, index) => index !== indexToRemove));
        setSelectedVcfAliases(selectedVcfAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'gambit_db') {
        // If successful, update the state to remove the file from the list
        setSelectedGambitDb(selectedGambitDb.filter((_, index) => index !== indexToRemove));
        setSelectedGambitDbAliases(selectedGambitDbAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'gambit_h5_dir') {
        // If successful, update the state to remove the file from the list
        setSelectedGambitH5Dir(selectedGambitH5Dir.filter((_, index) => index !== indexToRemove));
        setSelectedGambitH5DirAliases(selectedGambitH5DirAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'subtype_db') {
        // If successful, update the state to remove the file from the list
        setSelectedSubtypeDb(selectedSubtypeDb.filter((_, index) => index !== indexToRemove));
        setSelectedSubtypeDbAliases(selectedSubtypeDbAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'fasta') {
        // If successful, update the state to remove the file from the list
        setSelectedFastaFiles(selectedFastaFiles.filter((_, index) => index !== indexToRemove));
        setSelectedFastaAliases(selectedFastaAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'ref_dir') {
        // If successful, update the state to remove the file from the list
        setSelectedRefDir(selectedRefDir.filter((_, index) => index !== indexToRemove));
        setSelectedRefDirAliases(selectedRefDirAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'ref_masked_fasta') {
        // If successful, update the state to remove the file from the list
        setSelectedRefMaskedFasta(selectedRefMaskedFasta.filter((_, index) => index !== indexToRemove));
        setSelectedRefMaskedFastaAliases(selectedRefMaskedFastaAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'ref_fai') {
        // If successful, update the state to remove the file from the list
        setSelectedRefFai(selectedRefFai.filter((_, index) => index !== indexToRemove));
        setSelectedRefFaiAliases(selectedRefFaiAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'ref_dict') {
        // If successful, update the state to remove the file from the list
        setSelectedRefDict(selectedRefDict.filter((_, index) => index !== indexToRemove));
        setSelectedRefDictAliases(selectedRefDictAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'ref_bwa') {
        // If successful, update the state to remove the file from the list
        setSelectedRefBwa(selectedRefBwa.filter((_, index) => index !== indexToRemove));
        setSelectedRefBwaAliases(selectedRefBwaAliases.filter((_, index) => index !== indexToRemove));
      } else if (target === 'snpeffconfig') {
        // If successful, update the state to remove the file from the list
        setSelectedSnpeffconfig(selectedSnpeffconfig.filter((_, Index) => Index !== indexToRemove));
        setSelectedSnpeffconfigAliases(selectedSnpeffconfigAliases.filter((_, Index) => Index !== indexToRemove));
      } else if (target === 'skip_samples_file') {
        // If successful, update the state to remove the file from the list
        setSelectedSkipSamplesFile(selectedSkipSamplesFile.filter((_, index) => index !== indexToRemove));
        setSelectedSkipSamplesFileAliases(selectedSkipSamplesFileAliases.filter((_, index) => index !== indexToRemove));
      } else {
        // If successful, update the state to remove the file from the list
        setSelectedFiles(selectedFiles.filter((_, index) => index !== indexToRemove));
        setSelectedAliases(selectedAliases.filter((_, index) => index !== indexToRemove));
      }
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
    if (fileManagerTarget === 'add_sra_file') {
      setActiveOpenFiles(selectedSraFiles);
    } else if (fileManagerTarget === 'add_vcf_file') {
      setActiveOpenFiles(selectedVcfFiles);
    } else if (fileManagerTarget === 'gambit_db') {
      setActiveOpenFiles(selectedGambitDb);
    } else if (fileManagerTarget === 'gambit_h5_dir') {
      setActiveOpenFiles(selectedGambitH5Dir);
    } else if (fileManagerTarget === 'subtype_db') {
      setActiveOpenFiles(selectedSubtypeDb);
    } else if (fileManagerTarget === 'fasta') {
      setActiveOpenFiles(selectedFastaFiles);
    } else if (fileManagerTarget === 'ref_dir') {
      setActiveOpenFiles(selectedRefDir);
    } else if (fileManagerTarget === 'ref_masked_fasta') {
      setActiveOpenFiles(selectedRefMaskedFasta);
    } else if (fileManagerTarget === 'ref_fai') {
      setActiveOpenFiles(selectedRefFai);
    } else if (fileManagerTarget === 'ref_dict') {
      setActiveOpenFiles(selectedRefDict);
    } else if (fileManagerTarget === 'ref_bwa') {
      setActiveOpenFiles(selectedRefBwa);
    } else if (fileManagerTarget === 'snpeffconfig') {
      setActiveOpenFiles(selectedSnpeffconfig);
    } else if (fileManagerTarget === 'skip_samples_file') {
      setActiveOpenFiles(selectedSkipSamplesFile);
    } else {
      setActiveOpenFiles(selectedFiles);
    }
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

    if (fileManagerTarget === 'add_sra_file') {
      setSelectedSraFiles(nextFiles);
    } else if (fileManagerTarget === 'add_vcf_file') {
      setSelectedVcfFiles(nextFiles);
    } else if (fileManagerTarget === 'gambit_db') {
      setSelectedGambitDb(nextFiles);
    } else if (fileManagerTarget === 'gambit_h5_dir') {
      setSelectedGambitH5Dir(nextFiles);
    } else if (fileManagerTarget === 'subtype_db') {
      setSelectedSubtypeDb(nextFiles);
    } else if (fileManagerTarget === 'fasta') {
      setSelectedFastaFiles(nextFiles);
    } else if (fileManagerTarget === 'ref_dir') {
      setSelectedRefDir(nextFiles);
    } else if (fileManagerTarget === 'ref_masked_fasta') {
      setSelectedRefMaskedFasta(nextFiles);
    } else if (fileManagerTarget === 'ref_fai') {
      setSelectedRefFai(nextFiles);
    } else if (fileManagerTarget === 'ref_dict') {
      setSelectedRefDict(nextFiles);
    } else if (fileManagerTarget === 'ref_bwa') {
      setSelectedRefBwa(nextFiles);
    } else if (fileManagerTarget === 'snpeffconfig') {
      setSelectedSnpeffconfig(nextFiles);
    } else if (fileManagerTarget === 'skip_samples_file') {
      setSelectedSkipSamplesFile(nextFiles);
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
      title: "What is **MycoSNP** (Whole Genome Sequencing Analysis of Fungal Isolates)?", "prompt": "What is **MycoSNP** (Whole Genome Sequencing Analysis of Fungal Isolates)?" });
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
    formData.method = props.selectedWorkflow?.toUpperCase();
    if (formData['description'] == null || formData['description'].trim() === '') {
      setHasMessage(true);
      setMessage("Please provide a job description before submitting the form.");
      setIsError(true);
      return;
    }
    if (!formData['workflow']) {
      formData['workflow'] = 'NFCORE_MYCOSNP'; // Set a default workflow if not provided
    }
    if (!formData['assembler']) {
      formData['pre_mycosnp_params']['assembler'] = 'skesa'; // Set a default assembler if not provided
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
              formData.outdir = directory + "/mycosnp";
            } else {
              formData['outdir'] = updatedFiles[0] + "/mycosnp";
            }
          }

          if (selectedSraFiles.length === 1 && isCSVFile(selectedSraFiles[0]) ) {
            let updatedSraFiles = selectedSraFiles.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['add_sra_file'] = updatedSraFiles[0];
          } else if (selectedSraFiles.length > 1){
            setHasMessage(true);
            setMessage("Please select ONLY ONE CSV file containing SRA ids to download from NCBI.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE CSV file containing SRA ids to download from NCBI.");
            return;
          }

          if (selectedVcfFiles.length === 1 && isCSVFile(selectedVcfFiles[0])) {
            let updatedVcfFiles = selectedVcfFiles.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['add_vcf_file'] = updatedVcfFiles[0];
          } else if (selectedVcfFiles.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE CSV file containing a list of file paths to vcf files generated from previous runs of this workflow to include in this analysis.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE CSV file containing a list of file paths to vcf files generated from previous runs of this workflow to include in this analysis.");
            return;
          }

          if (selectedGambitDb.length === 1 && isCSVFile(selectedGambitDb[0])) {
            let updatedGambitDb = selectedGambitDb.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['gambit_db'] = updatedGambitDb[0];
          } else if (selectedGambitDb.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE GAMBIT '.gdb' or '.db' file containing metadata.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE GAMBIT '.gdb' or '.db' file containing metadata.");
            return;
          }

          if (selectedGambitH5Dir.length === 1 && isFolderPath(selectedGambitH5Dir[0])) {
            let updatedGambitH5Dir = selectedGambitH5Dir.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['gambit_h5_dir'] = updatedGambitH5Dir[0];
          } else if (selectedGambitH5Dir.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE directory containing a GAMBIT signature file, which can be split into one or more chunks.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE directory containing a GAMBIT signature file, which can be split into one or more chunks.");
            return;
          }

          if (selectedSubtypeDb.length === 1 && isFolderPath(selectedSubtypeDb[0])) {
            let updatedSubtypeDb = selectedSubtypeDb.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['subtype_db'] = updatedSubtypeDb[0];
          } else if (selectedSubtypeDb.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE directory containing the subtyper files.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE directory containing the subtyper files.");
            return;
          }

          if (selectedFastaFiles.length === 1) {
            let updatedFastaFiles = selectedFastaFiles.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['fasta'] = updatedFastaFiles[0];
          } else if (selectedFastaFiles.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE fasta file.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE fasta file.");
            return;
          }

          if (selectedRefDir.length === 1 && isFolderPath(selectedRefDir[0])) {
            let updatedRefDir = selectedRefDir.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['ref_dir'] = updatedRefDir[0];
          } else if (selectedRefDir.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE directory containing reference genome masked files/picard dict/samtools fai/bwa index from previous run.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE directory containing reference genome masked files/picard dict/samtools fai/bwa index from previous run.");
            return;
          }

          if (selectedRefMaskedFasta.length === 1) {
            let updatedRefMaskedFasta = selectedRefMaskedFasta.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['ref_masked_fasta'] = updatedRefMaskedFasta[0];
          } else if (selectedRefMaskedFasta.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE reference genome masked fasta file.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE reference genome masked fasta file.");
            return;
          }

          if (selectedRefFai.length === 1) {
            let updatedRefFai = selectedRefFai.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['ref_fai'] = updatedRefFai[0];
          } else if (selectedRefFai.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE reference genome samtools fai file.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE reference genome samtools fai file.");
            return;
          }

          if (selectedRefDict.length === 1) {
            let updatedRefDict = selectedRefDict.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['ref_dict'] = updatedRefDict[0];
          } else if (selectedRefDict.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE reference genome picard tools dict file.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE reference genome picard tools dict file.");
            return;
          }

          if (selectedRefBwa.length === 1 && isFolderPath(selectedRefBwa[0])) {
            let updatedRefBwa = selectedRefBwa.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['ref_bwa'] = updatedRefBwa[0];
          } else if (selectedRefBwa.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE reference genome bwatools bwa directory.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE reference genome bwatools bwa directory.");
            return;
          }

          if (selectedSnpeffconfig.length === 1) {
            let updatedSnpeffconfig = selectedSnpeffconfig.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['snpeffconfig'] = updatedSnpeffconfig[0];
          } else if (selectedSnpeffconfig.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE directory containing the snpeff config.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE directory containing the snpeff config file.");
            return;
          }

          if (selectedSkipSamplesFile.length === 1 && isFolderPath(selectedSkipSamplesFile[0])) {
            let updatedSkipSamplesFile = selectedSkipSamplesFile.map(file => "/usr/src/app/storage/" + authData.username + file);
            formData['pre_mycosnp_params']['skip_samples_file'] = updatedSkipSamplesFile[0];
          } else if (selectedSkipSamplesFile.length > 1) {
            setHasMessage(true);
            setMessage("Please select ONLY ONE file with new-line separated list of IDs to skip for variant calling and analysis.");
            setIsError(true);
            setFormErrors("Please select ONLY ONE file with new-line separated list of IDs to skip for variant calling and analysis.");
            return;
          }

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
                You will need to create a samplesheet with information about the samples you would like to analyse before running either workflow. Use the --input parameter to specify its location. It has to be a comma-separated file (csv) with at least 3 columns, and a header row as shown in the examples below. An example samplesheet could be downloaded from the following link: <a href="https://github.com/CDCgov/mycosnp-nf/blob/b4216fec5d1ec2e7d8b136539748c3623a02fb45/assets/samplesheet.csv" target="_blank" rel="noopener noreferrer">Download Example Samplesheet</a>
                  </p>
                    <pre><code>
                      sample,fastq_1,fastq_2 <br />
                      CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz <br />
                      CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz <br />
                      CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz <br />
                      TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz <br />
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
            selectedFiles={activeOpenFiles}
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
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile("input", item, index)}>
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
            You will need to create a samplesheet with information about the samples you would like to analyse before running either workflow. Use the --input parameter to specify its location. It has to be a comma-separated file (csv) with at least 3 columns, and a header row as shown in the examples below. An example samplesheet could be downloaded from the following link: <a href="https://github.com/CDCgov/mycosnp-nf/blob/b4216fec5d1ec2e7d8b136539748c3623a02fb45/assets/samplesheet.csv" target="_blank" rel="noopener noreferrer">Download Example Samplesheet</a>
          </p>
            <pre><code>
              sample,fastq_1,fastq_2 <br />
              CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz <br />
              CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz <br />
              CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz <br />
              TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz <br />
            </code></pre>
          <p>If you choose a Directory, the system will create a samplesheet from the directory of fastq files. The script will search 1 directory deep and attempt to determine sample id names and pairing/multilane information and will automatically create a samplesheet.</p>
        </div>

        <div>
          <div id="upload-sra-data-div">
            <b>add_sra_file</b> <br />
            {selectedSraFiles.length > 0 && (
              <div id="sra-files-selected">
                <ScrollableListContainer>
                  <List dense>
                    <>
                      {selectedSraFiles.map((item, index) => {
                        return (
                          <div key={index} className="file-selections">
                            <CustomListItem key={index}>
                              <IconButton edge="start" aria-label="delete" onClick={() => removeFile("add_sra_file", item, index)}>
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
                    Notice: Please select ONLY ONE CSV file containing SRA ids to download from NCBI.
                  </div>
                }
              </div>
            )}
            <div className="folder-button">
              <button type="button" onClick={() => openFileManager('add_sra_file')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
                <FontAwesomeIcon icon={faFolderOpen} />
              </button>
            </div>
          </div>
          <div>
            <p>Choose a comma-separated file containing SRA ids to download from NCBI. Format: Name, SRAID</p>
            <p>You may provide a list of SRA ids as additional input sequences into either workflow. Use the --add_sra_file parameter to specify its location. It has to be a comma-separated file (csv) with one or 2 columns, and NO header row as shown in the examples below. </p>
            <p>If two fields are provided, the first field will be used as the sequence name, and the second field will be used to specify the SRA id to download. If one field is provided, it must be the SRA id, and this SRA id will be used as the sequence name.</p>
            <pre><code>
              B12352,SRR7909282<br />
              SRR7909249<br />
              B13520,SRR7909394<br />
            </code></pre>
          </div>
        </div>

        <div>
          <div id="upload-vcf-data-div">
            <b>add_vcf_file</b> <br />
            {selectedVcfFiles.length > 0 && (
              <div id="vcf-files-selected">
                <ScrollableListContainer>
                  <List dense>
                    <>
                      {selectedVcfFiles.map((item, index) => {
                        return (
                          <div key={index} className="file-selections">
                            <CustomListItem key={index}>
                              <IconButton edge="start" aria-label="delete" onClick={() => removeFile("add_vcf_file", item, index)}>
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
                {selectedVcfFiles && selectedVcfFiles.length > 1 &&
                  <div style={{ color: 'red' }}>
                    Notice: Please select ONLY ONE CSV file.
                  </div>
                }
              </div>
            )}
            <div className="folder-button">
              <button type="button" onClick={() => openFileManager('add_vcf_file')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
                <FontAwesomeIcon icon={faFolderOpen} />
              </button>
            </div>
          </div>
          <div>
          <p>(<strong>Main MycoSNP workflow only</strong>) Path to .csv (only one column) containing a list of file paths to vcf files generated from previous runs of this workflow to include in this analysis. They must use the same exact reference. *.tbi file must be in same location. Each line of .csv file has format: /path/to/vcf/file.gz</p>
          </div>
        </div>

      </div>

      <div className="separator heading">
        <div className="stripe"></div>
        <h2 className="h-sm font-weight-bold">
          Pre-MycoSNP workflow options (Optional)
        </h2>
        <div className="stripe"></div>
      </div>

      <div>
        <div id="upload-gambit_db-data-div">
          <b>gambit_db</b> <br />
          {selectedGambitDb.length > 0 && (
            <div id="gambit_db-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedGambitDb.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('gambit_db', item, index)}>
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
              {selectedGambitDb && selectedGambitDb.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE GAMBIT '.gdb' or '.db' file containing metadata.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('gambit_db')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
        <p>Path to the GAMBIT '.gdb' or '.db' file containing metadata.</p>
        </div>
      </div>

      <div>
        <div id="upload-gambit_h5_dir-data-div">
          <b>gambit_h5_dir</b> <br />
          {selectedGambitH5Dir.length > 0 && (
            <div id="gambit_h5_dir-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedGambitH5Dir.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('gambit_h5_dir', item, index)}>
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
              {selectedGambitH5Dir && selectedGambitH5Dir.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE directory containing a GAMBIT signature file, which can be split into one or more chunks.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('gambit_h5_dir')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Full (absolute) path to the directory containing a GAMBIT signature file, which can be split into one or more chunks. Each chunk may be gzipped if preferred. Each chunk must have one of the following extensions: '.gs', '.h5', '.gs.gz' or '.h5.gz'. Refer to Pre-MycoSNP's default database in assets/gambit_db/signatures/ for an example (file is split so each chunk is below GitHub's individual file size limit of 100 MB and warning limit of 50 MB).</p>
        </div>
      </div>

      <div>
        <div id="upload-subtype_db-data-div">
          <b>subtype_db</b> <br />
          {selectedSubtypeDb.length > 0 && (
            <div id="subtype_db-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedSubtypeDb.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('subtype_db', item, index)}>
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
              {selectedSubtypeDb && selectedSubtypeDb.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE directory containing the subtyper files.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('subtype_db')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to a directory containing the subtyper files. This directory should contain sourmash signature files, each containing multiple sketches for the representative subtypes. This directory should also contain a csv file called sourmash_taxa.csv mapping each taxon name to a sourmash signature file. Taxon names must be the same as what is reported by GAMBIT.</p>
        </div>
      </div>
      

      <div className="separator heading">
        <div className="stripe"></div>
        <h2 className="h-sm font-weight-bold">
          Main MycoSNP workflow: reference genome options (Optional)
        </h2>
        <div className="stripe"></div>
      </div>

      <div>
        <div id="upload-fasta-data-div">
          <b>fasta</b> <br />
          {selectedFastaFiles.length > 0 && (
            <div id="fasta-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedFastaFiles.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('fasta', item, index)}>
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
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('fasta')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to FASTA formatted reference genome file.</p>
        </div>
      </div>

      <div>
        <div id="upload-ref_dir-data-div">
          <b>ref_dir</b> <br />
          {selectedRefDir.length > 0 && (
            <div id="ref_dir-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedRefDir.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('ref_dir', item, index)}>
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
              {selectedRefDir && selectedRefDir.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE directory containing reference genome masked files/picard dict/samtools fai/bwa index from previous run.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('ref_dir')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to reference genome masked files/picard dict/samtools fai/bwa index from previous run.</p>
        </div>
      </div>

      <div>
        <div id="upload-ref_masked_fasta-data-div">
          <b>ref_masked_fasta</b> <br />
          {selectedRefMaskedFasta.length > 0 && (
            <div id="ref_masked_fasta-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedRefMaskedFasta.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('ref_masked_fasta', item, index)}>
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
              {selectedRefMaskedFasta && selectedRefMaskedFasta.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE reference genome masked fasta file.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('ref_masked_fasta')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to reference genome masked fasta file.</p>
        </div>
      </div>

      <div>
        <div id="upload-ref_fai-data-div">
          <b>ref_fai</b> <br />
          {selectedRefFai.length > 0 && (
            <div id="ref_fai-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedRefFai.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('ref_fai', item, index)}>
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
              {selectedRefFai && selectedRefFai.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE reference genome samtools fai file.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('ref_fai')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to reference genome samtools fai file.</p>
        </div>
      </div>

      <div>
        <div id="upload-ref_dict-data-div">
          <b>ref_dict</b> <br />
          {selectedRefDict.length > 0 && (
            <div id="ref_dict-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedRefDict.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('ref_dict', item, index)}>
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
              {selectedRefDict && selectedRefDict.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE reference genome picard tools dict file.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('ref_dict')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to reference genome picard tools dict file.</p>
        </div>
      </div>

      <div>
        <div id="upload-ref_bwa-data-div">
          <b>ref_bwa</b> <br />
          {selectedRefBwa.length > 0 && (
            <div id="ref_bwa-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedRefBwa.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('ref_bwa', item, index)}>
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
              {selectedRefBwa && selectedRefBwa.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE reference genome bwatools bwa directory.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('ref_bwa')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to reference genome bwatools bwa directory.</p>
        </div>
      </div>

      <div>
        <div id="upload-snpeffconfig-data-div">
          <b>snpeffconfig</b> <br />
          {selectedSnpeffconfig.length > 0 && (
            <div id="snpeffconfig-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedSnpeffconfig.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('snpeffconfig', item, index)}>
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
              {selectedSnpeffconfig && selectedSnpeffconfig.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE snpeff config file.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('snpeffconfig')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>Path to snpeff config.</p>
        </div>
      </div>


      <div className="separator heading">
        <div className="stripe"></div>
        <h2 className="h-sm font-weight-bold">
          Main MycoSNP workflow: Run/Save/Skip options (Optional)
        </h2>
        <div className="stripe"></div>
      </div>

      <div>
        <div id="upload-skip_samples_file-data-div">
          <b>skip_samples_file</b> <br />
          {selectedSkipSamplesFile.length > 0 && (
            <div id="skip_samples_file-files-selected">
              <ScrollableListContainer>
                <List dense>
                  <>
                    {selectedSkipSamplesFile.map((item, index) => {
                      return (
                        <div key={index} className="file-selections">
                          <CustomListItem key={index}>
                            <IconButton edge="start" aria-label="delete" onClick={() => removeFile('skip_samples_file', item, index)}>
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
              {selectedSkipSamplesFile && selectedSkipSamplesFile.length > 1 &&
                <div style={{ color: 'red' }}>
                  Notice: Please select ONLY ONE file with new-line separated list of IDs to skip for variant calling and analysis.
                </div>
              }
            </div>
          )}
          <div className="folder-button">
            <button type="button" onClick={() => openFileManager('skip_samples_file')} style={{ fontSize: "1em", padding: "10px", borderRadius: "5px" }}>
              <FontAwesomeIcon icon={faFolderOpen} />
            </button>
          </div>
        </div>
        <div>
          <p>A file with new-line separated list of IDs to skip for variant calling and analysis.</p>
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