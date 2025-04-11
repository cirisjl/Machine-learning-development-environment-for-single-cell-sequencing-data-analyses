import React, { useState, useEffect } from 'react';
import InputDataComponent from '../../Tools/inputDataCollection';
import Schema from '../../../../schema/react-json-schema/Workflows/clusteringWorkflow.json'
import { uiSchema } from '../../../../schema/UI-schema/Workflows/clusteringWorkflow';
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
import {CELERY_BACKEND_API} from '../../../../constants/declarations';
import { useNavigate } from 'react-router-dom';
import { isUserAuth, getCookie } from '../../../../utils/utilFunctions';


export function ClusteringWorkFlowComponent(props) {
    const [selectedDatasets, setSelectedDatasets] = useState({});
    const [formErrors, setFormErrors] = useState("");
    const [value, setValue] = useState("");
    const [formData, setFormData] = useState({});

    const [ message, setMessage ] = useState('');
    const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);
    const [ isError, setIsError ] = useState(false);
    const [loading, setLoading] = useState(false);

    const navigate = useNavigate();

    const extractDir =  (inputFile) => {
      const fileLocParts = inputFile.split('/');
      fileLocParts.pop(); // Remove the file name from the array
      const output = fileLocParts.join('/'); // Join the remaining parts with '/'
      return output;
  };

    const filterCategoryMap = {
      clustering: '/workflows/clustering',
      // Add more filter categories and their corresponding URL paths as needed
    };

    const handleChange = (event) => {
        setValue(event.target.value);
      };

    const onDeleteDataset = (id) => {
    const currentSelectedDatasets = { ...selectedDatasets};
    
    if (currentSelectedDatasets[id]) {
        delete currentSelectedDatasets[id];
    }
    setSelectedDatasets(currentSelectedDatasets);
    };

    const onSelectDataset = (dataset) => {
        let datasetId = dataset.Id; 
        let currentSelectedDatasets = { ...selectedDatasets};
    
        if(currentSelectedDatasets[datasetId]) {
          delete currentSelectedDatasets[datasetId];
        } else {
          if(props.selectedWorkflow !== "integration") {
            currentSelectedDatasets = {};
          }
          currentSelectedDatasets[datasetId] = dataset;
        }
      setSelectedDatasets(currentSelectedDatasets)
    };

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

    // Function to handle selection of sub-items
  const onSelectSubItem = (mainItem, subItem) => {
    const mainItemId = mainItem.Id;
    let currentSelectedDatasets = { ...selectedDatasets };
  
    // Check if the main item is already selected
    if (currentSelectedDatasets[mainItemId]) {
        // If sub-item is already selected, deselect it
        if (currentSelectedDatasets[mainItemId].selectedSubItem?.process_id  === subItem.process_id ) {
            delete currentSelectedDatasets[mainItemId];
        } else {
            // Update the selected main item with the selected sub-item
            currentSelectedDatasets[mainItemId] = {
                ...mainItem,
                selectedSubItem: subItem
            };
        }
    } else {
        // Select the main item and the sub-item
        currentSelectedDatasets = {
            [mainItemId]: {
                ...mainItem,
                selectedSubItem: subItem
            }
        };
    }
  
    setSelectedDatasets(currentSelectedDatasets);
    };
    const handleSubmit = ({ formData }) => {
        formData = formData.parameters;
        console.log("In submit function");
        console.log("Submitted data:", formData);
        console.log(selectedDatasets);
        // Verify the authenticity of the user
        isUserAuth(getCookie('jwtToken'))
        .then((authData) => {
          if (authData.isAuth) {
            // Prepare the data to send. For workflow we only use a single dataset for clustering.
            //TODO - Change this logic accordingly based on other workflows in the future
            formData.userID = authData.username;
            console.log("Adding user data");
            if(Object.keys(selectedDatasets).length > 0) {
              const dataset = Object.values(selectedDatasets)[0]; // Assuming single dataset
              formData.dataset = dataset.Title;
              if (dataset.inputFiles.length > 1) {
                formData.input = extractDir(dataset.inputFiles[0]);
                formData.output = formData.input + "/Results";
              } else if (dataset.inputFiles.length === 1) {
                formData.input = dataset.inputFiles[0];
                const directory = extractDir(formData.input);
                formData.output = directory + "/Results";
              }
            }
            console.log(props.selectedWorkflow);
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
              navigate("/mydata/taskDetails", { state: { job_id: jobId, method: "", datasetURL: formData.input, description: "Clustering workflow", process: "Clustering" } });
              // let datasetName = "";
              // if (typeof formData.dataset === 'string') {
              //   datasetName = formData.dataset
              // } else if (Array.isArray(formData.dataset)) {
              //   if (formData.dataset.length > 1) {
              //     datasetName = formData.dataset.join('_');
              //   } else if (formData.dataset.length === 1) {
              //     datasetName = formData.dataset[0];
              //   }
              // }

              //   // After a successfull task creation, store the intermediate task information in the mongoDB task_results collection
              //   const taskId = response.task_id;
              //   let method = "";

              //   if(filterCategory === "reduction") {
              //     method = "Reduction";
              //   } else if(filterCategory === "formatting") {
              //     method = "Formatting";
              //   } else if(parametersKey[filterCategory]) {
              //     method = formData[parametersKey[filterCategory]].methods[0];
              //   } else {
              //     method = formData.methods[0];
              //   }

              //   const output = formData.output;

              //   // Make API call to store the task information
              //   const requestBody = {
              //     datasetTitle: formData.dataset,
              //     taskId: taskId,
              //     method: method,
              //     datasetURL: formData.input,
              //     tool: filterCategory,
              //     outputPath: output,
              //     Owner: authData.username,
              //     status: 'Processing'
              //   };

                // fetch(`${SERVER_URL}/createTask`, {
                //   method: 'POST',
                //   headers: {
                //     'Content-Type': 'application/json'
                //   },
                //   body: JSON.stringify(requestBody)
                // })
                // .then(response => {
                //   if (response.ok) {
                //     if (response.status === 200) {
                //       console.log('Task created successfully!');
                //       setLoading(false);
                //       setSuccessMessage('Form submitted successfully!');
                //       setErrorMessage('');
                //       navigate("/mydata/taskDetails", { state: { taskId: taskId, method: method, datasetURL: formData.input, datasetTitle: formData.dataset, tool: filterStaticCategoryMap[filterCategory] } });
                //     } else if (response.status === 400) {
                //       response.json().then(data => {
                //         console.error('Validation error:', data.error);
                //         setErrorMessage(data.error); // Set specific error message based on response
                //       });
                //     } else {
                //       throw new Error('Unexpected response from server');
                //     }
                //   } else {
                //     throw new Error('Failed to submit form: ' + response.statusText);
                //   }
                // })
                // .catch(error => {
                //   console.error('API request failed:', error);
                //   setLoading(false);
                //   setHasMessage(true);
                //   setMessage('An error occurred while submitting the form: ' + error.message);
                //   setIsError(true);
                // });
            })
            .catch(error => {
              // Handle any errors that occur during the API call
              console.error("Form submission error:", error);
              setLoading(false);
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

      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} isError={isError}/>}

      <div className="separator heading">
        <div className="stripe"></div> 
          <h2 className="h-sm font-weight-bold">
            Tool Parameters 
          </h2> 
        <div className="stripe"></div>
      </div>

      <div>
        <InputDataComponent formErrors={formErrors} filterCategory={props.selectedWorkflow} selectedDatasets={selectedDatasets}
            onSelectDataset={onSelectDataset} onDeleteDataset={onDeleteDataset} onSelectSubItem={onSelectSubItem}/>
      </div>

      <div className="form-component">
        <Form
            schema = {Schema}
            uiSchema={uiSchema}
            widgets={widgets}
            formData={formData}
            onChange={({ formData }) => setFormData(formData)}
            onSubmit={handleSubmit}
            onError={(errors) => console.log("Form errors:", errors)}
            validator={validator}
        />
      </div>

    </div>
  )
};