import React, { useState, useEffect, useRef } from 'react';
import { getCookie, isUserAuth} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Form from 'react-jsonschema-form';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import InputDataComponent from './inputDataCollection';
import { CELERY_BACKEND_API, SERVER_URL, defaultValues, WEB_SOCKET_URL, defaultQcParams,defaultNormalizationParams ,defaultReductionParams} from '../../../constants/declarations';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { faSpinner } from '@fortawesome/free-solid-svg-icons';
import GeneRangeSlider from './components/geneRangeSlider';
import RangeSlider from './components/sliderComponent';
import SwitchComponent from './components/switchComponent';
import UseDefaultSwitch from './components/useDefaultSwitch';
import MultiSelectComponent from './components/multiselectComponent';
import SelectComponent from './components/selectComponent';
import ClusterLabelInput from './components/customInputComponent';

export default function ToolsDetailsComponent(props) {
    const filterName = props.filter;
    const filterCategory = props.category;
    const [selectedDatasets, setSelectedDatasets] = useState({});
    const [shouldHideForSeurat, setShouldHideForSeurat] = useState(false);
    const [useDefaultValue, setUseDefaultValue] = useState(false);

    const filterCategoryMap = {
      quality_control: '/api/tools/qc',
      normalization: '/api/tools/normalize',
      imputation: '/api/tools/impute',
      integration: '/api/tools/integrate',
      evaluation: '/api/tools/evaluate',
      formatting: '/api/tools/convert',
      reduction: '/api/tools/reduce'
      // Add more filter categories and their corresponding URL paths as needed
    };

    const parametersKey = {
      quality_control: 'qc_params',
      normalization: 'normalization_params',
      imputation: 'imputation_params',
      reduction: 'reduction_params'
    };

    const filterStaticCategoryMap = {
      quality_control: 'Quality Control',
      normalization: 'Normalization',
      imputation: 'Imputation',
      integration: 'Integration',
      evaluation: 'Evaluation',
      formatting: 'Formatting',
      reduction: 'Reduction'
      // Add more filter categories and their corresponding Names as needed
    };

    console.log(filterName);

    let jwtToken = getCookie('jwtToken');
    const [formData, setFormData] = useState({});
    const [useDefault, setUseDefault] = useState(true);
    // let useDefault = true;
    const [filterSchema, setFilterSchema] = useState(null);
    const [UIfilterSchema, setUIFilterSchema] = useState(null);
    const [selectedDataset, setSelectedDataset] = useState([]);
    const [selectedOptions, setSelectedOptions] = useState([]);
    const [formErrors, setFormErrors] = useState("");

    const [loading, setLoading] = useState(false);
    const [successMessage, setSuccessMessage] = useState('');
    const [errorMessage, setErrorMessage] = useState('');

    const navigate = useNavigate();
    
    const onSelectDataset = (dataset) => {
      let datasetId = dataset.Id; 
      let currentSelectedDatasets = { ...selectedDatasets};
  
      if(currentSelectedDatasets[datasetId]) {
        delete currentSelectedDatasets[datasetId];
      } else {
        if(filterCategory !== "integration") {
          currentSelectedDatasets = {};
        }
        currentSelectedDatasets[datasetId] = dataset;
      }
      if(filterCategory === "quality_control") {
        // Check if any of the selected datasets should trigger hiding for Seurat
        const shouldHideForSeurat = Object.values(currentSelectedDatasets).some(dataset =>
          dataset.inputFiles.length === 1 &&
          (dataset.inputFiles[0].toLowerCase().endsWith('h5seurat') ||
          dataset.inputFiles[0].toLowerCase().endsWith('rds') ||
          dataset.inputFiles[0].toLowerCase().endsWith('robj'))
        );
        setShouldHideForSeurat(shouldHideForSeurat);
      }
    setSelectedDatasets(currentSelectedDatasets)
  };

  const onDeleteDataset = (id) => {
    const currentSelectedDatasets = { ...selectedDatasets};
  
    if (currentSelectedDatasets[id]) {
        delete currentSelectedDatasets[id];
    }
    setSelectedDatasets(currentSelectedDatasets);
  };

    const extractDir = (inputFile) => {
        const fileLocParts = inputFile.split('/');
        fileLocParts.pop(); // Remove the file name from the array
        const output = fileLocParts.join('/'); // Join the remaining parts with '/'
        return output;
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
        UseDefaultSwitch: UseDefaultSwitch,
        ClusterLabelInput: ClusterLabelInput
      };

      const onSubmit = ({ formData }) => {

        // Handle form submission here
        formData = formData.parameters;

        // Perform form validation and set formErrors accordingly
        if(filterCategory === "integration" && Object.keys(selectedDatasets).length < 2) {
          setFormErrors("Please select atleast two datasets for integration before submitting the form");
          console.log("Failed to submit the form");
        }
        else if(filterCategory !== "integration" && Object.keys(selectedDatasets).length === 0) {
          setFormErrors("Please select a dataset before submitting the form");
          console.log("Failed to submit the form");
        } else {
          setLoading(true);
          if(filterCategory === "integration") {
            const datasetsArray = Object.values(selectedDatasets);
            const titlesArray = datasetsArray.map(dataset => dataset.Title);
            formData.dataset = titlesArray;

             let inputArray = datasetsArray.map(dataset => {
              if (dataset.inputFiles.length > 1) {
                return extractDir(dataset.inputFiles[0]); // Assuming extractDir is a function defined elsewhere
              } else {
                return dataset.inputFiles[0];
              }
            });

            formData.input = inputArray;
            formData.output = "/IntegrationResults";

          } else {
              const dataset = Object.values(selectedDatasets)[0]; // Assuming single dataset for non-integration category
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
      
          // Verify the authenticity of the user
          isUserAuth(jwtToken)
          .then((authData) => {
            if (authData.isAuth) {
              formData.userID = authData.username;

              const RELATIVE_PATH = filterCategoryMap[filterCategory];
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
                
                let datasetName = "";
                if (typeof formData.dataset === 'string') {
                  datasetName = formData.dataset
                } else if (Array.isArray(formData.dataset)) {
                  if (formData.dataset.length > 1) {
                    datasetName = formData.dataset.join('_');
                  } else if (formData.dataset.length === 1) {
                    datasetName = formData.dataset[0];
                  }
                }

                  // After a successfull task creation, store the intermediate task information in the mongoDB task_results collection
                  const taskId = response.task_id;
                  let method = "";

                  if(filterCategory === "reduction") {
                    method = "Reduction";
                  } else if(filterCategory === "formatting") {
                    method = "Formatting";
                  } else if(parametersKey[filterCategory]) {
                    method = formData[parametersKey[filterCategory]].methods[0];
                  } else {
                    method = formData.methods[0];
                  }

                  const output = formData.output;

                  // Make API call to store the task information
                  const requestBody = {
                    datasetTitle: formData.dataset,
                    taskId: taskId,
                    method: method,
                    datasetURL: formData.input,
                    tool: filterCategory,
                    outputPath: output,
                    Owner: authData.username,
                    status: 'Processing'
                  };

                  fetch(`${SERVER_URL}/createTask`, {
                    method: 'POST',
                    headers: {
                      'Content-Type': 'application/json'
                    },
                    body: JSON.stringify(requestBody)
                  })
                  .then(response => {
                    if (response.ok) {
                      if (response.status === 200) {
                        console.log('Task created successfully!');
                        setLoading(false);
                        setSuccessMessage('Form submitted successfully!');
                        setErrorMessage('');
                        navigate("/mydata/taskDetails", { state: { taskId: taskId, method: method, datasetURL: formData.input, datasetTitle: formData.dataset, tool: filterStaticCategoryMap[filterCategory] } });
                      } else if (response.status === 400) {
                        response.json().then(data => {
                          console.error('Validation error:', data.error);
                          setErrorMessage(data.error); // Set specific error message based on response
                        });
                      } else {
                        throw new Error('Unexpected response from server');
                      }
                    } else {
                      throw new Error('Failed to submit form: ' + response.statusText);
                    }
                  })
                  .catch(error => {
                    console.error('API request failed:', error);
                    setLoading(false);
                    setSuccessMessage('');
                    setErrorMessage('An error occurred while submitting the form: ' + error.message);
                  });

                  setFormErrors("");
              })
              .catch(error => {
                // Handle any errors that occur during the API call
                console.error("Form submission error:", error);
                setLoading(false);
                setSuccessMessage('');
                setErrorMessage('An error occurred while submitting the form.');
              });
            } else {
              console.warn("Unauthorized - please login first to continue");
              navigate("/routing");
            }
          })
          .catch((error) => {
            console.error(error);
            setLoading(false);
            setSuccessMessage('');
            setErrorMessage('An error occurred while submitting the form.');
          } 
          );
          }
      };

      useEffect(() => {
        const timeoutId = setTimeout(() => {
          setSuccessMessage('');
        }, 3000);
        // Return a cleanup function to cancel the timeout when the component unmounts
        return () => clearTimeout(timeoutId);
      }, [successMessage]);

      useEffect(() => {
        const timeoutId = setTimeout(() => {
          setErrorMessage('');
        }, 3000);
        // Return a cleanup function to cancel the timeout when the component unmounts
        return () => clearTimeout(timeoutId);
      }, [errorMessage]);

      useEffect(() => {
        const timeoutId = setTimeout(() => {
          setFormErrors('');
        }, 3000);
        // Return a cleanup function to cancel the timeout when the component unmounts
        return () => clearTimeout(timeoutId);
      }, [formErrors]);


      const handleDatasetChange = event => {
        let value = event.target.value;
        if(value !== "") {
          setSelectedDataset(event.target.value);
        } else {
          setSelectedDataset([]);
        }
      };

      const handleMultipleDatasetChange = (event) => {
        const selectedValues = Array.from(event.target.selectedOptions)
            .filter(option => option.value !== "") // Filter out options with empty value
            .map(option => option.value);
    
        setSelectedOptions(selectedValues);
        console.log(selectedOptions);
    };


  useEffect(() => {
    import(`./../../../schema/react-json-schema/Tools/${filterCategory}/${filterName}.json`)
    .then((module) => {
      setFilterSchema(JSON.parse(JSON.stringify(module.default)));
      // setFilterSchema({...module.default});

      console.log("react json schema")
      console.log(filterSchema);
    })
    .catch((error) => {
      console.error('Error loading filter schema:', error);
      setFilterSchema(null);
    });

    import(`./../../../schema/UI-schema/Tools/${filterCategory}/${filterName}.js`)
    .then((module) => {
      setUIFilterSchema(JSON.parse(JSON.stringify(module.uiSchema)));
      // setUIFilterSchema({...module.uiSchema});
      console.log("react json UI schema")
      console.log(UIfilterSchema);
    })
    .catch((error) => {
      console.error('Error loading UI filter schema:', error);
      setUIFilterSchema(null);
    });
  },[filterName,filterCategory]);

  const handleChange = ({ formData }) => {

    // Choose the appropriate default parameters based on the category
    let defaultParams;
    let paramsKey;
    if (filterCategory === 'quality_control') {
        defaultParams = defaultQcParams;
        paramsKey = 'qc_params';

    } else if (filterCategory === 'normalization') {
        defaultParams = defaultNormalizationParams;
        paramsKey = 'normalization_params';

    } else if(filterCategory === 'reduction') {
        defaultParams = defaultReductionParams;
        paramsKey = 'reduction_params';
    } else {
        // For other categories, we do not reset parameters, just set formData
        setFormData(formData);
        return;
    }

    const currentToolParams = formData.parameters[paramsKey] || {};

    // Determine if there's a change in the use_default toggle
    const useDefaultChanged = useDefault !== currentToolParams.use_default;

    // Check for any changes in default parameters
    let defaultParamsChanged = Object.keys(defaultParams).some(key => {
        return JSON.stringify(currentToolParams[key]) !== JSON.stringify(defaultParams[key]);
    });

    if (useDefaultChanged) {
        if (currentToolParams.use_default) {
            // If use_default is toggled to true, reset only the default parameters
            const resetParams = {};
            Object.keys(defaultParams).forEach(key => {
                resetParams[key] = defaultParams[key];
            });
            formData.parameters[paramsKey] = {
              ...formData.parameters[paramsKey],
              ...resetParams
          };
        }
    } else if (defaultParamsChanged) {
        // If any default parameters have changed and use_default was previously true, set it to false
        formData.parameters[paramsKey].use_default = false;
    }

    setUseDefault(formData.parameters[paramsKey].use_default);
    setFormData(formData);
};

  return (
    <div className='tools-container common-class-tools-and-workflows'>
         {formErrors && (
            <div className='message-box' style={{ backgroundColor: '#f0c0c0' }}>
              <div style={{ textAlign: 'center' }}>
                <p>{formErrors}</p>    
              </div>
            </div>
          )}
         
          {loading && (
            <div className='message-box loadingIcon' style={{ backgroundColor: '#bdf0c0' }}>
              <div style={{ textAlign: 'center' }}>
                <FontAwesomeIcon icon={faSpinner} spin />
                <p>Loading...</p>       
              </div>
            </div>
          )}
          {successMessage && (  
            <div className='message-box success' id="tooltip" style={{ backgroundColor: '#bdf0c0' }}>
              <div style={{ textAlign: 'center' }}>
                <p>{successMessage}</p>       
              </div>
            </div>
          )}
          {errorMessage && (
            <div className='message-box error' id="tooltip" style={{ backgroundColor: '#f0c0c0' }}>
              <div style={{ textAlign: 'center' }}>
                <p>{errorMessage}</p>       
              </div>
            </div>
          )}
      <div className="separator heading">
        <div className="stripe"></div> 
          <h2 className="h-sm font-weight-bold">
            Tool Parameters 
          </h2> 
        <div className="stripe"></div>
      </div>
      <div>
        <InputDataComponent handleDatasetChange={handleDatasetChange} handleMultipleDatasetChange={handleMultipleDatasetChange} 
        formErrors={formErrors} filterCategory={filterCategory} filterName={filterName} selectedDatasets={selectedDatasets}
        onSelectDataset={onSelectDataset} onDeleteDataset={onDeleteDataset}/>
      </div>
            
        {filterSchema && UIfilterSchema ? (
          <div className="form-component">
            <Form
            schema={filterSchema}
            formData={formData}
            widgets={widgets}
            onChange={handleChange}
            uiSchema={UIfilterSchema}
            onSubmit={onSubmit}
            key={JSON.stringify(formData)} // Helps in re-rendering the form with updated data
        />
        </div>
          ) : (
            <div>No Schema for this tool.</div>
          )}
    </div>
  )
};
