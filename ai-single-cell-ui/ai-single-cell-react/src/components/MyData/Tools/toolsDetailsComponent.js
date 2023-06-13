import React, { useState, useEffect } from 'react';
import { getCookie, isUserAuth} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
// import schema from '../../../react-json-schema/Tools/normalizeUsingScanpySchema.json';
import Form from 'react-jsonschema-form';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import InputDataComponent from './inputDataCollection';
import { CELERY_BACKEND_API, SERVER_URL } from '../../../constants/declarations';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { faSpinner } from '@fortawesome/free-solid-svg-icons';


export default function ToolsDetailsComponent(props) {
    const filterName = props.filter;

    console.log(filterName);

    let jwtToken = getCookie('jwtToken');
    const [formData, setFormData] = useState({});
    const [filterSchema, setFilterSchema] = useState(null);
    const [selectedDataset, setSelectedDataset] = useState([]);
    const [formErrors, setFormErrors] = useState("");

    const [loading, setLoading] = useState(false);
    const [successMessage, setSuccessMessage] = useState('');
    const [errorMessage, setErrorMessage] = useState('');

    const navigate = useNavigate();

    const extractDir =  (inputFile) => {
        const fileLocParts = inputFile.split('/');
        fileLocParts.pop(); // Remove the file name from the array
        const output = fileLocParts.join('/'); // Join the remaining parts with '/'
        return output;
    };

    const uiSchema = {
      "parameters": {
        "classNames": "category",
          "output_format": {
            "classNames": "sub-category",
            "ui:widget": "select",
            "ui:placeholder": "Select file format"
          },
          "methods": {
            "classNames": "sub-category",
            // "ui:widget": "select",
            "ui:placeholder": "Select a method",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <select>
                  <option>Scanpy</option>
                </select>
          </div>
            ),
          },
          "default_assay": {
            "classNames": "sub-category",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <span data-v-22825496="" class="ui-form-title-message warning"> * Optional </span>
                <input type='text' />
          </div>
            ),
          },
          "layer": {
            "classNames": "sub-category"
          },
          "path_of_scrublet_calls": {
            "classNames": "sub-category"
          },
          "species": {
            "classNames": "sub-category",
            "ui:placeholder": "Select species type"
          },
          "idtype": {
            "classNames": "sub-category"
          },
          "genes": {
            "classNames": "sub-category",
          },
          "ncores": {
            "classNames": "sub-category",
            "ui:widget": "range",
          },
          "show_umap": {
            "classNames": "sub-category",
            "ui:widget": "toggle"
          },
          "show_error": {
            "classNames": "sub-category",
            "ui:widget": "toggle"
          }
      }
    };
    
      const widgets = {
        toggle: (props) => (
          <Toggle
            checked={props.value}
            onChange={(e) => props.onChange(e.target.checked)}
          />
        ),
      };

      const onSubmit = ({ formData }) => {

        // Handle form submission here
        formData = formData.parameters;

        // Perform form validation and set formErrors accordingly
        if(selectedDataset.length === 0) {
          setFormErrors("Please select a dataset before submitting the form");
          console.log("Failed to submit the form");
        } else {
          setLoading(true);

          console.log(selectedDataset);
            const parsedSelectedDataset = JSON.parse(selectedDataset);
            formData.dataset = parsedSelectedDataset.title;

            if (parsedSelectedDataset.files.length > 1) {
              formData.input = extractDir(parsedSelectedDataset.files[0].file_loc)
              formData.output = formData.input + "/Results";
            } else if(parsedSelectedDataset.files.length === 1) {
              formData.input = parsedSelectedDataset.files[0].file_loc;
              const directory = extractDir(formData.input)
              formData.output = directory + "/Results";
            }
            console.log(formData);

          // Verify the authenticity of the user
          isUserAuth(jwtToken)
          .then((authData) => {
            if (authData.isAuth) {
              formData.userID = authData.username;
              fetch(CELERY_BACKEND_API + '/tools/qc', {
                method: 'POST',
                headers: {
                  'Content-Type': 'application/json',
                },
                body: JSON.stringify(formData),
              })
              .then(response => {
                // Check the status code
                if (response.ok) {
                  console.log('Request succeeded');
                  console.log('Status code:', response.status);
                  return response.json();
                } else {
                  console.log('Request failed');
                  console.log('Status code:', response.status);
                  throw new Error('Error while making a call to the celery API');
                }
              })
              .then(response => {

                // Handle the successful response from the API
                console.log(JSON.stringify(response)); // Log the response data to the console

                // After a successfull task creation, store the intermediate task information in the database
                const taskId = response.task_id;
                const datasetId = parsedSelectedDataset.dataset_id;
                const method = formData.methods[0];
                const output = formData.output;
                      // Make API call to store the task information

                      const requestBody = {
                        taskId: taskId,
                        datasetId: datasetId,
                        method: method,
                        authToken:jwtToken,
                        outputPath: output
                      };
                      
                      fetch(`${SERVER_URL}/createTask`, {
                        method: 'POST',
                        headers: {
                          'Content-Type': 'application/json'
                        },
                        body: JSON.stringify(requestBody)
                      })
                      .then(response => {
                        if (response.ok && response.status === 201) {
                          // Response is successful (status code in the 200-299 range)
                            console.log('Task created successfully!');
                            setLoading(false);
                            setSuccessMessage('Form submitted successfully!');
                            setErrorMessage('');
                        } else if (response.status === 400) {
                          // Response is not successful
                          throw new Error('Please log in first');
                        }
                      })
                      .catch(error => {
                          if (error.message === 'Please log in first') {
                            navigate('/routing');
                            return;
                        } else {
                            console.error(error);
                        }
                        setLoading(false);
                        setSuccessMessage('');
                        setErrorMessage('An error occurred while submitting the form.');
                      });
                    setFormErrors("");
              })
              .catch(error => {
                // Handle any errors that occur during the API call
                console.error(error);
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

      const handleDatasetChange = event => {
        let value = event.target.value;
        if(value !== "") {
          setSelectedDataset(event.target.value);
        } else {
          setSelectedDataset([]);
        }
      };


  useEffect(() => {
    import(`./../../../schema/react-json-schema/Tools/${filterName}.json`)
    .then((module) => {
      setFilterSchema(module.default);
      console.log(filterSchema);
    })
    .catch((error) => {
      console.error('Error loading filter schema:', error);
      setFilterSchema(null);
    });
  }, [filterName, filterSchema]);

  return (
    <div className='tools-container common-class-tools-and-workflows'>
      <div className="separator heading">
        <div className="stripe"></div> 
          <h2 className="h-sm font-weight-bold">
            Tool Parameters 
          </h2> 
        <div className="stripe"></div>
      </div>
      {/* {formErrors && <span className="error">{formErrors}</span>} */}
      <div>
        <InputDataComponent handleDatasetChange={handleDatasetChange} formErrors={formErrors}/>
      </div>
            
        {filterSchema ? (
            <Form
            schema={filterSchema}
            formData={formData}
            widgets={widgets}
            onChange={({ formData }) => setFormData(formData)}
            uiSchema={uiSchema}
            onSubmit={onSubmit}
        />
          ) : (
            <div>No Schema for this tool.</div>
          )}

          {loading && (
            <div id="loadingIcon">
              {/* Replace with your loading icon */}
              <FontAwesomeIcon icon={faSpinner} spin />
              <p>Loading...</p>
            </div>
          )}
          {successMessage && (
            <div id="tooltip" className="success">
              {successMessage}
            </div>
          )}
          {errorMessage && (
            <div id="tooltip" className="error">
              {errorMessage}
            </div>
          )}
    </div>
  )
};
