import React, { useState, useEffect } from 'react';
import { getCookie} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
// import schema from '../../../react-json-schema/Tools/normalizeUsingScanpySchema.json';
import Form from 'react-jsonschema-form';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import InputDataComponent from './inputDataCollection';

export default function ToolsDetailsComponent(props) {
    const filterName = props.filter;

    console.log(filterName);

    let jwtToken = getCookie('jwtToken');
    const [formData, setFormData] = useState({});
    const [filterSchema, setFilterSchema] = useState(null);
    const [selectedDataset, setSelectedDataset] = useState([]);
    const [formErrors, setFormErrors] = useState("");

    const navigate = useNavigate();

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
            const parsedSelectedDataset = JSON.parse(selectedDataset);
            formData.dataset = parsedSelectedDataset.title;

            if (parsedSelectedDataset.files.length > 1) {
              const fileLocParts = parsedSelectedDataset.files[0].file_loc.split('/');
              fileLocParts.shift(); // Remove the first empty element
              formData.input = '/' + fileLocParts[0];
            } else if(parsedSelectedDataset.files.length === 1) {
              formData.input = parsedSelectedDataset.files[0].file_loc;
            }
            console.log(formData);
            setFormErrors("");
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
    </div>
  )
};
