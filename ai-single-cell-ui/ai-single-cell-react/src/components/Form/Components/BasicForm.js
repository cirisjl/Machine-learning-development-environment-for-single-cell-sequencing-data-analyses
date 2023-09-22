import React, { useState, useEffect } from 'react';
// import Select from 'react-select';
import Form from 'react-jsonschema-form';
import { SERVER_URL } from '../../../constants/declarations';
// import createUISchema from '../../../schema/UI-schema/basicFormUISchema';
import formSchema from '../../../schema/react-json-schema/basicFormSchema.json'
import MyCreatableSelect from '../Components/Creatable'


const BasicFormComponent = () => {
    const [formData, setFormData] = useState({});
    const [commonOptions, setCommonOptions] = useState({});
    // const [uiSchema, setUiSchema] = useState(null); // Initialize uiSchema as null

    useEffect(() => {
     console.log(formData)
      }, [formData]);

    const onSubmit = ({ formData }) => {

        console.log(formData);
      };

    useEffect(() => {
      // Fetch common options from MongoDB API once when the component mounts
      const fetchCommonOptions = async () => {
        try {
          const optionsURL = `${SERVER_URL}/mongoDB/api/options`;
          const response = await fetch(optionsURL);
          const data = await response.json();
          setCommonOptions(data);
        
          // // Create the uiSchema based on commonOptions
          // const generatedUiSchema = createUISchema(data, formData, setFormData);
          // setUiSchema(generatedUiSchema);

        } catch (error) {
          console.error(`Error fetching common options:`, error);
        }
      };
  
      fetchCommonOptions();
    }, []);

    const uiSchema = {
      "Dataset": {
        "ui:placeholder": "Enter the Dataset name"
      },
      "Task": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Task" options={commonOptions}/>
          </div>
        ),
      },
      "Downloads": {
        "ui:placeholder": "http://"
      },
      "Title": {
        "ui:placeholder": "Enter title"
      },
      "Author": {
          "ui:placeholder": "Select/Create an Option",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <MyCreatableSelect fieldName="Author" options={commonOptions} />
            </div>
          ),
      },
      "Reference (paper)": {
        "ui:placeholder": "Enter Reference"
      },
      "Abstract": {
        "ui:widget": "textarea"
      },
      "DOI": {
        "ui:placeholder": "http://"
      },
      "Species": {
          "ui:placeholder": "Select/Create an Option",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <MyCreatableSelect fieldName="Species" options={commonOptions} />
            </div>
          ),
      },
      "Sample Type": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Sample Type" options={commonOptions}/>
          </div>
        ),
      },
      "Anatomical Entity": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Anatomical Entity" options={commonOptions}/>
          </div>
        ),
      },
      "Organ Part": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Organ Part" options={commonOptions}/>
          </div>
        ),
      },
      "Model Organ": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Model Organ" options={commonOptions}/>
          </div>
        ),
      },
      "Selected Cell Types": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Selected Cell Types" options={commonOptions}/>
          </div>
        ),
      },
      
      "Library Construction Method": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Library Construction Method" options={commonOptions}/>
          </div>
        ),
      },
      
      "Nucleic Acid Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Nucleic Acid Source" options={commonOptions}/>
          </div>
        ),
      },
  
      
      "Paired End": {
        "ui:widget": "radio" // could also be "select"
      },
      
      "Analysis Protocol": {
        "ui:placeholder": "Enter ...",
      },
      "Disease Status (Specimen)": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Specimen)" options={commonOptions} />
          </div>
        ),
      },
      
      "Disease Status (Donor)": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Donor)" options={commonOptions} />
          </div>
        ),
      },
      
      "Development Stage": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Development Stage" options={commonOptions}/>
          </div>
        ),
      },
      "Donor Count": {
        "ui:placeholder": "Enter ...",
      },
      "Cell Count Estimate": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Cell Count Estimate" options={commonOptions}/>
          </div>
        ),
      },
      "Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Source" options={commonOptions}/>
          </div>
        ),
      },
      "Source Key": {
        "ui:placeholder": "Enter ...",
      },
      "Submission Date": {
        "className":"date"
      },
  };
  
  
    return (
      <div className='tools-container common-class-tools-and-workflows'>
            
      {formSchema && uiSchema ? (
          <Form
          schema={formSchema}
          formData={formData}
          onChange={({ formData }) => setFormData(formData)}
          uiSchema={uiSchema}
          onSubmit={onSubmit}
      />
        ) : (
          <div>No Schema for this tool.</div>
        )}
  </div>
    );
};

export default BasicFormComponent;
