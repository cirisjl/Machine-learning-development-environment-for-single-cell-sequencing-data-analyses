import React, { useState, useEffect } from 'react';
// import Select from 'react-select';
import Form from 'react-jsonschema-form';
import { SERVER_URL } from '../../../constants/declarations';
import createUISchema from '../../../schema/UI-schema/basicFormUISchema';
import formSchema from '../../../schema/react-json-schema/basicFormSchema.json'


const BasicFormComponent = () => {
    const [formData, setFormData] = useState({});
    const [commonOptions, setCommonOptions] = useState({});
    const [uiSchema, setUiSchema] = useState(null); // Initialize uiSchema as null

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
        
          // Create the uiSchema based on commonOptions
          const generatedUiSchema = createUISchema(data);
          setUiSchema(generatedUiSchema);

        } catch (error) {
          console.error(`Error fetching common options:`, error);
        }
      };
  
      fetchCommonOptions();
    }, []);
  
  
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
