import React, { useState, useEffect } from 'react';
// import Select from 'react-select';
import Form from 'react-jsonschema-form';
import formSchema from "../components/Form/basicFormSchema.json"
import {uiSchema} from "../schema/UI-schema/Tools/evaluation/basicFormUISchema.js"
import { SERVER_URL } from '../constants/declarations';


const DynamicForm = () => {
    const [formData, setFormData] = useState({});
    const [commonOptions, setCommonOptions] = useState({});
    // const [uiSchema, setUiSchema] = useState(null); // Initialize uiSchema as null

    // useEffect(() => {
    //   // Fetch common options from MongoDB API once when the component mounts
    //   const fetchCommonOptions = async () => {
    //     try {
    //       const optionsURL = `${SERVER_URL}/mongoDB/api/options`;
    //       const response = await fetch(optionsURL);
    //       const data = await response.json();
    //       setCommonOptions(data);
        
    //       // Create the uiSchema based on commonOptions
    //       // const generatedUiSchema = createUISchema(data);
    //       // setUiSchema(uiSchema);

    //     } catch (error) {
    //       console.error(`Error fetching common options:`, error);
    //     }
    //   };
  
    //   fetchCommonOptions();
    // }, []);
  
  
    return (
      <div>
        {uiSchema && (
          <Form
            schema={formSchema}
            onChange={({ formData }) => setFormData(formData)}
            uiSchema={uiSchema}
          />
        )}
      </div>
    );
};

export default DynamicForm;
