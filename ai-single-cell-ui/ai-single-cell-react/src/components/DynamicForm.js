import React, { useState, useEffect } from 'react';
// import Select from 'react-select';
import Form from 'react-jsonschema-form';
import formSchema from "../components/Form/basicFormSchema.json"
import createUISchema from "../components/Form/basicFormUISchema.js"
import { SERVER_URL } from '../constants/declarations';


const DynamicForm = () => {
    const [formData, setFormData] = useState({});
    const [commonOptions, setCommonOptions] = useState({});
    const [UISchema, setUISchema] = useState(null);

    // const uiSchema = createUISchema(commonOptions); // Pass commonOptions to create the uiSchema

    useEffect(() => {
      // Fetch common options from MongoDB API once when the component mounts
      const fetchCommonOptions = async () => {
        try {
          const optionsURL = `${SERVER_URL}/mongoDB/api/options`;
          const response = await fetch(optionsURL);
          const data = await response.json();
          setCommonOptions(data);
          setUISchema(createUISchema(commonOptions))
        } catch (error) {
          console.error(`Error fetching common options:`, error);
        }
      };
  
      fetchCommonOptions();
    }, []);

  return (
    <div className='tools-container common-class-tools-and-workflows'>
         <Form
            schema={formSchema}
            // formData={formData}
            // widgets={widgets}
            onChange={({ formData }) => setFormData(formData)}
            uiSchema={UISchema}
            // onSubmit={onSubmit}
        />
    </div>
  );
};

export default DynamicForm;
