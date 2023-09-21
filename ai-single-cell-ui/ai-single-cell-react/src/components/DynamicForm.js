import React, { useState } from 'react';
// import Select from 'react-select';
import Form from 'react-jsonschema-form';
import formSchema from "../components/Form/basicFormSchema.json"
import uiSchema from "../components/Form/basicFormUISchema.js"

const DynamicForm = () => {
    const [formData, setFormData] = useState({});


  return (
    <div className='tools-container common-class-tools-and-workflows'>
         <Form
            schema={formSchema}
            // formData={formData}
            // widgets={widgets}
            onChange={({ formData }) => setFormData(formData)}
            uiSchema={uiSchema}
            // onSubmit={onSubmit}
        />
    </div>
  );
};

export default DynamicForm;
