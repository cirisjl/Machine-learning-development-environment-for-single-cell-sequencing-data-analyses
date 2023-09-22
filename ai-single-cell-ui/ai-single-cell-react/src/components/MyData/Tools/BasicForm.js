import React, { useState, useEffect } from 'react';
import { getCookie, isUserAuth} from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import Form from 'react-jsonschema-form';
import Toggle from 'react-toggle';
import 'react-toggle/style.css';
import InputDataComponent from './inputDataCollection';
import { CELERY_BACKEND_API, SERVER_URL } from '../../../constants/declarations';
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import { faSpinner } from '@fortawesome/free-solid-svg-icons';
import { uiSchema } from '../../../schema/UI-schema/Tools/evaluation/basicFormUISchema';
import schema from '../../../schema/react-json-schema/Tools/evaluation/basicFormSchema.json'


export default function BasicFormComponent(props) {
    const [formData, setFormData] = useState({});

    const onSubmit = ({ formData }) => {

        console.log(formData);
      };

  return (
    <div className='tools-container common-class-tools-and-workflows'>
            
        {schema && uiSchema ? (
            <Form
            schema={schema}
            formData={formData}
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
