import { faFile, faFolder } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React, { useState, useEffect } from 'react';
import { getCookie} from '../../../utils/utilFunctions';
import { PREVIEW_DATASETS_API } from '../../../constants/declarations';
import axios from 'axios';
import { useNavigate } from 'react-router-dom';
import { FormHelperText } from '@material-ui/core';


export default function InputDataComponent(props) {
    let jwtToken = getCookie('jwtToken');
    const [datasets, setDatasets] = useState([]);
    const navigate = useNavigate();

    const handleDatasetChange = props.handleDatasetChange;

    useEffect(() => {
        if (jwtToken) {
          const fetchData = async () => {
            const response = await axios.get(PREVIEW_DATASETS_API + "/preview/datasets?authToken=" + jwtToken);
            setDatasets(Object.values(response.data));
            console.log(Object.values(response.data));
          };
          fetchData();
        } else {
          navigate("/routing");
        }
      }, [jwtToken, navigate]);

  return (
    <div className='inputData-management'>
        <div id="upload-data-div">
            <div className='common-row-wrap'>
                <b>Choose the input dataset</b> <span data-v-22825496="" class="ui-form-title-message warning"> * required </span>
                <br/>
            </div>
            <div className='wrap-div-row common-row-wrap'>
                <div className='icons-wrap-input'>
                <button type="button" fdprocessedid="s89ftr">
                    <FontAwesomeIcon icon={faFile} />
                </button>

                <button type="button" fdprocessedid="s89ftr">
                    <FontAwesomeIcon icon={faFile} />
                </button>

                <button type="button" fdprocessedid="s89ftr">
                    <FontAwesomeIcon icon={faFolder} />
                </button>
                </div>
                <div className='datasets-input-select'>
                    <select onChange={handleDatasetChange}>
                        <option value="">Select the dataset from the list</option>
                        {datasets.map((dataset) => (
                            <option value={JSON.stringify(dataset)}>{dataset.title}</option>
                        ))}
                    </select>
                    {props.formErrors && <FormHelperText>{props.formErrors}</FormHelperText>}
                </div>
                <div>
                    <button type="button" fdprocessedid="s89ftr" className='input-path-button' title='Create Dataset'>
                        <svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="folder-open" class="svg-inline--fa fa-folder-open " role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 576 512">
                            <path fill="currentColor" d="M88.7 223.8L0 375.8V96C0 60.7 28.7 32 64 32H181.5c17 0 33.3 6.7 45.3 18.7l26.5 26.5c12 12 28.3 18.7 45.3 18.7H416c35.3 0 64 28.7 64 64v32H144c-22.8 0-43.8 12.1-55.3 31.8zm27.6 16.1C122.1 230 132.6 224 144 224H544c11.5 0 22 6.1 27.7 16.1s5.7 22.2-.1 32.1l-112 192C453.9 474 443.4 480 432 480H32c-11.5 0-22-6.1-27.7-16.1s-5.7-22.2 .1-32.1l112-192z"></path>
                        </svg>
                    </button>
                </div>
            </div>
        </div>
    </div>
  )
};
