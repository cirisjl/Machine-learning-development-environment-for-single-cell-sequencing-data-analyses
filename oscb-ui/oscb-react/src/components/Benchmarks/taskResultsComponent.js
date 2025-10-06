import React, { useState, useEffect } from 'react';
import ReactMarkdown from "react-markdown"
import 'github-markdown-css';
import axios from 'axios';
import gfm from "remark-gfm";
import remarkImgToJsx from "remark-unwrap-images";
// import {NODE_API_URL} from '../../constants/declarations'
// import { getCookie } from "../../utils/utilFunctions";
// import { useNavigate } from 'react-router-dom';
import { DIRECTUS_URL } from '../../constants/declarations'
import RightRail from '../RightNavigation/rightRail';
import SearchTasks from './components/taskResults';


export default function TaskResultsComponent(task_type) {
    const [markdownText, setMarkdownText] = useState('');
    const title = task_type.task_type;
    // const navigate = useNavigate();
    // let jwtToken = getCookie('jwtToken');

    // useEffect(() => {
    //     let jwtToken = getCookie('jwtToken');
    //     if(jwtToken===undefined || jwtToken === '') {
    //         navigate('/routing');
    //     }
    // },[]);
    useEffect(() => {
            async function fetchFileData() {
            try {
                const response = await axios.get(DIRECTUS_URL + "/items/filemappings?filter[filename]=" + title.replace(" ", "_"));
                const data = response.data.data;
                
                if(data.length === 1) {
                    const fileMappingObject = data[0];
                    const fileID = fileMappingObject.fileID;
                    if(fileID !== null) {
                        fetch(DIRECTUS_URL + "/assets/" + fileID)
                        .then(response => response.text())
                        .then(data => setMarkdownText(data))
                        .catch(error => console.error('Error retrieving markdown:', error));
                    }
                }
                
              } catch (error) {
                console.error('Error retrieving data:', error);
              }
            }
            
            fetchFileData();
          }, []);

    return (
        <div className="task-results-container eighty-twenty-grid">
            <div className="main-content task-builder-task">
                <h1 style={{ textAlign: "left" }}>{title}</h1>
                <h2>Datasets</h2>
                <p><SearchTasks taskType={task_type} /></p>
                <hr/>
                <p><ReactMarkdown plugins={[gfm, remarkImgToJsx]} children={markdownText} /></p>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    );
}

// export default TaskResultsComponent;
