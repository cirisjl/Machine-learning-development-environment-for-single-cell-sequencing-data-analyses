import React, { useState, useEffect } from 'react';
import axios from 'axios';
import {SERVER_URL} from '../../constants/declarations'
import { getCookie } from "../../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';
import RightRail from '../RightNavigation/rightRail';
import SearchTasks from './components/taskResults';

function TaskResultsComponent(task_type) {

    const navigate = useNavigate();

    useEffect(() => {
        let jwtToken = getCookie('jwtToken');
        if(jwtToken===undefined || jwtToken === '') {
            navigate('/routing');
        }
    },[]);

    return (
        <div className="task-results-container eighty-twenty-grid">
            <div className="main-content task-builder-task">
                <SearchTasks taskType={task_type} />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    );
}

export default TaskResultsComponent;
