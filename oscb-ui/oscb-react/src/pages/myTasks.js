// import leftNav from "../components/leftNav";
import React, { useState, useEffect } from 'react';
import TaskTable from "../components/MyData/taskTable";
import RightRail from "../components/RightNavigation/rightRail";
import { getCookie } from "../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';


export default function MyTasks() {
    const navigate = useNavigate();
    let jwtToken = getCookie('jwtToken');

    useEffect(() => {
        let jwtToken = getCookie('jwtToken');
        if(jwtToken===undefined || jwtToken === '') {
            navigate('/routing');
        }
    },[]);

    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <TaskTable/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}