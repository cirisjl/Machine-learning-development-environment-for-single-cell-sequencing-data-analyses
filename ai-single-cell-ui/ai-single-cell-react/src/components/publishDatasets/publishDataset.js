import RightRail from "../RightNavigation/rightRail";
import LeftNav from "./components/leftNav";
import MiddleContent from "./components/mainContent";
import React, { useState } from 'react';
import './publishDatasets.css'; // Import a CSS file for styles

export default function PublishDataset() {

    const [taskStatus, setTaskStatus] = useState({
        1: false, // Task 1 is initially not completed
        2: false,
        3: false,
        4: false,
        5: false,
        6: false,
        7: false
        // Add other tasks here
      });
    

    const [activeTask, setActiveTask] = useState(1); // Initialize with the first task

    return(
        <div className="page-container">
            <div className="left-nav">
            <LeftNav activeTask={activeTask} setActiveTask={setActiveTask} taskStatus={taskStatus} />
            </div>
            <div className="main-content">
                <MiddleContent activeTask={activeTask} setTaskStatus={setTaskStatus} />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}