import RightRail from "../RightNavigation/rightRail";
import LeftNav from "./components/leftNav";
import MiddleContent from "./components/mainContent";
import React, { useState } from 'react';

export default function PublishDataset() {

    const [activeTask, setActiveTask] = useState(1); // Initialize with the first task

    return(
        <div className="page-container">
            <div className="left-nav">
                <LeftNav activeTask={activeTask} setActiveTask={setActiveTask} />
            </div>
            <div className="main-content">
                <MiddleContent activeTask={activeTask} />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}