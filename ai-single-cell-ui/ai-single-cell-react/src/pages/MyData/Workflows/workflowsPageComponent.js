import LeftNavComponent from "../../../components/MyData/Workflows/components/leftnavComponent";
import { WorkflowsComponent } from "../../../components/MyData/Workflows/workflowsComponent"
import RightRail from "../../../components/RightNavigation/rightRail"
import useUserAuthCheck from "../../../components/common_components/userAuthCheckComponent";
import React, {useState} from 'react'

export default function WorkflowsPageComponent() {

    const [selectedWorkflow, setSelectedWorkflow] = useState("");
    const [uniqueFilter, setUniqueFilter] = useState("");

    const handleFilterSelection = (category, filter) => {
        setUniqueFilter(category+ "_" + filter);
        setSelectedWorkflow(category);
        console.log(category + filter);
    };


    useUserAuthCheck();

    return(
        <div className="page-container">
            <div className="left-nav">
                <LeftNavComponent selectedWorkflow={selectedWorkflow} setSelectedWorkflow={setSelectedWorkflow} handleFilterSelection={handleFilterSelection}/>
            </div>
            <div className="main-content">
                <WorkflowsComponent selectedWorkflow={selectedWorkflow} uniqueFilter={uniqueFilter}/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}