// import leftNav from "../components/leftNav";
import IntermediateFiles from "../components/MyData/intermediatefiles";
import RightRail from "../components/RightNavigation/rightRail";
import { useLocation } from 'react-router-dom';

export default function ResultFiles() {
    const location = useLocation();
    const taskId = location.state?.taskId || '';
    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <IntermediateFiles taskId={taskId}/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}