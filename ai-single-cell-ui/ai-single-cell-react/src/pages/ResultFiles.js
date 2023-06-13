// import leftNav from "../components/leftNav";
import IntermediateFiles from "../components/MyData/intermediatefiles";
import RightRail from "../components/RightNavigation/rightRail";
import { useLocation } from 'react-router-dom';

export default function ResultFiles() {
    const location = useLocation();
    const taskId = location.state?.taskId || '';
    const results_path= location.state?.results_path || '';
    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <IntermediateFiles taskId={taskId} results_path={results_path}/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}