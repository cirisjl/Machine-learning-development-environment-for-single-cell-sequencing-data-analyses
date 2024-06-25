// import leftNav from "../components/leftNav";
import TaskTable from "../components/MyData/taskTable";
import RightRail from "../components/RightNavigation/rightRail";

export default function MyTasks() {
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