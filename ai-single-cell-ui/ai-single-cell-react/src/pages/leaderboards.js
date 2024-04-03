// import LeftNav from "../components/leftNav";
import LeaderCharts from "../components/LeftNavigation/leaderChart";
import RightRail from "../components/RightNavigation/rightRail";

export default function Leaderboards() {
    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <LeaderCharts />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}