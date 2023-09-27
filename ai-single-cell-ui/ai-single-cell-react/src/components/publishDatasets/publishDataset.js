import RightRail from "../components/RightNavigation/rightRail";
import LeftNav from "./components/leftNav";

export default function PublishDataset() {
    return(
        <div className="page-container">
            <div className="left-nav">
                <LeftNav />
            </div>
            <div className="main-content">
                
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}