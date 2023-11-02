// import LeftNav from "../components/leftNav";
import RightRail from "../components/RightNavigation/rightRail";

export default function Docs() {
    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <h1>I'm inside Docs page</h1>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}