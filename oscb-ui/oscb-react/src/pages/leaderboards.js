// import LeftNav from "../components/leftNav";
import RightRail from "../components/RightNavigation/rightRail";
import { getCookie } from "../utils/utilFunctions";


export default function Leaderboards() {
    let jwtToken = getCookie('jwtToken');

    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <h1>I'm inside Leaderboards page</h1>
            </div>
            {(jwtToken != undefined && jwtToken != '') && (<div className="right-rail">
                <RightRail />
            </div>)}
        </div>
    )
}