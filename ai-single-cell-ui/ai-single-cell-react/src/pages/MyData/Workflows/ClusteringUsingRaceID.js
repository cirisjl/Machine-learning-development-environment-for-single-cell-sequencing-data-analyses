import { ClusteringUsingRaceIDComponent } from "../../../components/MyData/Workflows/clusteringUsingRaceID"
import RightRail from "../../../components/RightNavigation/rightRail"

export default function ClusteringUsingRaceID() {
    return(
        <div className="page-container">
            <div className="left-nav">
            </div>
            <div className="main-content">
                <ClusteringUsingRaceIDComponent />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}