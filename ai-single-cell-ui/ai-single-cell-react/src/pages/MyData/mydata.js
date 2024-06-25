// import LeftNav from "../components/leftNav";
import RightRail from "../../components/RightNavigation/rightRail";
import StorageChart from "../../components/MyData/storageChart";
import { getCookie } from "../../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';
import { useEffect } from 'react';

export default function MyData() {

    const navigate = useNavigate();

    useEffect(() => {
        let jwtToken = getCookie('jwtToken');
        if(jwtToken===undefined || jwtToken === '') {
            navigate('/routing');
        }
    },[]);

    return(
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
            <StorageChart/>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}