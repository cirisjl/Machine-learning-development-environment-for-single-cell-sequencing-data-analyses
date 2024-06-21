// import LeftNav from "../components/leftNav";
import DatasetSelectionDialog from "../../components/publishDatasets/components/datasetsDialog";
import RightRail from "../../components/RightNavigation/rightRail";
import StorageChart from "../../components/MyData/storageChart";
import { getCookie } from "../../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';
import { useState, useEffect } from 'react';

export default function MyData() {

    const navigate = useNavigate();

    const [selectedFilter, setSelectedFilter] = useState(null);
    const [category, setCategory] = useState(null);
    const handleFilterSelection = (category, filter) => {
        setSelectedFilter(category + "_" + filter);
        setCategory(category);
    };

    useEffect(() => {
        let jwtToken = getCookie('jwtToken');
        if (jwtToken === undefined || jwtToken === '') {
            navigate('/routing');
        }
    }, []);

    return (
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <StorageChart />
                    <DatasetSelectionDialog
                        onSelect={null}
                        multiple={true}
                        onClose={null}
                        isVisible={true}
                        selectedDatasets={null}
                        fromToolsPage={true}
                    />
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}