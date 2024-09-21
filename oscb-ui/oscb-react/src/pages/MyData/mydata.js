// import LeftNav from "../components/leftNav";
import DatasetTable from "../../components/MyData/datasetTable";
import RightRail from "../../components/RightNavigation/rightRail";
import StorageChart from "../../components/MyData/storageChart";
import { getCookie } from "../../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';
import { useState, useEffect } from 'react';

export default function MyData() {

    const navigate = useNavigate();
    const filterCategory = null;
    const [shouldHideForSeurat, setShouldHideForSeurat] = useState(false);
    const [selectedDatasets, setSelectedDatasets] = useState({});

    useEffect(() => {
        let jwtToken = getCookie('jwtToken');
        if (jwtToken === undefined || jwtToken === '') {
            navigate('/routing');
        }
    }, []);

    const onSelectDataset = (dataset) => {
        let datasetId = dataset.Id;
        let currentSelectedDatasets = { ...selectedDatasets };

        if (currentSelectedDatasets[datasetId]) {
            delete currentSelectedDatasets[datasetId];
        } else {
            if (filterCategory !== "integration") {
                currentSelectedDatasets = {};
            }
            currentSelectedDatasets[datasetId] = dataset;
        }
        if (filterCategory === "quality_control") {
            // Check if any of the selected datasets should trigger hiding for Seurat
            const shouldHideForSeurat = Object.values(currentSelectedDatasets).some(dataset =>
                dataset.inputFiles.length === 1 &&
                (dataset.inputFiles[0].toLowerCase().endsWith('h5seurat') ||
                    dataset.inputFiles[0].toLowerCase().endsWith('rds') ||
                    dataset.inputFiles[0].toLowerCase().endsWith('robj'))
            );
            setShouldHideForSeurat(shouldHideForSeurat);
        }
        setSelectedDatasets(currentSelectedDatasets)
    };

  // Function to handle selection of sub-items
  const onSelectSubItem = (mainItem, subItem) => {
    const mainItemId = mainItem.Id;
    let currentSelectedDatasets = { ...selectedDatasets };
  
    // Check if the main item is already selected
    if (currentSelectedDatasets[mainItemId]) {
        // If sub-item is already selected, deselect it
        if (currentSelectedDatasets[mainItemId].selectedSubItem?.process_id  === subItem.process_id ) {
            delete currentSelectedDatasets[mainItemId];
        } else {
            // Update the selected main item with the selected sub-item
            currentSelectedDatasets[mainItemId] = {
                ...mainItem,
                selectedSubItem: subItem
            };
        }
    } else {
        // Select the main item and the sub-item
        currentSelectedDatasets = {
            [mainItemId]: {
                ...mainItem,
                selectedSubItem: subItem
            }
        };
    }
  
    setSelectedDatasets(currentSelectedDatasets);
  };

    return (
        <div className="page-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>
            <div className="main-content">
                <StorageChart />
                <div className="task-builder-task">
                    <DatasetTable
                        onSelect={onSelectDataset}
                        onClose={null}
                        isVisible={true}
                        selectedDatasets={selectedDatasets}
                        fromToolsPage={true}
                        onSelectSubItem = {onSelectSubItem}
                    />
                </div>
            </div>
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}