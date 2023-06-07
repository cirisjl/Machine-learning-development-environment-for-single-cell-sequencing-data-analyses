import LeftNav from "../../../components/LeftNavigation/leftNav"
import { NormalizeUsingScanpyComponent } from "../../../components/MyData/Tools/normalizeWithScanpy"
import ToolsDetailsComponent from "../../../components/MyData/Tools/toolsDetailsComponent";
import RightRail from "../../../components/RightNavigation/rightRail"
import React, { useState } from 'react';


export default function NormalizeUsingScanpy() {

    const [selectedFilter, setSelectedFilter] = useState(null);
  
    const handleFilterSelection = (filterName) => {
      setSelectedFilter(filterName);
  
    //   // Dynamically import the schema file based on the filter name
    //   import(`./schemas/${filterName}.json`)
    //     .then((module) => {
    //       setFilterSchema(module.default);
    //     })
    //     .catch((error) => {
    //       console.error('Error loading filter schema:', error);
    //     });
    };
  
    return(
        <div className="page-container">
            <div className="left-nav border-r left-nav-background">
                <LeftNav handleFilterSelection={handleFilterSelection}/>
            </div>
      {/* Render the selected filter details in the middle of the page */}
      {selectedFilter && (
        <div className="filter-details-tools main-content">
          <ToolsDetailsComponent filter={selectedFilter} />
        </div>
      )}
       {!selectedFilter && (
            <div className="tool-message">
            <p>Please select a tool to run</p>
            </div>
      )}
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}