import React from "react";
  
const ToggleSwitch = ({ label , toggleSwitchForPublicDatasets}) => {
  return (
    <div>
      {label}{" "}
      <div className="toggle-switch">
        <input type="checkbox" className="checkbox" 
               name={label} id={label} onClick={toggleSwitchForPublicDatasets} />
        <label className="label" htmlFor={label}>
          <span className="inner" />
          <span className="switch" />
        </label>
      </div>
    </div>
  );
};
  
export default ToggleSwitch;