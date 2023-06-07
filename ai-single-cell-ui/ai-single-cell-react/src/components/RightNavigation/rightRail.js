import React, { useState, useEffect } from "react";
import { FontAwesomeIcon } from "@fortawesome/react-fontawesome";
import {getCookie } from "../../utils/utilFunctions";
import { getStorageDetails } from '../../utils/utilFunctions';
import {faArrowRightArrowLeft , faPlus, faCaretDown, faPen, faDatabase, faLocationDot, faArrowsRotate, faSquareCheck, faCompress, faGear, faChevronDown, faCross, faXmark, faAngleDoubleDown, faInfoCircle} from "@fortawesome/free-solid-svg-icons";
import MyTasksSideNav from "../MyData/myTasksSideNav";

function RightRail() {
    const [isUserLoggedIn, setIsUserLoggedIn] = useState(false);
    const [usedStorage, setUsedStorage] = useState(0);
    const [totalStorage, setTotalStorage] = useState(5);

    useEffect(() => {
        const userCookie = getCookie('jwtToken');

        if(userCookie!== '') {
            setIsUserLoggedIn(true);
        }
        
      }, [isUserLoggedIn]);

      useEffect(() => {
        const jwtToken = getCookie('jwtToken'); // Assuming you have a function to get the JWT token
    
        getStorageDetails(jwtToken)
          .then(data => {
            setUsedStorage(data.usedStorage);
            setTotalStorage(data.totalStorage);
          });
      }, []);

  return (
    <div>
            {isUserLoggedIn && (
    <div className="right-container border-l border-gray-100 right-rail-container"> 
        <div className="rightpane-1">
            <div className="search-window">
                <div className="row1-search-window rows-search-window">
                    <div className="left-side-of-container">
                        <span className="rightpane-text rightpane-text-history">History</span>
                    </div>
                    <div className="right-side-of-container">
                        <FontAwesomeIcon icon={faPlus} title="Create new history"/>
                        <FontAwesomeIcon icon={faArrowRightArrowLeft} title="Switch to history"/>
                        <FontAwesomeIcon icon={faCaretDown} title="History options"/>
                    </div>
                </div>
                {/* <div className="row2-search-window rows-search-window">
                    <input
                        type="text"
                        placeholder="Search Datasets..."
                        className="w-full dark:bg-gray-950 pl-8 form-input-alt h-9 pr-3 focus:shadow-xl"
                    />
                    <FontAwesomeIcon className="right-rail-search-svg" icon={faAngleDoubleDown} />
                    <FontAwesomeIcon className="right-rail-search-svg" icon={faXmark} />
                </div> */}

<div class="row2-search-window rows-search-window"><div role="group" class="input-group"><input type="text" placeholder="search datasets" class="form-control form-control-sm" data-description="filter text input" fdprocessedid="zujxpl" id="__BVID__1228"/> <div class="input-group-append"><button data-description="show advanced filter toggle" aria-label="Show advanced filter" type="button" aria-pressed="false" autocomplete="off" class="btn btn-secondary btn-sm" fdprocessedid="n5md6r"><svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="angle-double-down" role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 320 512" class="svg-inline--fa fa-angle-double-down search-icon-svg"><path fill="currentColor" d="M143 256.3L7 120.3c-9.4-9.4-9.4-24.6 0-33.9l22.6-22.6c9.4-9.4 24.6-9.4 33.9 0l96.4 96.4 96.4-96.4c9.4-9.4 24.6-9.4 33.9 0L313 86.3c9.4 9.4 9.4 24.6 0 33.9l-136 136c-9.4 9.5-24.6 9.5-34 .1zm34 192l136-136c9.4-9.4 9.4-24.6 0-33.9l-22.6-22.6c-9.4-9.4-24.6-9.4-33.9 0L160 352.1l-96.4-96.4c-9.4-9.4-24.6-9.4-33.9 0L7 278.3c-9.4 9.4-9.4 24.6 0 33.9l136 136c9.4 9.5 24.6 9.5 34 .1z" class=""></path></svg></button> <button aria-label="Clear filters" data-description="clear filters" type="button" class="btn btn-secondary btn-sm" fdprocessedid="6mmb9"><svg aria-hidden="true" focusable="false" data-prefix="fas" data-icon="times" role="img" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 352 512" class="svg-inline--fa fa-times search-icon-svg"><path fill="currentColor" d="M242.72 256l100.07-100.07c12.28-12.28 12.28-32.19 0-44.48l-22.24-22.24c-12.28-12.28-32.19-12.28-44.48 0L176 189.28 75.93 89.21c-12.28-12.28-32.19-12.28-44.48 0L9.21 111.45c-12.28 12.28-12.28 32.19 0 44.48L109.28 256 9.21 356.07c-12.28 12.28-12.28 32.19 0 44.48l22.24 22.24c12.28 12.28 32.2 12.28 44.48 0L176 322.72l100.07 100.07c12.28 12.28 32.2 12.28 44.48 0l22.24-22.24c12.28-12.28 12.28-32.19 0-44.48L242.72 256z" class=""></path></svg></button></div></div></div>
                <div className="row3-search-window rows-search-window">
                    <div className="left-side-of-container">
                        <span className="rightpane-text rightpane-text-unnamed-history">Unnamed history</span>
                    </div>
                <div className="right-side-of-container">
                        <FontAwesomeIcon icon={faPen} title="Edit"/>
                </div>
                </div>
            </div>
            <div className="row4-search-window">
                    <div className="left-side-of-container">
                        <FontAwesomeIcon icon={faDatabase} title="History size"/><p className="storage-span">{usedStorage} GB/{totalStorage} GB</p>
                    </div>
                    <div className="right-side-of-container">
                        <FontAwesomeIcon icon={faLocationDot} title="Show active"/>
                        <FontAwesomeIcon icon={faArrowsRotate} title="Last Refreshed" />
                    </div>
                </div>
                <div className="row5-search-window border-b">
                    <div className="left-side-of-container">
                        <FontAwesomeIcon icon={faSquareCheck} />
                        <FontAwesomeIcon icon={faCompress} />
                    </div>
                    <div className="right-side-of-container">
                        <FontAwesomeIcon icon={faGear} />
                    </div>
                </div>
        </div>
        <div className="results-window">
        <MyTasksSideNav/>
        <div><div><div role="alert" aria-live="polite" aria-atomic="true" class="alert m-2 alert-info"><h4 class="mb-1"><FontAwesomeIcon icon={faInfoCircle} /><span>This history is empty.</span></h4> <p><a href="#">You can load your own data</a> <span>or</span> <a href="#">get data from an external source</a>.
            </p>
        </div>
        </div>
        </div>

        </div>
        </div>
  )
  }
    </div>
  );
}

export default RightRail;
