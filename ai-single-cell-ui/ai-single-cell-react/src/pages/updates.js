import React, { useState, useEffect } from "react";
import axios from "axios";
import LeftNav from "../components/LeftNavigation/leftNav";
import RightRail from "../components/RightNavigation/rightRail";

const UPDATES_PAGE_API = `http://${process.env.REACT_APP_HOST_URL}:8055`

export default function Updates() {

    const [updates,setUpdates] = useState([]);

    const fetchData = async () => {
        try {
          const response = await axios.get(UPDATES_PAGE_API + '/items/updates');
          console.log(response.data.data);
          const sortedData = response.data.data.sort((a, b) => {
            return new Date(b.date_published) - new Date(a.date_published);
          });
          
          setUpdates(sortedData);
        } catch (error) {
          console.error(error);
        }
      }

    useEffect(() => {
        fetchData();
    },[])
    
    function formatDate(date_given) {
        const months = ['January', 'February','March','April','May','June','July','August','September','October','November','December']
        const splitDate = date_given.split("-");
        const month = months[parseInt(splitDate[1]) - 1]
        const date  = splitDate[2]
        const year = splitDate[0]
        let formattedDate = month.concat(" ", date, ", ",year)
        return formattedDate;
    }

    return(
        <div className="updates-container">
            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>

            <div className="main-content">
                <h1>Updates from OSGB</h1>
                {/* <h2>{user[user.length-1].Title} - {formatDate(user[user.length-1].Date_published)}</h2> */}

                {updates.map((item,index)=>(
                    <div className="updates-div border-gray-100">
                        <h4>{formatDate(item.date_published)} : {item.title}</h4>
                        <div dangerouslySetInnerHTML={{ __html: item.update_description }}></div>
                    </div>

                ))}
            </div>
            
            <div className="right-rail">
                <RightRail />
            </div>
        </div>
    )
}