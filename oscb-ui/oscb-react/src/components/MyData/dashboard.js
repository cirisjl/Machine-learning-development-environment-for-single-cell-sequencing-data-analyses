import React, { useState, useEffect } from 'react';
import ReactDOM from 'react-dom';
import { getCookie, isUserAuth} from '../../utils/utilFunctions';
import { useNavigate, useLocation } from 'react-router-dom';
import Iframe from 'react-iframe';
import { FLASK_BACKEND_API } from '../../constants/declarations'


export default function FlaskDashboard (props) {
  const [dashApp, setDashApp] = useState(null);
  const [flaskURL, setFlaskURL] = useState(null);

  const location = useLocation();
  const state = location.state; // This will contain the state passed through navigate

  // Now you can use the 'state' object to access the passed data
  const message = state.message;
  
  const navigate = useNavigate();

  useEffect(() => {

    const jwtToken = getCookie('jwtToken');

    isUserAuth(jwtToken)
    .then((authData) => {
      if (authData.isAuth) {
        const userID = authData.username;

        // Construct the query parameters
        const queryParams = new URLSearchParams({
            authToken: jwtToken,
            username: userID,
            title: state.title
        });

        setFlaskURL(FLASK_BACKEND_API)
      } else {
        console.warn("Unauthorized - please login first to continue");
        navigate("/routing");
      }
    })
    .catch((error) => {
      console.error(error);
    } 
    );
}, []);

  return (
    <div>
      {/* {dashApp && (
        <div dangerouslySetInnerHTML={{ __html: dashApp }} />
      )} */}
{flaskURL && (
    <Iframe url= {flaskURL} // Replace this with the Flask app URL
        width="100%"
        height="100vh"
        id="dashFrame"
        display="initial"
        position="relative"
    />
    )}
    </div>
  );
};