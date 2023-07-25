import React, { useState, useEffect } from 'react';
import ReactDOM from 'react-dom';
import { getCookie, isUserAuth} from '../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';


export default function FlaskDashboard  () {
  const [dashApp, setDashApp] = useState(null);
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
        });
        const FLASK_BACKEND_API = `http://${process.env.REACT_APP_HOST_URL}:5003/dashboard?${queryParams}`

    fetch(FLASK_BACKEND_API) 
      .then(response => response.text())
      .then(html => setDashApp(html))
      .catch(error => console.error('Error fetching Dash app:', error));  

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
      {dashApp && (
        <div dangerouslySetInnerHTML={{ __html: dashApp }} />
      )}
    </div>
  );
};