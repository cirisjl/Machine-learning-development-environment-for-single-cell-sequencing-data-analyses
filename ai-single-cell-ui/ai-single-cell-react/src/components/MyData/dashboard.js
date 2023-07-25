import React, { useState, useEffect } from 'react';
import ReactDOM from 'react-dom';
import DashHtmlComponents from 'dash-html-components';
import { FLASK_BACKEND_API} from '../../constants/declarations';
import { getCookie} from '../../utils/utilFunctions';


export const FlaskDashboard = () => {
  const [dashApp, setDashApp] = useState(null);

  useEffect(() => {

    const jwtToken = getCookie('jwtToken');

    isUserAuth(jwtToken)
    .then((authData) => {
      if (authData.isAuth) {
        userID = authData.username;

        const queryParams = {
            authToken: jwtToken,
            username: userID,
        };

    fetch(FLASK_BACKEND_API + "/dashboard", {
        method: 'GET',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(queryParams),
      }) 
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
      setLoading(false);
      setSuccessMessage('');
      setErrorMessage('An error occurred while rendering the dashboard');
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