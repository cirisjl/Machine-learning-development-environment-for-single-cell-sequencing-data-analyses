import React, { useState, useEffect } from 'react';
import ReactDOM from 'react-dom';
import { getCookie, isUserAuth} from '../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import DashRenderer from 'dash-renderer';


export default function FlaskDashboard  () {
    const [dashLayout, setDashLayout] = useState(null);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState(null);  const navigate = useNavigate();


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
    .then((response) => {
        if (!response.ok) {
          throw new Error('Failed to fetch the Dash layout.');
        }
        return response.json();
      })
      .then((layoutJson) => {
        setDashLayout(JSON.parse(layoutJson)); // Parse the JSON string into a JavaScript object
        setLoading(false);
      })
      .catch((error) => {
        setError(error.message);
        setLoading(false);
      });

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
      {loading && <div>Loading...</div>}
      {error && <div>Error: {error}</div>}
      {dashLayout && <DashRenderer dashJson={dashLayout} />}
      {/* Your other React components */}
    </div>
  );
};