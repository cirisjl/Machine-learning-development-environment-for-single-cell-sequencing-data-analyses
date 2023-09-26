import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';

function OptionsList() {
  const [options, setOptions] = useState({});

  useEffect(() => {
    // Define the API URL (update with your actual API endpoint)
    const apiUrl = `${SERVER_URL}/mongoDB/api/addNewOption`;

    // Make an API call to retrieve options
    axios.get(apiUrl)
      .then((response) => {
        setOptions(response.data);
      })
      .catch((error) => {
        console.error('Error fetching data:', error);
      });
  }, []);

    // Function to handle option deletion
    const handleDeleteOption = (optionId) => {
        // Make an API call to delete the option
        // You can use Axios to send a DELETE request to the server
        // Include code to update the UI or handle any errors as needed
      };

  return (
    <div>
      <h2>Options Grouped by Field</h2>
      {Object.keys(options).map((field) => (
        <div key={field}>
          <h3>{field}</h3>
          <ul>
            {options[field].map((option) => (
              <li key={option._id}>
                {option.name}
                {option.username === 'kbcfh' && (
                  <button
                    onClick={() => handleDeleteOption(option._id)}
                    style={{ marginLeft: '10px' }}
                  >
                    Delete
                  </button>
                )}
              </li>
            ))}
          </ul>
        </div>
      ))}
    </div>
  );
}

export default OptionsList;
