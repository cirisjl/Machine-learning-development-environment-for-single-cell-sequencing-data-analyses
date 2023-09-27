import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';

function ManageOptions() {
  const [options, setOptions] = useState({});
  const [selectedOptions, setSelectedOptions] = useState({}); // Track selected options
  const navigate = useNavigate();

  useEffect(() => {
    let jwtToken = getCookie('jwtToken');

    if (!jwtToken) {
      navigate("/routing");
    }

    isUserAuth(jwtToken).then((authData) => {
      if (authData.isAdmin) {
        const username = authData.username;
        const apiUrl = `${SERVER_URL}/mongoDB/api/groupedUserOptions?username=${username}`;

        axios.get(apiUrl)
          .then((response) => {
            setOptions(response.data);
          })
          .catch((error) => {
            console.error('Error fetching data:', error);
          });
      } else {
        console.warn("Unauthorized - you must be an admin to access this page");
        navigate("/accessDenied");
      }
    }).catch((error) => console.error(error));
  }, []);

  // Function to handle option selection
  const handleSelectOption = (field, optionId) => {
    setSelectedOptions((prevSelectedOptions) => ({
      ...prevSelectedOptions,
      [field]: prevSelectedOptions[field]
        ? [...prevSelectedOptions[field], optionId]
        : [optionId],
    }));
    console.log(selectedOptions);
  };

  // Function to handle option deletion
  const handleDeleteSelectedOptions = (field) => {
    // Implement the logic to delete the selected options within the specified field
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
                  <input
                    type="checkbox"
                    checked={selectedOptions[field]?.includes(option._id)}
                    onChange={() => handleSelectOption(field, option._id)}
                    style={{ marginLeft: '10px' }}
                  />
                )}
              </li>
            ))}
          </ul>
          <button onClick={() => handleDeleteSelectedOptions(field)}>
            Delete Selected Options
          </button>
        </div>
      ))}
    </div>
  );
}

export default ManageOptions;
