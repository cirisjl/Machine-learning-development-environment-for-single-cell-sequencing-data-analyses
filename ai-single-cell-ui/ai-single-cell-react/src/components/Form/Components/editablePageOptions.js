import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';

function ManageOptions() {
  const [options, setOptions] = useState({});
  const [username, setUserName] = useState('');
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
        setUserName(username)
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
    setSelectedOptions((prevSelectedOptions) => {
      const selected = prevSelectedOptions[field] || [];
      if (selected.includes(optionId)) {
        // If the option is already selected, unselect it.
        return {
          ...prevSelectedOptions,
          [field]: selected.filter((id) => id !== optionId),
        };
      } else {
        // If the option is not selected, select it.
        return {
          ...prevSelectedOptions,
          [field]: [...selected, optionId],
        };
      }
    });
    console.log(selectedOptions);
  };

// Function to handle option deletion
const handleDeleteSelectedOptions = (field) => {
  const selectedOptionIds = selectedOptions[field] || [];

  // Check if there are selected options to delete
  if (selectedOptionIds.length > 0) {
    const deleteApiUrl = `${SERVER_URL}/mongoDB/api/deleteOptions`;

    // Send a DELETE request with the array of selected option IDs to delete from MongoDB
    axios
      .delete(deleteApiUrl, { data: { optionIds: selectedOptionIds } })
      .then((response) => {
        // Handle success response, e.g., update the UI to reflect the deleted options
        // After successful deletion, re-fetch the updated options and set the state
        const updatedApiUrl = `${SERVER_URL}/mongoDB/api/groupedUserOptions?username=${username}`;
        axios
          .get(updatedApiUrl)
          .then((response) => {
            setOptions({ ...options, [field]: response.data[field] });
          })
          .catch((error) => {
            console.error('Error fetching updated options:', error);
          });
      })
      .catch((error) => {
        console.error('Error deleting options:', error);
      });

    // Clear the selected options after deletion
    setSelectedOptions((prevSelectedOptions) => ({
      ...prevSelectedOptions,
      [field]: [],
    }));
  }
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
                {option.username === username && (
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
