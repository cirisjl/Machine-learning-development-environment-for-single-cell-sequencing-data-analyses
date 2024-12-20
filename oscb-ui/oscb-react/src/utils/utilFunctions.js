import { NODE_API_URL } from '../constants/declarations'
import axios from 'axios';
import LZString from 'lz-string';

// Function to compress data
export function compressData(data) {
    // Convert data to JSON string
    const jsonString = JSON.stringify(data);
    // Compress and encode in Base64 format
    const compressedString = LZString.compressToBase64(jsonString);
    return compressedString;
}

// Function to decompress data
export function decompressData(compressedString) {
    // Decode and decompress the Base64 string
    const jsonString = LZString.decompressFromBase64(compressedString);
    // Parse JSON string to get original data
    return JSON.parse(jsonString);
}

// Get the value of a cookie with a given name
export function getCookie(name) {
    const cookies = document.cookie.split(';');
    for (let i = 0; i < cookies.length; i++) {
      const cookie = cookies[i].trim();
      if (cookie.startsWith(`${name}=`)) {
        return cookie.substring(`${name}=`.length, cookie.length);
      }
    }
    return '';
  }
  
  // Set a cookie with a given name, value, and expiration time (in days)
  export function setCookie(name, value, expiration) {
    const date = new Date();
    date.setTime(date.getTime() + (expiration * 60 * 1000));
    const expires = "expires=" + date.toUTCString();
    document.cookie = name + "=" + value + ";" + expires + ";path=/";
  }


  export function isUserAuth(jwtToken) {
    return new Promise((resolve, reject) => {
      if (jwtToken) {
        fetch(NODE_API_URL + "/protected", {
          method: 'GET',
          credentials: 'include', // send cookies with the request
          headers: { 'Authorization': `Bearer ${jwtToken}`},
        }) 
        .then((response) => response.json())
        .then((data) => {
          if(data.authData !== null) {
              if (data.authData.username !== null && data.authData.username !== undefined) {
                console.log("Heeloo from isUserAuth ::: " + data.authData.isAdmin);
                resolve({isAuth: true, username: data.authData.username, isAdmin: data.authData.isAdmin});
              } else {
                resolve({isAuth: false, username: null, isAdmin: false});
              }
          }
        })
        .catch((error) => {
          console.error(error);
          // reject(error);
          resolve({isAuth: false, username: null, isAdmin: false});
        });
      }
    });
  }


// Delete a cookie with a given name
export function deleteCookie(name) {
  document.cookie = `${name}=; expires=Thu, 01 Jan 1970 00:00:01 GMT; path=/`;
  return true;
}


export async function getStorageDetails(jwtToken) {
  try {
    const response = await fetch(`${NODE_API_URL}/getStorageDetails?authToken=${jwtToken}`);

    if (response.status === 403) {
      throw new Error('Please log in first');
    }

    const data = await response.json();

    return {
      usedStorage: data.used,
      totalStorage: data.allowed
    };
  } catch (error) {
    if (error.message === 'Please log in first') {
      window.alert('Please log in first');
    } else {
      console.error(error);
    }

    return {
      usedStorage: 0,
      totalStorage: 0
    };
  }
}


export function createUniqueFolderName(title) {
  // Sanitize the title by removing spaces and special characters
  const sanitizedTitle = title
    .replace(/\s+/g, '-') // Replace spaces with hyphens
    .replace(/[^a-zA-Z0-9-]/g, '') // Remove special characters and non-alphanumeric characters

  // Generate a unique identifier (timestamp)
  const timestamp = Date.now();

  // Combine the sanitized title and timestamp to create a unique folder name
  const folderName = `${sanitizedTitle}_${timestamp}`;

  return folderName;
}

export function moveFilesToNewDirectory(newDirectoryPath, isBenchmarks=false) {
  axios
    .post(`${NODE_API_URL}/move-files`, { newDirectoryPath, isBenchmarks, jwtToken:getCookie('jwtToken')})
    .then((response) => {
      console.log('Files moved successfully');
    })
    .catch((error) => {
      // Handle errors if the API call fails.
      console.error('Error moving files', error);
      throw error; // Re-throw the error so it can be caught in the calling code
    });
}


// Data transformation function
export function prepareTableData(cellMetadataObs) {
  const identifiers = new Set();
  Object.values(cellMetadataObs).forEach(values => {
    Object.keys(values).forEach(key => identifiers.add(key));
  });

  // Convert the Set of identifiers into an array and map over it to create rows
  const rows = Array.from(identifiers).map(identifier => {
    const row = { identifier }; // Start each row with the identifier
    Object.entries(cellMetadataObs).forEach(([category, values]) => {
      row[category] = values[identifier] || 'N/A'; // Use 'N/A' for missing values
    });
    return row;
  });

  // Limit the number of rows to 5
  return rows.slice(0, 5);
};

export async function copyFilesToPrivateStorage(selectedFiles, userId){
  try {
    const response = await axios.post(`${NODE_API_URL}/copyFiles`, {
      selectedFiles,
      userId,
    });

    // Check if the server responded with a non-200 status code
    if (response.data.status !== 200) {
      console.error('Server error:', response.data.message);
      return { success: false, message: response.data.message };
    }

    console.log('Success:', response.data.message);
    return { success: true, message: response.data.message };

} catch (error) {
    // Handle error depending on if it's an Axios error or a different error
    let errorMessage = error.response ? error.response.data.message : error.message;
    return { success: false, message: errorMessage };  }
};


export function downloadFile(fileUrl) {
  const apiUrl = `${NODE_API_URL}/download`;
  const pwd = "jobResults";

  if (fileUrl) {
    const filename = fileUrl.substring(fileUrl.lastIndexOf('/') + 1);

    fetch(`${apiUrl}?fileUrl=${fileUrl}&authToken=${getCookie("jwtToken")}&pwd=${pwd}`)
      .then(response => {
        return response.blob();
      })
      .then(blob => {
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;

        document.body.appendChild(link);
        link.click();
        // Remove the link from the DOM
        // document.body.removeChild(link);
      })
      .catch(error => {
        console.error('Error downloading file:', error);
      });
  }
};

export function getFileNameFromURL(fileUrl){
  if (fileUrl) {
    return fileUrl.substring(fileUrl.lastIndexOf('/') + 1);
  } else{
    return '';
  }
};