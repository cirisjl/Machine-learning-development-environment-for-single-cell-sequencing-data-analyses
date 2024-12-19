import { NODE_API_URL } from '../constants/declarations'
import axios from 'axios';
import LZString from 'lz-string';
import pako from 'pako'; // Import pako, a zlib-compatible library for browsers

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



export function plotUmapObs(cellMetadata, umap, clusteringPlotType, selectedCellIntersection = [], annotation = null, nDim = 2) {
  // Validate if the clustering ID exists
  cellMetadata = gunzipDict(cellMetadata);
  umap = gunzipDict(umap);
  console.log("Cell Metadata : " + cellMetadata);
  console.log(umap);
  if (!cellMetadata[clusteringPlotType]) {
    const validClusterIds = ['cluster.ids', 'leiden', 'louvain', 'seurat_clusters'];
    clusteringPlotType = validClusterIds.find((id) => cellMetadata[id]) || null;

    if (!clusteringPlotType) {
      throw new Error(`Clustering type ${clusteringPlotType} does not exist in cell metadata.`);
    }
  }

  const coords = umap.map((row) => ({ x: row[0], y: row[1], z: row[2] || null })); // Convert UMAP array to an array of coordinates
  const clusters = [...new Set(cellMetadata[clusteringPlotType])]; // Extract unique clusters

  const traces = clusters.map((cluster, index) => {
    // Filter metadata and coordinates for the current cluster
    const clusterIndices = cellMetadata[clusteringPlotType]
      .map((val, idx) => (val === cluster ? idx : null))
      .filter((idx) => idx !== null);

    const filteredCoords = clusterIndices.map((idx) => coords[idx]);
    const selectedPoints = selectedCellIntersection.length
      ? selectedCellIntersection
          .map((cell) => cellMetadata.indexOf(cell))
          .filter((idx) => clusterIndices.includes(idx))
      : clusterIndices;

    const textAnnotations = annotation
      ? clusterIndices.map((idx) => String(cellMetadata[annotation][idx]))
      : clusterIndices.map((idx) => `Cell ID: ${cellMetadata.index[idx]}`);

    // Generate trace for 2D or 3D plot
    if (nDim === 2) {
      return {
        type: 'scattergl',
        x: filteredCoords.map((point) => point.x),
        y: filteredCoords.map((point) => point.y),
        text: textAnnotations,
        selectedpoints: selectedPoints,
        mode: 'markers',
        marker: {
          size: 5, // Customize size
          line: { width: 1, color: 'grey' },
          color: discrete_colors_3[index % discrete_colors_3.length],
        },
        unselected: { marker: { opacity: 0.3 } },
        selected: { marker: { opacity: 1 } },
        name: `Cluster ${cluster}`,
      };
    } else if (nDim === 3) {
      return {
        type: 'scatter3d',
        x: filteredCoords.map((point) => point.x),
        y: filteredCoords.map((point) => point.y),
        z: filteredCoords.map((point) => point.z),
        text: textAnnotations,
        selectedpoints: selectedPoints,
        mode: 'markers',
        marker: {
          size: 3, // Customize size
          line: { width: 1, color: 'grey' },
          color: discrete_colors_3[index % discrete_colors_3.length],
        },
        name: `Cluster ${cluster}`,
      };
    }

    throw new Error(`Unsupported dimension: ${nDim}`);
  });

  // Return data and layout for the plot
  return {
    data: traces,
    layout: {
      xaxis: { title: 'UMAP 1' },
      yaxis: { title: 'UMAP 2' },
      ...(nDim === 3 && { zaxis: { title: 'UMAP 3' } }),
      margin: { l: 40, r: 40, t: 40, b: 40 },
      hovermode: 'closest',
      autosize: true,
      width: 800,
      height: 600,
    },
  };
};

const gunzipDict = (toUngzip) => {
  try {
    // Decode base64 string to a Uint8Array (binary data)
    const compressedBuffer = Uint8Array.from(atob(toUngzip), c => c.charCodeAt(0));

    // Decompress the data using pako.ungzip (equivalent to zlib.gunzip in Node.js)
    const decompressedBuffer = pako.ungzip(compressedBuffer, { to: 'string' });

    // Parse the decompressed JSON string into a JavaScript object
    return JSON.parse(decompressedBuffer);
  } catch (error) {
    console.error('Error during decompression:', error);
    throw error;
  }
};

// Scale of plot sizes
// All plot geometry is expressed as multiples
// of these parameters
const scale = 250;
const scaleratio = 1.0;
const pt_expression_scaleratio = 0.5;
const violin_expression_scaleratio = 1.5;

// Margins on plots
const margin = { r: 50, l: 50, t: 50, b: 50 };

// Point sizes
const point_line_width_2d = 0.5;
const point_line_width_3d = 0.5;
const point_size_2d = 7;
const point_size_3d = 2.5;
const point_size_pt_trend = 2;

// Min and max opacity of points in scatter plots
const min_opacity = 0.15;
const max_opacity = 1;

// Discrete colors
const discrete_colors_0 = [
  "#e28e31", "#8a9bde", "#9f5036", "#8a5ad2", "#a2b937", "#59c8b5",
  "#e07d93", "#406caa", "#4ab4dd", "#9d4564", "#38977f", "#65c14c",
  "#d288c3", "#d175df", "#c1303c", "#bdb466", "#7a81e0", "#dd3d72",
  "#93a95f", "#6b8627", "#e26d69", "#3f9335", "#8e6e2e", "#caa637",
  "#72b879", "#377945", "#636c29", "#a439a6", "#db976c", "#82559d",
  "#4a61d1", "#d0499c", "#c6662e", "#39c685", "#dd4f2e"
];

const discrete_colors_1 = [
  "#d1ff89", "#808fbb", "#ce729b", "#00d9cf", "#9cd2ff", "#b077db",
  "#ffbb5c", "#02dd95", "#ffb6d3", "#709c4e", "#bcffeb", "#ff97a1",
  "#65acff", "#ff8a77", "#b5b600", "#b77bb5", "#6bffd8", "#a5ff9d",
  "#00af4e", "#ff8c5d", "#ffdcaa", "#03cbe3", "#af8a3b", "#bd8900",
  "#cdffc6", "#d0bbff", "#00be93", "#d364d0", "#f893ff", "#99ff64",
  "#6e9e1b", "#85bca7", "#31b600", "#ff94e1", "#1a99d3"
];

// Colors from Plotly's qualitative color scales
const discrete_colors_2 = [
  "#F0F0F0", "#D4D4D4", "#B8B8B8", "#9C9C9C", "#808080", "#646464", 
  "#484848", "#2C2C2C", "#101010", "#F9F9F9", "#B1B1B1", "#9E9E9E",
  "#D9D9D9", "#A5A5A5", "#8C8C8C", "#707070", "#565656", "#3A3A3A", 
  "#1E1E1E", "#8F8F8F", "#D1D1D1", "#B2B2B2", "#9C9C9C", "#919191"
];

// Colors from multiple Plotly qualitative color scales combined (D3, Set3, T10, Plotly, Alphabet)
const discrete_colors_3 = [
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#8dd3c7", "#ffffb3", 
  "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", 
  "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#ffff99", "#c4e1ff", 
  "#eb8d2b", "#d7d7d7", "#ff94b5", "#a6a6a6", "#c5b0b0", "#8d9b9e"
];