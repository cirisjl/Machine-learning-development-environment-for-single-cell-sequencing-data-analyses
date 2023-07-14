import { useState, useEffect } from 'react';
import { getCookie } from '../../utils/utilFunctions';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faEye, faDownload , faSquarePollVertical } from '@fortawesome/free-solid-svg-icons';
import FilePreviewModal from './filePreviewModal';
import { useNavigate } from 'react-router-dom';

export default function IntermediateFiles({ taskId, results_path, task_title }) {
  const [files, setFiles] = useState([]);
  const [selectedFile, setSelectedFile] = useState(null);
  const [previewBoxOpen, setPreviewBoxOpen] = useState(false);
  const SERVER_URL = "http://" + process.env.REACT_APP_HOST_URL + ":3001";
  const navigate = useNavigate();
  let jwtToken = getCookie('jwtToken');

  useEffect(() => {
    async function fetchFiles() {
      const dirPath = results_path + "/" + taskId
      const response = await fetch(`${SERVER_URL}/getDirContents?dirPath=${dirPath}&authToken=${jwtToken}&usingFor=resultFiles`);
      const data = await response.json();
      setFiles(data.Files);
    }
    fetchFiles();
  }, [taskId]);

  useEffect(() => {
    if (!jwtToken || jwtToken === '')
      navigate('/routing');
  }, [jwtToken]);

  const handlePreviewClick = async (fileName) => {
    setSelectedFile(fileName);
    setPreviewBoxOpen(true);
  };

  function downloadFile(fileUrl) {

    const apiUrl = `${SERVER_URL}/download`;

    const filename = fileUrl.substring(fileUrl.lastIndexOf('/') + 1);
    fetch(`${apiUrl}?fileUrl=${results_path}/${fileUrl}&authToken=${jwtToken}&forResultFile=Yes`)
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
        // document.body.removeChild(link);
      })
      .catch(error => {
        console.error('Error downloading file:', error);
      });
  }

  if (jwtToken)
    return (
      <div>
        <h2>Analysis Results for {task_title}</h2>
        {previewBoxOpen && (
          <FilePreviewModal
            selectedFile={selectedFile}
            setPreviewBoxOpen={setPreviewBoxOpen}
            taskId={taskId}
            jwtToken={jwtToken}
            forResultFile={true}
          />
        )}
        {files.map((file, index) => (
          <div
            key={index}
            style={{
              display: 'flex',
              alignItems: 'center',
              padding: '5px',
              margin: '5px 0', // Add margin between each file option
              cursor: 'pointer',
              backgroundColor: index % 2 === 0 ? '#f9f9f9' : '#e0e0e0', // Alternating background colors
            }}
          >
            <span style={{ marginRight: '10px', flex: '1' }}>{file.name}</span>
            {file.name.endsWith('.h5ad') && (
              <a
                onClick={() => {
                  handlePreviewClick(file.name);
                }}
                target="_blank"
                rel="noopener noreferrer"
                style={{ marginLeft: '10px', textAlign: 'center' }}
              >
                <FontAwesomeIcon icon={faSquarePollVertical} style={{ marginRight: '5px' }} />
                Visualize
              </a>
            )}
            <a
              onClick={() => {
                handlePreviewClick(file.name);
              }}
              target="_blank"
              rel="noopener noreferrer"
              style={{ marginLeft: '10px', textAlign: 'center' }}
            >
              <FontAwesomeIcon icon={faEye} style={{ marginRight: '5px' }} />
              Preview
            </a>
            <a
              download
              onClick={() => {
                downloadFile(file.name);
              }}
              style={{ marginLeft: '10px', textAlign: 'center' }}
            >
              <FontAwesomeIcon icon={faDownload} style={{ marginRight: '5px' }} />
              Download
            </a>
          </div>
        ))}
      </div>
    );
}