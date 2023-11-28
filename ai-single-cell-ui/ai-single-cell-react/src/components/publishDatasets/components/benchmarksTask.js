import React, {useEffect, useState} from 'react';
import {CELERY_BACKEND_API} from '../../../constants/declarations';
import {  ScaleLoader } from 'react-spinners';
import AlertMessageComponent from './alertMessageComponent';

function BenchmarksTaskComponent({ setTaskStatus, taskData, setTaskData, setActiveTask, activeTask  }) {

  const [loading, setLoading] = useState(false);
  const [ message, setMessage ] = useState('');
  const [hasMessage, setHasMessage] = useState(message !== '' && message !== undefined);

  const handleTaskCompletion = () => {
    // After Task 6 is successfully completed, update the task status
    setTaskStatus((prevTaskStatus) => ({
      ...prevTaskStatus,
      6: true, // Mark Task 6 as completed
    }));
    //The current task is finished, so make the next task active
    setActiveTask(7);
  };

  useEffect(() => {

      setLoading(true);

      // Extract data needed for the backend API call
      const { task_builder, validation } = taskData;
      const { task_type } = task_builder;
      const { fileMappings } = validation;
      const adataPaths = fileMappings.map((fileDetails) => fileDetails.adata_path);

      // Prepare data for the API call
      const requestData = {
        task_type,
        adata_paths: adataPaths,
        // Add any other data you need to send
      };

      // Perform the API call
      fetch(`${CELERY_BACKEND_API}/convert/publishDatasets/benchmarks`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(requestData),
      })
        .then((response) => {
          if (response.ok) {
            console.log('Backend API call successful');
          } else {
            console.error('Backend API call failed');
            setMessage('Backend API call failed');
          setHasMessage(true);
          }
        })
        .catch((error) => {
          console.error('Error during API call:', error);
          setMessage('Error during API call:', error);
          setHasMessage(true);
        });

        setLoading(false);

  }, [setTaskStatus, taskData]);

  return (
    <div className='benchmarks-task'>
      
      {hasMessage && <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage = {setMessage} />}

      {loading ? (
              <div className="spinner-container">
                <ScaleLoader color="#36d7b7" loading={loading} />
              </div>
            ) : (
      <div className='navigation-buttons'>
            <div className="previous">
              <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                Previous
              </button>
            </div>
            <div className="next-upon-success">
              <button type="submit" className="btn btn-info button" onClick={handleTaskCompletion}>
                Next
              </button>
            </div>
          </div>
        )}
    </div>
  );
}

export default BenchmarksTaskComponent;
