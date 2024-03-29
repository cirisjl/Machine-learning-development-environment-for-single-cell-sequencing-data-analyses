import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { SERVER_URL } from '../../../constants/declarations';

function TaskDetailsComponent() {
  const [tasks, setTasks] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const fetchTasks = async () => {
      try {
        const response = await axios.get(`${SERVER_URL}/api/tools/getTasks`);
        setTasks(response.data);
        setLoading(false);
      } catch (error) {
        console.error('Error fetching tasks:', error);
        setLoading(false);
      }
    };

    fetchTasks();
  }, []);

  return (
    <div>
      <h2>Task Details</h2>
      {loading ? (
        <p>Loading tasks...</p>
      ) : (
        <ul>
          {tasks.map(task => (
            <li key={task._id}>
              <strong>Task Title:</strong> {task.taskTitle}
              <br />
              <strong>Task ID:</strong> {task.taskId}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}

export default TaskDetailsComponent;
