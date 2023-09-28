import React from 'react';

function LeftNav({ activeTask, setActiveTask, taskStatus }) {
  const tasks = [
    { id: 1, name: 'Upload File', completed: taskStatus[1] },
    { id: 2, name: 'Validate File', completed: taskStatus[2] },
    { id: 3, name: 'Run APIs', completed: taskStatus[3] },
    // Add other tasks here
  ];

  const handleTaskClick = (task) => {
    const currentIndex = tasks.findIndex((t) => t.id === task.id);

    if (currentIndex === 0 || tasks[currentIndex - 1].completed) {
      setActiveTask(task.id);
    } else {
      alert("Complete the previous task first.");
    }
  };

  return (
    <nav>
      <ul>
        {tasks.map((task) => (
          <li key={task.id} className={task.id === activeTask ? 'active' : ''}>
            <a
              href="#"
              onClick={() => handleTaskClick(task)}
              className={task.completed ? 'completed' : ''}
            >
              {task.name}
            </a>
          </li>
        ))}
      </ul>
    </nav>
  );
}

export default LeftNav;
