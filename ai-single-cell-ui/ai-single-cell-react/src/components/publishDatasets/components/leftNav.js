import React from 'react';

function LeftNav({ activeTask, setActiveTask }) {
  const tasks = [
    { id: 1, name: 'Upload File', completed: false },
    { id: 2, name: 'Validate File', completed: false },
    { id: 3, name: 'Run APIs', completed: false },
    // Add other tasks here
  ];

  const handleTaskClick = (task) => {
    if (task.completed || task.id === activeTask) {
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
            <a href="#" onClick={() => handleTaskClick(task)}>
              {task.name}
            </a>
          </li>
        ))}
      </ul>
    </nav>
  );
}

export default LeftNav;
