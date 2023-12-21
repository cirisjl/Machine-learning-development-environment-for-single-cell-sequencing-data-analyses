import React from 'react';

function LeftNav({ activeTask, setActiveTask, taskStatus, taskData, setTaskData, flow }) {

  const uploadTasks = [
    { id: 1, name: 'Upload', completed: taskStatus[1] },
    { id: 2, name: 'Validation', completed: taskStatus[2] },
    { id: 3, name: 'QC', completed: taskStatus[3] },
    { id: 4, name: 'Metadata', completed: taskStatus[4] },
    // Add other upload tasks here
  ];

  const tbTasks = [
    { id: 5, name: 'Task Builder', completed: taskStatus[5] },
    { id: 6, name: 'Benchmarks', completed: taskStatus[6] },
    { id: 7, name: 'Review', completed: taskStatus[7] },
    // Add other task builder tasks here
  ];

  const tasks = flow === 'taskBuilder' ? tbTasks : uploadTasks;

  const handleTaskClick = (task) => {
    const currentIndex = tasks.findIndex((t) => t.id === task.id);

    if (currentIndex === 0 || tasks[currentIndex - 1].completed) {
      setActiveTask(task.id);
    } else {
      alert("Complete the previous task first.");
    }
  };

  return (
    <nav className='benchmarks-left-nav'>
      <ul>
        {tasks.map((task) => (
          <li key={task.id}
          className={`${
            task.id === activeTask ? 'active' : ''
          } ${task.completed ? 'completed' : 'disable-click'}`}
          >
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
