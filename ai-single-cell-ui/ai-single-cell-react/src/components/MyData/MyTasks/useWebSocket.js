// import { useEffect, useRef } from 'react';
// import {WEB_SOCKET_URL} from '../../../constants/declarations';

// function useWebSocket(taskId, onStatusMessage, onLogMessage) {
//   const webSocketStatus = useRef(null);
//   const webSocketLog = useRef(null);

//   useEffect(() => {
//     if (!taskId) {
//       console.log("No task ID available for WebSocket connection.");
//       return; // Don't proceed if taskId is not set
//     }

//     // Setup WebSocket for task status updates
//     const statusUrl = `${WEB_SOCKET_URL}/taskStatus/${taskId}`;
//     console.log("Connecting to status WebSocket:", statusUrl);
//     webSocketStatus.current = new WebSocket(statusUrl);
//     webSocketStatus.current.onopen = () => console.log('WebSocket Status Connected:', taskId);
//     webSocketStatus.current.onmessage = onStatusMessage;
//     webSocketStatus.current.onerror = error => console.error('WebSocket Status Error:', error);
//     webSocketStatus.current.onclose = () => console.log('WebSocket for status closed:', taskId);

//     // Setup WebSocket for log messages
//     const logUrl = `${WEB_SOCKET_URL}/log/${taskId}`;
//     console.log("Connecting to log WebSocket:", logUrl);
//     webSocketLog.current = new WebSocket(logUrl);
//     webSocketLog.current.onopen = () => console.log('WebSocket Log Connected:', taskId);
//     webSocketLog.current.onmessage = onLogMessage;
//     webSocketLog.current.onerror = error => console.error('WebSocket Log Error:', error);
//     webSocketLog.current.onclose = () => console.log('WebSocket for logs closed:', taskId);

//     return () => {
//       // Cleanup on unmount or taskId change
//       console.log('Cleaning up WebSockets for:', taskId);
//       if (webSocketStatus.current) {
//         webSocketStatus.current.close();
//       }
//       if (webSocketLog.current) {
//         webSocketLog.current.close();
//       }
//     };
//   }, [taskId]); // Ensure WebSocket is only re-initialized if taskId changes
// }


// export default useWebSocket;


import { useEffect, useRef } from 'react';
import { WEB_SOCKET_URL } from '../../../constants/declarations';

function useWebSocket(taskId, onStatusMessage, onLogMessage) {
  const webSocketStatus = useRef(null);
  const webSocketLog = useRef(null);

  useEffect(() => {
    if (!taskId) {
      console.log("No task ID available for WebSocket connection.");
      return; // Don't proceed if taskId is not set
    }

    // Setup WebSocket for task status updates
    const statusUrl = `${WEB_SOCKET_URL}/taskCurrentStatus/${taskId}`;
    console.log("Connecting to status WebSocket:", statusUrl);
    webSocketStatus.current = new WebSocket(statusUrl);
    webSocketStatus.current.onopen = () => console.log('WebSocket Status Connected:', taskId);
    webSocketStatus.current.onmessage = onStatusMessage;
    webSocketStatus.current.onerror = error => console.error('WebSocket Status Error:', error);
    webSocketStatus.current.onclose = () => console.log('WebSocket for status closed:', taskId);

    // Setup WebSocket for log messages
    const logUrl = `${WEB_SOCKET_URL}/log/${taskId}`;
    console.log("Connecting to log WebSocket:", logUrl);
    webSocketLog.current = new WebSocket(logUrl);
    webSocketLog.current.onopen = () => console.log('WebSocket Log Connected:', taskId);
    webSocketLog.current.onmessage = onLogMessage;
    webSocketLog.current.onerror = error => console.error('WebSocket Log Error:', error);
    webSocketLog.current.onclose = () => console.log('WebSocket for logs closed:', taskId);

    return () => {
      // Cleanup on unmount or taskId change
      console.log('Cleaning up WebSockets for:', taskId);
      if (webSocketStatus.current) {
        webSocketStatus.current.close();
      }
      if (webSocketLog.current) {
        webSocketLog.current.close();
      }
    };
  }, [taskId]);

  const closeWebSockets = () => {
    // Function to manually close WebSockets from the component
    if (webSocketStatus.current) {
      webSocketStatus.current.close();
    }
    if (webSocketLog.current) {
      webSocketLog.current.close();
    }
  };

  return { closeWebSockets };
}

export default useWebSocket;
