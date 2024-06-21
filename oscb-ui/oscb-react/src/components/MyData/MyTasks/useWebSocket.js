import { useEffect, useRef } from 'react';
import { WEB_SOCKET_URL } from '../../../constants/declarations';

function useWebSocket(jobId, onStatusMessage, onLogMessage, setLoading) {
  const webSocketStatus = useRef(null);
  const webSocketLog = useRef(null);

  useEffect(() => {
    if (!jobId) {
      console.log("No task ID available for WebSocket connection.");
      return; // Don't proceed if jobId is not set
    }

    // Setup WebSocket for task status updates
    const statusUrl = `${WEB_SOCKET_URL}/taskCurrentStatus/${jobId}`;
    console.log("Connecting to status WebSocket:", statusUrl);
    webSocketStatus.current = new WebSocket(statusUrl);
    webSocketStatus.current.onopen = () => console.log('WebSocket Status Connected:', jobId);
    webSocketStatus.current.onmessage = onStatusMessage;
    webSocketStatus.current.onerror = error => {
      setLoading(false);
      console.error('WebSocket Status Error:', error);
  };
      webSocketStatus.current.onclose = () => console.log('WebSocket for status closed:', jobId);

    // Setup WebSocket for log messages
    const logUrl = `${WEB_SOCKET_URL}/log/${jobId}`;
    console.log("Connecting to log WebSocket:", logUrl);
    webSocketLog.current = new WebSocket(logUrl);
    webSocketLog.current.onopen = () => console.log('WebSocket Log Connected:', jobId);
    webSocketLog.current.onmessage = onLogMessage;
    webSocketLog.current.onerror = error => {
      setLoading(false);
      console.error('WebSocket Log Error:', error);
    }
    webSocketLog.current.onclose = () => console.log('WebSocket for logs closed:', jobId);

    return () => {
      // Cleanup on unmount or jobId change
      console.log('Cleaning up WebSockets for:', jobId);
      if (webSocketStatus.current) {
        webSocketStatus.current.close();
      }
      if (webSocketLog.current) {
        webSocketLog.current.close();
      }
    };
  }, [jobId]);

  const closeWebSockets = () => {
    // Function to manually close WebSockets from the component
    if (webSocketStatus.current) {
      webSocketStatus.current.close();
    }
    // if (webSocketLog.current) {
    //   webSocketLog.current.close();
    // }
  };

  return { closeWebSockets };
}

export default useWebSocket;