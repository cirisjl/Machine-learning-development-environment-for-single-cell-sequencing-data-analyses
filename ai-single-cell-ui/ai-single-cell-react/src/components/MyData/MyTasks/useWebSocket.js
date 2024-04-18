import { useEffect, useState, useRef } from 'react';

function useWebSocket(url, onMessage) {
  const websocket = useRef(null);

  useEffect(() => {
    // Initialize WebSocket connection
    websocket.current = new WebSocket(url);
    websocket.current.onopen = () => {
      console.log('WebSocket Connected');
    };
    websocket.current.onerror = (error) => {
      console.error('WebSocket Error', error);
    };
    websocket.current.onmessage = onMessage;

    return () => {
      // Clean up the WebSocket connection
      websocket.current.close();
      console.log('WebSocket Disconnected');
    };
  }, [url]);
}


export default useWebSocket;