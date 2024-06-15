import { useEffect, useRef } from 'react';

function useWebSockets(taskId, handleStatusMessage, WEB_SOCKET_URL) {
    const webSocketsRef = useRef({});

    useEffect(() => {
        if (!taskId) return;

        const statusUrl = `${WEB_SOCKET_URL}/taskCurrentStatus/${taskId}`;
        if (!webSocketsRef.current[taskId]) {
            console.log("Opening new WebSocket connection for:", taskId);
            const ws = new WebSocket(statusUrl);
            ws.onopen = () => console.log('WebSocket Connected:', taskId);
            ws.onmessage = handleStatusMessage;
            ws.onerror = error => console.error('WebSocket Error:', taskId, error);
            ws.onclose = () => console.log('WebSocket Closed:', taskId);
            webSocketsRef.current[taskId] = ws;
        }

        return () => {
            // Cleanup function to close the WebSocket when component unmounts or taskId changes
            const ws = webSocketsRef.current[taskId];
            if (ws) {
                console.log('Cleaning up WebSocket:', taskId);
                ws.close();
                delete webSocketsRef.current[taskId];
            }
        };
    }, [taskId]);

    const closeWebSocket = (id) => {
        const ws = webSocketsRef.current[id];
        if (ws) {
            console.log('Manually closing WebSocket:', id);
            ws.close();
            delete webSocketsRef.current[id];
        }
    };

    return { closeWebSocket };
}


export default useWebSockets;