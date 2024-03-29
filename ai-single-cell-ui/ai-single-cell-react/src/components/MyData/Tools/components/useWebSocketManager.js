import { useState, useEffect } from 'react';
import { WEB_SOCKET_URL } from '../../../../constants/declarations';

export function useWebSocketManager(taskId, setLiveLogs) {
    const [logs, setLogs] = useState('');
    const [taskStatus, setTaskStatus] = useState('Processing');

    useEffect(() => {
        if (!taskId) return;

        const statusWs = new WebSocket(`${WEB_SOCKET_URL}/taskStatus/${taskId}`);
        const logWs = new WebSocket(`${WEB_SOCKET_URL}/log/${taskId}`);

        const handleStatusMessage = (event) => {
            const data = JSON.parse(event.data);
            if (data[taskId]) {
                setTaskStatus(data[taskId]);
                if (['Success', 'Failed'].includes(data[taskId])) {
                    statusWs.close();
                    logWs.close();
                    console.log("Web socket disconnected");
                }
            }
        };

        const handleLogMessage = (event) => {
            setLiveLogs(event.data);
        };

        statusWs.onmessage = handleStatusMessage;
        logWs.onmessage = handleLogMessage;

        statusWs.onopen = () => console.log('Status WebSocket Connected');
        logWs.onopen = () => console.log('Log WebSocket Connected');

        // Error handling can be more nuanced based on your app's requirements
        const handleError = (error) => console.error('WebSocket Error:', error);
        statusWs.onerror = handleError;
        logWs.onerror = handleError;

        // Cleanup function
        return () => {
            statusWs.close();
            logWs.close();
            console.log("websockets disconnected");
        };
    }, [taskId]); // Reconnect WebSockets if taskId changes

    return { taskStatus};
}
