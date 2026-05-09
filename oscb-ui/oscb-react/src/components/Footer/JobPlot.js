import React, { useEffect, useMemo, useState } from 'react';
import Plot from 'react-plotly.js';

const MAX_POINTS = 60;
const POLL_INTERVAL_MS = 10000;

const formatLocalDateTime = (date) => date.toLocaleString([], {
    month: 'short',
    day: 'numeric',
    hour: '2-digit',
    minute: '2-digit',
    second: '2-digit',
});

const createSnapshot = (summary = {}) => {
    const receivedAt = new Date();

    return {
        time: receivedAt,
        hoverTime: formatLocalDateTime(receivedAt),
        running: Number(summary.running || 0),
        queued: Number(summary.queued || 0),
    };
};

const getCeleryJobsUrl = () => {
    if (window.location.port === '3000') {
        return `${window.location.protocol}//${window.location.hostname}:5005/api/celery/jobs`;
    }

    return '/api/celery/jobs';
};

const JobPlot = () => {
    const [snapshots, setSnapshots] = useState(() => [createSnapshot()]);
    const [error, setError] = useState(null);
    const [hasLoaded, setHasLoaded] = useState(false);

    useEffect(() => {
        let isMounted = true;

        const fetchJobSnapshot = async () => {
            try {
                const response = await fetch(getCeleryJobsUrl());

                if (!response.ok) {
                    throw new Error(`Request failed with status ${response.status}`);
                }

                const data = await response.json();
                const summary = data.summary || {};
                const snapshot = createSnapshot(summary);

                if (!isMounted) return;

                setSnapshots((prevSnapshots) => [...prevSnapshots, snapshot].slice(-MAX_POINTS));
                setError(null);
                setHasLoaded(true);
            } catch (err) {
                if (!isMounted) return;
                setError(err.message);
            }
        };

        fetchJobSnapshot();
        const intervalId = window.setInterval(fetchJobSnapshot, POLL_INTERVAL_MS);

        return () => {
            isMounted = false;
            window.clearInterval(intervalId);
        };
    }, []);

    const plotData = useMemo(() => {
        return {
            times: snapshots.map((snapshot) => snapshot.time),
            hoverTimes: snapshots.map((snapshot) => snapshot.hoverTime),
            running: snapshots.map((snapshot) => snapshot.running),
            queued: snapshots.map((snapshot) => snapshot.queued),
        };
    }, [snapshots]);
    const latestRunning = plotData.running[plotData.running.length - 1] || 0;
    const latestQueued = plotData.queued[plotData.queued.length - 1] || 0;

    const traces = [
        {
            x: plotData.times,
            y: plotData.running,
            text: plotData.hoverTimes,
            name: 'running',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#34a853', width: 2 },
            fill: 'tozeroy',
            fillcolor: 'rgba(52, 168, 83, 0.08)',
            hovertemplate: '%{text}<br>running: %{y}<extra></extra>',
        },
        {
            x: plotData.times,
            y: plotData.queued,
            text: plotData.hoverTimes,
            name: 'queued',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#f9ab00', width: 2 },
            fill: 'tozeroy',
            fillcolor: 'rgba(249, 171, 0, 0.08)',
            hovertemplate: '%{text}<br>queued: %{y}<extra></extra>',
        },
    ];

    const layout = {
        autosize: true,
        height: 130,
        margin: { l: 36, r: 12, t: 4, b: 28 },
        xaxis: {
            type: 'date',
            showgrid: false,
            tickformat: '%H:%M',
            tickfont: { size: 10, color: '#9e9e9e' },
            linecolor: '#e0e0e0',
            tickcolor: '#e0e0e0',
        },
        yaxis: {
            showgrid: true,
            gridcolor: '#f5f5f5',
            tickfont: { size: 10, color: '#9e9e9e' },
            zeroline: false,
            rangemode: 'nonnegative',
            linecolor: '#e0e0e0',
        },
        showlegend: true,
        legend: {
            orientation: 'h',
            x: 0,
            y: -0.35,
            font: { size: 10, color: '#757575' },
            traceorder: 'normal',
        },
        hovermode: 'x unified',
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
    };

    return (
        <div style={{
            border: '1px solid #e0e0e0',
            borderRadius: '4px',
            backgroundColor: '#fff',
            overflow: 'hidden',
            width: '100%',
        }}>
            {/* Card Header */}
            <div style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                padding: '8px 12px',
                borderBottom: '1px solid #f0f0f0',
                backgroundColor: '#fafafa',
            }}>
                <span style={{
                    fontSize: '13px',
                    fontWeight: 600,
                    color: '#333',
                    letterSpacing: 0,
                }}>
                    Currently Running and Queued Jobs
                </span>
                <span style={{
                    fontSize: '12px',
                    color: error ? '#b91c1c' : '#666',
                    lineHeight: 1,
                    padding: '0 4px',
                }}>
                    {error && !hasLoaded ? 'offline' : `${latestRunning} running / ${latestQueued} queued`}
                </span>
            </div>
            {/* Chart */}
            <div style={{ padding: '4px 8px 0 8px' }}>
                <Plot
                    data={traces}
                    layout={layout}
                    style={{ width: '100%', height: '170px' }}
                    useResizeHandler={true}
                    config={{ responsive: true, displayModeBar: false }}
                />
            </div>
        </div>
    );
};

export default JobPlot;
