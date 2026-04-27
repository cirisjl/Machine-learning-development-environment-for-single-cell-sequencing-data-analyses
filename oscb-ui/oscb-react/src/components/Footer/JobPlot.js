import React, { useMemo } from 'react';
import Plot from 'react-plotly.js';

const JobPlot = () => {
    // Generate realistic time-series mock data
    const plotData = useMemo(() => {
        const now = new Date();
        const points = 60;
        const times = [];

        for (let i = points; i >= 0; i--) {
            const t = new Date(now.getTime() - i * 3 * 60000); // 3-min intervals
            times.push(t.toISOString());
        }

        // Generate realistic job count curves
        const generateTrace = (base, variance, spike = false) => {
            return times.map((_, i) => {
                let val = base + Math.sin(i * 0.15) * variance + (Math.random() - 0.5) * variance * 0.5;
                if (spike && i > 30 && i < 40) val += variance * 2;
                return Math.max(0, Math.round(val));
            });
        };

        return {
            times,
            completed: generateTrace(480, 30),
            held: generateTrace(8, 4),
            idle: generateTrace(520, 25),
            suspended: generateTrace(350, 20),
        };
    }, []);

    const traces = [
        {
            x: plotData.times,
            y: plotData.completed,
            name: 'completed',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#34a853', width: 1.5 },
            fill: 'tozeroy',
            fillcolor: 'rgba(52, 168, 83, 0.08)',
        },
        {
            x: plotData.times,
            y: plotData.held,
            name: 'held',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#f9ab00', width: 1.5 },
        },
        {
            x: plotData.times,
            y: plotData.idle,
            name: 'idle',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#4285f4', width: 1.5 },
            fill: 'tozeroy',
            fillcolor: 'rgba(66, 133, 244, 0.05)',
        },
        {
            x: plotData.times,
            y: plotData.suspended,
            name: 'suspended',
            type: 'scatter',
            mode: 'lines',
            line: { color: '#1a237e', width: 1.5 },
            fill: 'tozeroy',
            fillcolor: 'rgba(26, 35, 126, 0.05)',
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
                    letterSpacing: '-0.01em',
                }}>
                    Currently Running and Queued Jobs
                </span>
                <span style={{
                    fontSize: '16px',
                    color: '#999',
                    cursor: 'pointer',
                    lineHeight: 1,
                    padding: '0 4px',
                }}>⋮</span>
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
