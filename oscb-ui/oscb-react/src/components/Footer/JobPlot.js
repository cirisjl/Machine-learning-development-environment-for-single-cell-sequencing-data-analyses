import React from 'react';
import Plot from 'react-plotly.js';

const JobPlot = () => {
    return (
        <div className="w-full flex justify-center items-center my-4">
            <div className="w-full shadow-sm border border-gray-200 bg-white rounded-md">
                <div className="bg-gray-50 px-4 py-2 border-b border-gray-200 rounded-t-md">
                    <h3 className="text-lg font-semibold text-gray-700 m-0">Currently Running and Queued Jobs</h3>
                </div>
                <div className="p-4" style={{ width: '100%', height: '240px', position: 'relative' }}>
                    <Plot
                        data={[]} // Empty data, as requested
                        layout={{
                            autosize: true,
                            margin: { l: 40, r: 20, t: 20, b: 40 },
                            xaxis: { title: 'Time', showgrid: true, zeroline: false },
                            yaxis: { title: 'Job Count', showgrid: true, zeroline: false },
                            hovermode: 'closest',
                            showlegend: true
                        }}
                        style={{ width: '100%', height: '100%' }}
                        useResizeHandler={true}
                        config={{ responsive: true, displayModeBar: false }}
                    />
                </div>
            </div>
        </div>
    );
};

export default JobPlot;
