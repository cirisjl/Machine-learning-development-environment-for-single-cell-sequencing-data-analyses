import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import PlotConstants from '../plotConstants';

function UmapPlot({ traces}) {
  const [umapPlot, setUmapPlot] = useState(null);

  useEffect(() => {
    if (traces) {
      const umapPlotData = {
        data: traces,
        layout: {
          xaxis: { title: 'UMAP 1' },
          yaxis: { title: 'UMAP 2' },
          margin: PlotConstants.margin,
          legend: { x: 0, y: 1 },
          hovermode: 'closest',
          transition: { duration: 250 },
        },
      };
      setUmapPlot(umapPlotData);
    }
  }, [traces]);

  return (
    <div>
      {umapPlot && <Plot data={umapPlot.data} layout={umapPlot.layout} />}
    </div>
  );
}

export default UmapPlot;
