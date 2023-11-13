import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';

function UmapPlot ({ umapCoords, obs }) {
  const [umapPlot, setUmapPlot] = useState(null);

  useEffect(() => {
    if (umapCoords && obs) {
      const traces = [];
      obs.forEach((val, i) => {
        const a = obs[val === undefined ? 'undefined' : val];
        const b = umapCoords[val === undefined ? 'undefined' : val];
        const s = Array.from({ length: a.length }, (_, i) => i); // all points are selected by default

        traces.push({
          x: b[0],
          y: b[1],
          text: a.map((id) => `Cell ID: ${id}`),
          mode: 'markers',
          selectedpoints: s,
          marker: {
            size: 7,
            line: { width: 0.5, color: 'grey' },
            color: `rgb(${discrete_colors_3[i % discrete_colors_3.length].substring(1)})`,
          },
          unselected: { marker: { opacity: min_opacity } },
          selected: { marker: { opacity: max_opacity } },
          name: `Cluster ${val}`,
        });
      });

      const umapPlotData = {
        data: traces,
        layout: {
          xaxis: { title: 'UMAP 1' },
          yaxis: { title: 'UMAP 2' },
          margin: { r: 50, l: 50, t: 50, b: 50 },
          legend: { x: 0, y: 1 },
          hovermode: 'closest',
          transition: { duration: 250 },
        },
      };

      setUmapPlot(umapPlotData);
    }
  }, [umapCoords, obs]);

  return (
    <div>
      {umapPlot && <Plot data={umapPlot.data} layout={umapPlot.layout} />}
    </div>
  );
};

export default UmapPlot;
