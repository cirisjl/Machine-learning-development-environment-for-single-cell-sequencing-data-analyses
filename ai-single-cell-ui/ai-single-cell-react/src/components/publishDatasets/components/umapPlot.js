import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import PlotConstants from '../plotConstants';

function UmapPlot({ umapCoords, obs }) {
  const [umapPlot, setUmapPlot] = useState(null);

  useEffect(() => {
    if (umapCoords && obs) {
      const traces = [];
      obs.forEach((val, i) => {
        const a = obs[val === undefined ? 'undefined' : val];
        const b = umapCoords[val === undefined ? 'undefined' : val];
        const s = Array.from({ length: a.length }, (_, i) => i);

        traces.push({
          x: b[0],
          y: b[1],
          text: a.map((id) => `Cell ID: ${id}`),
          mode: 'markers',
          selectedpoints: s,
          marker: {
            size: PlotConstants.point_size_2d,
            line: { width: PlotConstants.point_line_width_2d, color: 'grey' },
            color: `rgb(${PlotConstants.discrete_colors_3[i % PlotConstants.discrete_colors_3.length].substring(1)})`,
          },
          unselected: { marker: { opacity: PlotConstants.min_opacity } },
          selected: { marker: { opacity: PlotConstants.max_opacity } },
          name: `Cluster ${val}`,
        });
      });

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
  }, [umapCoords, obs]);

  return (
    <div>
      {umapPlot && <Plot data={umapPlot.data} layout={umapPlot.layout} />}
    </div>
  );
}

export default UmapPlot;
