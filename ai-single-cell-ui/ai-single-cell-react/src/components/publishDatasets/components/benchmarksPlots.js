import React from 'react';
import ReactPlotly from './reactPlotly';

function BenchmarksPlots({ barPlot, linePlot }) {
  return (
    <div>
      <h2>Bar Plot</h2>
      <ReactPlotly plotData={barPlot} />

      <h2>Line Plot</h2>
      <ReactPlotly plotData={linePlot} />
    </div>
  );
}

export default BenchmarksPlots;