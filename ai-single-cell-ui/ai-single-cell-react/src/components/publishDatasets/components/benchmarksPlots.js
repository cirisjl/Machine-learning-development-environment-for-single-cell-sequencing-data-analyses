import React, {useEffect, useState} from 'react';
import Plot from 'react-plotly.js';

function BenchmarksPlots({ barPlot, linePlot }) {

    const [barPlotData, setBarPlotData] = useState(null);
    const [linePlotData, setLinePlotData] = useState(null);

    useEffect(() => {
      if (barPlot) {
        setBarPlotData(barPlot);
      }

      if(linePlot) {
        setLinePlotData(linePlot);
      }
    }, [barPlot, linePlot]);
  return (
    <div>
      <h2>Bar Plot</h2>
      {barPlotData && <Plot data={barPlotData.data} layout={barPlotData.layout} />}

      <h2>Line Plot</h2>
      {linePlotData && <Plot data={linePlotData.data} layout={linePlotData.layout} />}
    </div>
  );
}

export default BenchmarksPlots;