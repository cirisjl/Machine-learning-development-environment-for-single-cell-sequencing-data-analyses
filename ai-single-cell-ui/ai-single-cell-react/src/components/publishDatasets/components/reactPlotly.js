import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';

function ReactPlotly({plot_data}) {
  const [plotData, setPlotData] = useState(null);

  useEffect(() => {
    if (plot_data) {
      setPlotData(JSON.parse(plot_data));
    }
  }, [plot_data]);

  return (
    <div>
      {plotData && <Plot data={plotData.data} layout={plotData.layout} />}
    </div>
  );
}

export default ReactPlotly;
