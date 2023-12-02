import Plot from 'react-plotly.js';

function ReactPlotly({plot_data}) {

  const parsedData = JSON.parse(plot_data);

  return (
    <div>
      {parsedData && <Plot data={parsedData.data} layout={parsedData.layout} />}
    </div>
  );
}

export default ReactPlotly;
