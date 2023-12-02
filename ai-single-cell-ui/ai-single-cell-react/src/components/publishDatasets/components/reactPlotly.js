import Plot from 'react-plotly.js';

function ReactPlotly({plot_data}) {

  return (
    <div>
      {plot_data && <Plot data={JSON.parse(plot_data).data} layout={JSON.parse(plot_data).layout} />}
    </div>
  );
}

export default ReactPlotly;
