import React from 'react';
import { ZAxis, Scatter, ScatterChart, Bar, BarChart, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

const Circle = ({ radiusProp }) => {
  const radius = (radiusProp > 0 && radiusProp < 1 ? radiusProp * 50 : radiusProp) / 2; // Adjust the multiplier as needed for the desired range
  return (
    <svg width={100} height={100}>
      <circle cx={50} cy={50} r={radius} fill="#ff8a4f" />
    </svg>
  );
};

const LeaderCharts = () => {
  const headerStyle = {
    fontSize: '12px',
    padding: '5px',
    transform: 'rotate(-45deg)',
    whiteSpace: 'nowrap',
    maxWidth: '80px'
  };

  // const gridStyles = {
  //   parent: {
  //     display: 'grid',
  //     gridTemplateColumns: 'repeat(12, 1fr)',
  //     gridTemplateRows: 'repeat(5, 1fr)',
  //     gridColumnGap: '5px',
  //     gridRowGap: '5px',
  //   },
  //   leftCol: {
  //     gridArea: '2 / 1 / 3 / 2',
  //   },
  //   rightCol: {
  //     gridArea: '2 / 2 / 3 / 13',
  //   },
  //   header: {
  //     gridArea: '1 / 1 / 2 / 13',
  //   },
  // };

  const renderTooltip = (props) => {
    const { active, payload } = props;

    if (active && payload && payload.length) {
      const data = payload[0] && payload[0].payload;

      return (
        <div
          style={{
            backgroundColor: '#fff',
            border: '1px solid #999',
            margin: 0,
            padding: 10,
          }}
        >
          <p>{data.keyName}</p>
          <p>
            <span>value: </span>
            {data.value}
          </p>
        </div>
      );
    }

    return null;
  };

  const data = [
    {
      id: "CL-h-heart-2500-Wang-2024",
      TaskType: {
        label: "Clustering",
        value: "CL"
      },
      TaskLabel: {
        label: "cluster.ids",
        value: "cluster.ids"
      },
      DatasetId: "h-heart-2500-Wang-2024",
      TrainFraction: 0.8,
      ValidationFraction: 0.1,
      TestFraction: 0.1,
      ArchivePath: "/usr/src/app/storage/leijiang/projects/Test-dataset_1712074783936/droplet_Bladder_seurat_tiss_qc_result_data_split.zip",
      BenchmarkResults: {
        metrics: [
          "ARI",
          "Silhouette",
          "NMI"
        ],
        Scanpy: { asw_score: 0.1646, nmi_score: 0.581, ari_score: 0.403, time_points: [1712077886.772374, 1712077887.7800019, 1712077888.782236, 1712077889.7847345, 1712077890.787469, 1712077891.79025, 1712077892.7942264, 1712077893.7964902, 1712077894.7982328, 1712077895.8022242, 1712077896.8045743, 1712077897.8074322, 1712077898.8096998, 1712077899.8120008, 1712077900.8142674, 1712077901.8165302, 1712077902.8287234, 1712077903.832718, 1712077904.8369772, 1712077905.8383331, 1712077906.8408623, 1712077907.8422284, 1712077908.8445282, 1712077909.8492625, 1712077910.8522403, 1712077911.855315, 1712077912.8582358, 1712077913.8620312, 1712077914.8866768, 1712077915.8957734, 1712077916.9014661, 1712077917.9552882, 1712077919.0097656, 1712077920.0185902, 1712077921.0232737, 1712077922.0252178, 1712077923.0269039, 1712077924.0297074, 1712077925.0315673, 1712077926.0385509, 1712077927.0824695, 1712077928.1678581, 1712077929.2606177, 1712077930.3620815, 1712077931.4783, 1712077932.5845902, 1712077933.6845415, 1712077934.7940705, 1712077935.8933537, 1712077937.048438], cpu_usage: [0, 80.2, 98.7, 99.6, 99.6, 98.7, 98.5, 98.3, 98.9, 98.3, 97.8, 99.5, 99.4, 99, 99.3, 99.4, 99.4, 98.3, 99.2, 99.7, 98.6, 97.6, 98.3, 96.8, 97.8, 97.9, 98.1, 85.9, 22.1, 18.5, 1.9, 1.9, 2, 31.7, 32.1, 32.2, 30.5, 31.5, 32.7, 31.2, 1.8, 1.9, 1.9, 1.9, 13.1, 1.9, 1.9, 1.9, 1.8, 1.9], mem_usage: [17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3], gpu_mem_usage: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] },
        bar_plot: { data: [{ type: "bar", x: ["ARI", "Silhouette", "NMI"], y: [0.403, 0.1646, 0.581], text: [0.403, 0.1646, 0.581], textposition: "auto", name: "Scanpy", marker: { opacity: 0.5 } }], layout: { title: { text: "Benchmarks" }, xaxis: { tickangle: 0, type: "category", range: [-0.5, 2.5], autorange: true }, margin: { r: 50, l: 50, t: 50, b: 50 }, barmode: "group", hovermode: "closest", transition: { duration: 100 }, autosize: true, width: 1000, height: 750, yaxis: { type: "linear", range: [0, 0.611578947368421], autorange: true } } },
        line_plot: { data: [{ type: "scatter", x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49], y: [0, 80.2, 98.7, 99.6, 99.6, 98.7, 98.5, 98.3, 98.9, 98.3, 97.8, 99.5, 99.4, 99, 99.3, 99.4, 99.4, 98.3, 99.2, 99.7, 98.6, 97.6, 98.3, 96.8, 97.8, 97.9, 98.1, 85.9, 22.1, 18.5, 1.9, 1.9, 2, 31.7, 32.1, 32.2, 30.5, 31.5, 32.7, 31.2, 1.8, 1.9, 1.9, 1.9, 13.1, 1.9, 1.9, 1.9, 1.8, 1.9], name: "Scanpy_CPU", mode: "lines+markers", line: { dash: "solid" }, marker: { size: 7 }, connectgaps: true }, { type: "scatter", x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49], y: [17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3], name: "Scanpy_Memory", mode: "lines+markers", line: { dash: "dot" }, marker: { size: 7 }, connectgaps: true }], layout: { title: { text: "Computing assessments" }, xaxis: { title: { text: "Time points (s)" }, type: "linear", range: [-3.0676067570506778, 52.067606757050676], autorange: true }, yaxis: { title: { text: "Utilization (%)" }, type: "linear", range: [-6.379934924078091, 106.07993492407809], autorange: true }, margin: { r: 50, l: 50, t: 50, b: 50 }, hovermode: "closest", transition: { duration: 100 }, autosize: true, width: 1000, height: 750 } }
      }
    },
    {
      id: "CL-h-heart-2500-Wang-2024",
      TaskType: {
        label: "Clustering",
        value: "CL"
      },
      TaskLabel: {
        label: "cluster.ids",
        value: "cluster.ids"
      },
      DatasetId: "h-heart-2500-Wang-2024",
      TrainFraction: 0.8,
      ValidationFraction: 0.1,
      TestFraction: 0.1,
      ArchivePath: "/usr/src/app/storage/leijiang/projects/Test-dataset_1712074783936/droplet_Bladder_seurat_tiss_qc_result_data_split.zip",
      BenchmarkResults: {
        metrics: [
          "ARI",
          "Silhouette",
          "NMI"
        ],
        Scanpy: { asw_score: 0.1646, nmi_score: 0.581, ari_score: 0.403, time_points: [1712077886.772374, 1712077887.7800019, 1712077888.782236, 1712077889.7847345, 1712077890.787469, 1712077891.79025, 1712077892.7942264, 1712077893.7964902, 1712077894.7982328, 1712077895.8022242, 1712077896.8045743, 1712077897.8074322, 1712077898.8096998, 1712077899.8120008, 1712077900.8142674, 1712077901.8165302, 1712077902.8287234, 1712077903.832718, 1712077904.8369772, 1712077905.8383331, 1712077906.8408623, 1712077907.8422284, 1712077908.8445282, 1712077909.8492625, 1712077910.8522403, 1712077911.855315, 1712077912.8582358, 1712077913.8620312, 1712077914.8866768, 1712077915.8957734, 1712077916.9014661, 1712077917.9552882, 1712077919.0097656, 1712077920.0185902, 1712077921.0232737, 1712077922.0252178, 1712077923.0269039, 1712077924.0297074, 1712077925.0315673, 1712077926.0385509, 1712077927.0824695, 1712077928.1678581, 1712077929.2606177, 1712077930.3620815, 1712077931.4783, 1712077932.5845902, 1712077933.6845415, 1712077934.7940705, 1712077935.8933537, 1712077937.048438], cpu_usage: [0, 80.2, 98.7, 99.6, 99.6, 98.7, 98.5, 98.3, 98.9, 98.3, 97.8, 99.5, 99.4, 99, 99.3, 99.4, 99.4, 98.3, 99.2, 99.7, 98.6, 97.6, 98.3, 96.8, 97.8, 97.9, 98.1, 85.9, 22.1, 18.5, 1.9, 1.9, 2, 31.7, 32.1, 32.2, 30.5, 31.5, 32.7, 31.2, 1.8, 1.9, 1.9, 1.9, 13.1, 1.9, 1.9, 1.9, 1.8, 1.9], mem_usage: [17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3], gpu_mem_usage: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] },
        bar_plot: { data: [{ type: "bar", x: ["ARI", "Silhouette", "NMI"], y: [0.403, 0.1646, 0.581], text: [0.403, 0.1646, 0.581], textposition: "auto", name: "Scanpy", marker: { opacity: 0.5 } }], layout: { title: { text: "Benchmarks" }, xaxis: { tickangle: 0, type: "category", range: [-0.5, 2.5], autorange: true }, margin: { r: 50, l: 50, t: 50, b: 50 }, barmode: "group", hovermode: "closest", transition: { duration: 100 }, autosize: true, width: 1000, height: 750, yaxis: { type: "linear", range: [0, 0.611578947368421], autorange: true } } },
        line_plot: { data: [{ type: "scatter", x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49], y: [0, 80.2, 98.7, 99.6, 99.6, 98.7, 98.5, 98.3, 98.9, 98.3, 97.8, 99.5, 99.4, 99, 99.3, 99.4, 99.4, 98.3, 99.2, 99.7, 98.6, 97.6, 98.3, 96.8, 97.8, 97.9, 98.1, 85.9, 22.1, 18.5, 1.9, 1.9, 2, 31.7, 32.1, 32.2, 30.5, 31.5, 32.7, 31.2, 1.8, 1.9, 1.9, 1.9, 13.1, 1.9, 1.9, 1.9, 1.8, 1.9], name: "Scanpy_CPU", mode: "lines+markers", line: { dash: "solid" }, marker: { size: 7 }, connectgaps: true }, { type: "scatter", x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49], y: [17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.5, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.1, 17.3, 17.3, 17.3, 17.3, 17.3, 17.3], name: "Scanpy_Memory", mode: "lines+markers", line: { dash: "dot" }, marker: { size: 7 }, connectgaps: true }], layout: { title: { text: "Computing assessments" }, xaxis: { title: { text: "Time points (s)" }, type: "linear", range: [-3.0676067570506778, 52.067606757050676], autorange: true }, yaxis: { title: { text: "Utilization (%)" }, type: "linear", range: [-6.379934924078091, 106.07993492407809], autorange: true }, margin: { r: 50, l: 50, t: 50, b: 50 }, hovermode: "closest", transition: { duration: 100 }, autosize: true, width: 1000, height: 750 } }
      }
    }
  ]

  const rowData = data.map(item => (
    {
      taskType: item.TaskType.label,
      datasetId: item.DatasetId,
      rowItems: [
        // { keyName: 'Train Fraction', index: 1, value: item.TrainFraction },
        // { keyName: 'Validation Fraction', index: 1, value: item.ValidationFraction },
        // { keyName: 'Test Fraction', index: 1, value: item.TestFraction },
        { keyName: 'Scanpy ASW Score', index: 1, value: item.BenchmarkResults.Scanpy.asw_score },
        { keyName: 'Scanpy NMI Score', index: 1, value: item.BenchmarkResults.Scanpy.ari_score },
        { keyName: 'Scanpy ARI Score', index: 1, value: item.BenchmarkResults.Scanpy.nmi_score },
        { keyName: 'Scanpy CPU Usage', index: 1, value: item.BenchmarkResults.Scanpy.cpu_usage.reduce((a, b) => a + b) / item.BenchmarkResults.Scanpy.cpu_usage.length },
        { keyName: 'Scanpy Mem Usage', index: 1, value: item.BenchmarkResults.Scanpy.mem_usage.reduce((a, b) => a + b) / item.BenchmarkResults.Scanpy.mem_usage.length },
        { keyName: 'Scanpy GPU Mem Usage', index: 1, value: item.BenchmarkResults.Scanpy.gpu_mem_usage.reduce((a, b) => a + b) / item.BenchmarkResults.Scanpy.gpu_mem_usage.length },
        ...item.BenchmarkResults.bar_plot.data[0].x.map((xVal, idx) => ({ keyName: xVal, index: 1, value: item.BenchmarkResults.bar_plot.data[0].y[idx] }))
      ]
    }))

  console.log(rowData)

  return (
    <>
      <div>
        <h1> AI Single Cell Dimensionality Reduction Benchmarks Leaderboard</h1>
        <h3>Benchmarking the Efficiency and Accuracy of AI Models for Single Cell Data Analysis</h3>
        <p>
          Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed aliquam, arcu in tincidunt mollis, velit nisl mattis erat, sed iaculis magna justo sed risus. Suspendisse potenti. Sed ac fringilla libero. Proin efficitur eros non ipsum tristique, vitae feugiat eros ullamcorper. Vivamus ullamcorper dui at lectus fermentum, ut condimentum libero vestibulum. Vivamus nec rhoncus sem. Phasellus pretium odio ut lorem placerat, nec viverra sapien pretium. Donec fringilla arcu non dolor dignissim, eget fermentum dui consequat. Integer nec leo vestibulum, aliquam nulla ac, lacinia lectus. Integer iaculis, enim vel vehicula suscipit, ipsum urna cursus urna, eget interdum turpis urna et eros. Sed in orci nec nulla commodo varius. Sed ultrices leo et dui tristique, id scelerisque sapien tincidunt. Nunc posuere metus nec ligula vestibulum, vel placerat nunc varius. Sed feugiat nibh nec erat consequat convallis. Sed ultricies, lacus at finibus finibus, velit elit ultrices ante, nec scelerisque risus nisl eu nunc.
          <br />
          <br />
          Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer convallis dolor nec leo vehicula, vel gravida ipsum lacinia. Nulla facilisi. Nullam non pharetra nunc. Duis sit amet justo a libero tristique tristique vel at justo. In non libero eu nisi varius elementum nec sed ligula. Sed bibendum nec leo vel suscipit. Cras id metus id urna lacinia bibendum nec ut nulla. Curabitur aliquam efficitur arcu, eu finibus libero convallis nec. Nunc non aliquam nisi. Integer sed urna vitae ligula pharetra efficitur ac vel sem. Sed efficitur est nec ante gravida eleifend. Curabitur sed turpis hendrerit, posuere turpis sed, luctus nisi. Sed quis vestibulum sapien. Sed ut purus in est aliquet hendrerit eget eu odio.
        </p>
      </div>
      <div style={{ maxWidth: "100%", overflow: "scroll" }}>
        <table style={{ maxHeight: "80vh", maxWidth: "100%", margin: "50px 0" }}>
          <thead>
            <tr>
              <td>Dataset Name</td>
              {rowData[0].rowItems.map(rowItem => (
                <td>{rowItem.keyName}</td>
              ))}
            </tr>
          </thead>
          <tbody>
            {rowData.map(row => (
              <tr>
                <td>{row.datasetId}</td>
                {row.rowItems.map(tdVal => { console.log(tdVal.value); return (<td><Circle radiusProp={tdVal.value} /></td>) })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div>
        <h2>Metrics</h2>
        <ul>
          <li>
            continuity (Zhang, Shang, and Zhang 2021): Continuity measures error of hard extrusions based on nearest neighbor coranking.
          </li>
          <li>
            Density preservation (Narayan, Berger, and Cho 2021): Similarity between local densities in the high-dimensional data and the reduced data.
          </li>
          <li>
            Density preservation (Narayan, Berger, and Cho 2021): Similarity between local densities in the high-dimensional data and the reduced data.
          </li>
          <li>
            Distance correlation (Schober, Boer, and Schwarte 2018): Spearman correlation between all pairwise Euclidean distances in the original and dimension-reduced data.
          </li>
          <li>
            Distance correlation (spectral) (Coifman and Lafon 2006): Spearman correlation between all pairwise diffusion distances in the original and dimension-reduced data.
          </li>
          <li>
            local continuity meta criterion (Zhang, Shang, and Zhang 2021): The local continuity meta criterion is the co-KNN size with baseline removal which favors locality.
          </li>
          <li>
            global property (Zhang, Shang, and Zhang 2021): The global property metric is a summary of the global co-KNN.
          </li>
        </ul>
      </div>
      <div>
        <h2>Results</h2>
        <table>
          <thead>
            <tr style={{ height: '0px' }}>
              <th className="sorting" style={{ width: '65.9167px', padding: '0px', height: '0px' }} aria-label="Method: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>Method</div>
              </th>
              <th className="sorting" style={{ width: '66.375px', padding: '0px', height: '0px' }} aria-label="Dataset: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>Dataset</div>
              </th>
              <th className="dt-right sorting" style={{ width: '47.5208px', padding: '0px', height: '0px' }} aria-label="Mean score: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>Mean score</div>
              </th>
              <th className="dt-right sorting" style={{ width: '87.5938px', padding: '0px', height: '0px' }} aria-label="continuity: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>continuity</div>
              </th>
              <th className="dt-right sorting" style={{ width: '108.948px', padding: '0px', height: '0px' }} aria-label="Density preservation: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>Density preservation</div>
              </th>
              <th className="dt-right sorting" style={{ width: '94.4792px', padding: '0px', height: '0px' }} aria-label="Distance correlation: activate to sort column ascending">
                <div className="dataTables_sizing" style={{ height: '0px', overflow: 'hidden' }}>Distance correlation</div>
              </th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>densMAP (logCP10k)</td>
              <td>Overall mean</td>
              <td>0.43</td>
              <td>0.65</td>
              <td>0.65</td>
              <td>0.56</td>
            </tr>
            <tr>
              <td>densMAP PCA (logCP10k)</td>
              <td>Overall mean</td>
              <td>0.37</td>
              <td>0.62</td>
              <td>0.47</td>
              <td>0.35</td>
            </tr>
            <tr>
              <td>PyMDE Preserve Distances (logCP10k)</td>
              <td>Overall mean</td>
              <td>0.35</td>
              <td>0.59</td>
              <td>0.29</td>
              <td>0.59</td>
            </tr>
            <tr>
              <td>densMAP PCA (logCP10k, 1kHVG)</td>
              <td>Overall mean</td>
              <td>0.34</td>
              <td>0.60</td>
              <td>0.33</td>
              <td>0.29</td>
            </tr>
            <tr>
              <td>t-SNE (logCP10k)</td>
              <td>Overall mean</td>
              <td>0.34</td>
              <td>0.62</td>
              <td>0.12</td>
              <td>0.36</td>
            </tr>
            <tr>
              <td>densMAP PCA (logCP10k, 1kHVG)</td>
              <td>Overall mean</td>
              <td>0.33</td>
              <td>0.60</td>
              <td>0.38</td>
              <td>0.28</td>
            </tr>
            <tr>
              <td>NeuralEE (CPU) (logCP10k, 1kHVG)</td>
              <td>Overall mean</td>
              <td>0.33</td>
              <td>0.62</td>
              <td>0.15</td>
              <td>0.36</td>
            </tr>
            <tr>
              <td>PyMDE Preserve Neighbors (logCP10k, 1kHVG)</td>
              <td>Overall mean</td>
              <td>0.32</td>
              <td>0.61</td>
              <td>0.01</td>
              <td>0.34</td>
            </tr>
            <tr>
              <td>t-SNE (logCP10k, 1kHVG)</td>
              <td>Overall mean</td>
              <td>0.32</td>
              <td>0.61</td>
              <td>0.07</td>
              <td>0.29</td>
            </tr>
            <tr>
              <td>PyMDE Preserve Neighbors (logCP10k)</td>
              <td>Overall mean</td>
              <td>0.32</td>
              <td>0.62</td>
              <td>0.01</td>
              <td>0.36</td>
            </tr>
          </tbody>
        </table>
      </div>
      <div>
        <h2>References</h2>
        <ul>
          <li>
            10x Genomics. 2019. “5k Peripheral Blood Mononuclear Cells (PBMCs) from a Healthy Donor with a Panel of TotalSeq-b Antibodies (V3 Chemistry).” https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-v-3-chemistry-3-1-standard-3-1-0.
          </li>
          <li>
            Agrawal, Akshay, Alnur Ali, and Stephen Boyd. 2021. “Minimum-Distortion Embedding.” Foundations and Trends in Machine Learning 14 (3): 211–378. https://doi.org/10.1561/2200000090.
          </li>
          <li>
            Coifman, Ronald R., and Stéphane Lafon. 2006. “Diffusion Maps.” Applied and Computational Harmonic Analysis 21 (1): 5–30. https://doi.org/10.1016/j.acha.2006.04.006.
          </li>
          <li>
            McInnes, Leland, John Healy, and James Melville. 2018. “UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.” arXiv. https://doi.org/10.48550/arxiv.1802.03426.
          </li>
          <li>
            Moon, Kevin R., David van Dijk, Zheng Wang, Scott Gigante, Daniel B. Burkhardt, William S. Chen, Kristina Yim, et al. 2019. “Visualizing Structure and Transitions in High-Dimensional Biological Data.” Nature Biotechnology 37 (12): 1482–92. https://doi.org/10.1038/s41587-019-0336-3.
          </li>
          <li>
            Narayan, Ashwin, Bonnie Berger, and Hyunghoon Cho. 2021. “Assessing Single-Cell Transcriptomic Variability Through Density-Preserving Data Visualization.” Nature Biotechnology 39 (6): 765–74. https://doi.org/10.1038/s41587-020-00801-7.
          </li>
          <li>
            Nestorowa, Sonia, Fiona K. Hamey, Blanca Pijuan Sala, Evangelia Diamanti, Mairi Shepherd, Elisa Laurenti, Nicola K. Wilson, David G. Kent, and Berthold Göttgens. 2016. “A Single-Cell Resolution Map of Mouse Hematopoietic Stem and Progenitor Cell Differentiation.” Blood 128 (8): e20–31. https://doi.org/10.1182/blood-2016-05-716480.
          </li>
          <li>
            Nestorowa, Sonia, Fiona K. Hamey, Blanca Pijuan Sala, Evangelia Diamanti, Mairi Shepherd, Elisa Laurenti, Nicola K. Wilson, David G. Kent, and Berthold Göttgens. 2016. “A Single-Cell Resolution Map of Mouse Hematopoietic Stem and Progenitor Cell Differentiation.” Blood 128 (8): e20–31. https://doi.org/10.1182/blood-2016-05-716480.
          </li>
        </ul>
      </div>
    </>
  );
};

export default LeaderCharts;