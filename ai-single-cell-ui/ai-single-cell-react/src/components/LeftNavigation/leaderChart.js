import React from 'react';
import { Bar, BarChart, LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

const LeaderCharts = () => {

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
    }
  ]

  const barData = data[0].BenchmarkResults.bar_plot.data[0].x.map((item, index) => ({
    x: item,
    y: data[0].BenchmarkResults.bar_plot.data[0].y[index]
  }))

  const scanpyData = data[0].BenchmarkResults.Scanpy.time_points.map((time, index) => ({
    time_point: time,
    cpu_usage: data[0].BenchmarkResults.Scanpy.cpu_usage[index],
    gpu_mem_usage: data[0].BenchmarkResults.Scanpy.gpu_mem_usage[index]
  }));

  console.log(barData)
  return (
    <div>
      <ResponsiveContainer width="100%" height={400}>
        <LineChart data={scanpyData} margin={{ top: 20, right: 30, left: 20, bottom: 10 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="time_point" type="number" domain={['dataMin', 'dataMax']} tickCount={10} tickFormatter={(time) => new Date(time * 1000).toLocaleTimeString()} />
          <YAxis />
          <Tooltip labelFormatter={(time) => new Date(time * 1000).toLocaleTimeString()} />
          <Legend />
          <Line type="monotone" dataKey="cpu_usage" name="CPU Usage" stroke="#8884d8" />
          <Line type="monotone" dataKey="gpu_mem_usage" name="GPU Memory Usage" stroke="#82ca9d" />
        </LineChart>
      </ResponsiveContainer>

      <ResponsiveContainer width="100%" height={400}>
        <BarChart data={barData} margin={{ top: 20, right: 30, left: 20, bottom: 10 }}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="x" />
          <YAxis />
          <Tooltip />
          <Bar dataKey="y" fill="#8884d8" />
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
};

export default LeaderCharts;