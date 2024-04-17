import React from 'react';
import { Tooltip } from 'react-tooltip';
import { useState, useEffect } from 'react'
import { SERVER_URL } from '../../constants/declarations';
import axios from 'axios'

const Circle = ({ radiusProp }) => {
  const radius = (radiusProp > 0 && radiusProp < 1 ? radiusProp * 50 : radiusProp) / 5; // Adjust the multiplier as needed for the desired range
  return (
    <>
      <svg data-tooltip-id="my-tooltip" data-tooltip-content={radiusProp} width={50} height={50}>
        <circle cx={25} cy={25} r={radius} fill="#ff8a4f" />
      </svg>
      <Tooltip id="my-tooltip" />
    </>
  );
};

const LeaderCharts = () => {
  const [rowData, setRowData] = useState(null);
  const fetchData = async () => {
    const response = await axios.get(`${SERVER_URL}/mongoDB/api/getLeaderboards`);
    setRowData(response.data.map(item => (
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
      })));
  };

  useEffect(() => {
    fetchData();
  }, []);

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
        {rowData && (<table style={{ maxHeight: "80vh", maxWidth: "100%", margin: "50px 0" }}>
          <thead>
            <tr>
              <td>Dataset Name</td>
              {rowData[0].rowItems.map(rowItem => (
                <td>{rowItem.keyName}</td>
              ))}
            </tr>
          </thead>
          <tbody>
            {rowData.map((row, index) => (
              <tr key={index} style={{ fontSize: '14px', backgroundColor: index % 2 !== 0 ? '#f0f0f0' : '#ffffff' }}>
                <td>{row.datasetId}</td>
                {row.rowItems.map(tdVal => (<td><Circle radiusProp={tdVal.value} /></td>))}
              </tr>
            ))}
          </tbody>
        </table>)}
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
        {rowData && <table>
          <thead>
            <tr>
              <td>Dataset Name</td>
              {rowData[0].rowItems.map(rowItem => (
                <td>{rowItem.keyName}</td>
              ))}
            </tr>
          </thead>
          <tbody>
            {rowData.map((row, index) => (
              <tr key={index} style={{ fontSize: '14px', backgroundColor: index % 2 !== 0 ? '#f0f0f0' : '#ffffff' }}>
                <td>{row.datasetId}</td>
                {row.rowItems.map(tdVal => (<td>{tdVal.value}</td>))}
              </tr>
            ))}
          </tbody>
        </table>}
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