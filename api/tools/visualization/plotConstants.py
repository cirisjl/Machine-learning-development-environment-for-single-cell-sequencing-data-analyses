# Scale of plot sizes
# All plot geometry is expressed as multiples
# of these parameters
scale = 250
scaleratio = 1.0
pt_expression_scaleratio = 0.5
violin_expression_scaleratio = 1.5

# margins on plots
margin = {"r": 50, "l": 50, "t": 50, "b": 50}

# point sizes
point_line_width_2d = 0.5
point_line_width_3d = 0.5
point_size_2d = 7
point_size_3d = 2.5
point_size_pt_trend = 2

# min and max opacity of points in scatter plots
min_opacity = 0.15
max_opacity = 1

import plotly.express as px

discrete_colors_0 = ["#e28e31",
"#8a9bde",
"#9f5036",
"#8a5ad2",
"#a2b937",
"#59c8b5",
"#e07d93",
"#406caa",
"#4ab4dd",
"#9d4564",
"#38977f",
"#65c14c",
"#d288c3",
"#d175df",
"#c1303c",
"#bdb466",
"#7a81e0",
"#dd3d72",
"#93a95f",
"#6b8627",
"#e26d69",
"#3f9335",
"#8e6e2e",
"#caa637",
"#72b879",
"#377945",
"#636c29",
"#a439a6",
"#db976c",
"#82559d",
"#4a61d1",
"#d0499c",
"#c6662e",
"#39c685",
"#dd4f2e"]

discrete_colors_1 = ["#d1ff89",
"#808fbb",
"#ce729b",
"#00d9cf",
"#9cd2ff",
"#b077db",
"#ffbb5c",
"#02dd95",
"#ffb6d3",
"#709c4e",
"#bcffeb",
"#ff97a1",
"#65acff",
"#ff8a77",
"#b5b600",
"#b77bb5",
"#6bffd8",
"#a5ff9d",
"#00af4e",
"#ff8c5d",
"#ffdcaa",
"#03cbe3",
"#af8a3b",
"#bd8900",
"#cdffc6",
"#d0bbff",
"#00be93",
"#d364d0",
"#f893ff",
"#99ff64",
"#6e9e1b",
"#85bca7",
"#31b600",
"#ff94e1",
"#1a99d3"]

discrete_colors_2 = px.colors.qualitative.Alphabet

discrete_colors_3 = px.colors.qualitative.D3 + px.colors.qualitative.Set3 + px.colors.qualitative.T10 + px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet