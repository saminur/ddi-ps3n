# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 18:31:09 2021

@author: samin
"""

path = 'G:\\Study\\Research WVU\\PDB\\drug_similarity\\Input Data\\Pre Generated Outputs\\drug_interactions.txt'


import pandas as pd

dictionary=dict()
drugList=[]
with open(path,'r') as token:
    count = 0
    for line in token:
        if count == 0:
            count = count+1
            continue
        split_string = line.split('\t')
        DrugId=split_string[1]
        relation_id=split_string[2]
        # drugList=drugList+[DrugId]
        if dictionary.get(DrugId):
            listVal=dictionary.get(DrugId)
            if type(listVal)==str:
                plusVal=[listVal]+[relation_id]
            else:
                plusVal=listVal+[relation_id]
        
            removeDupl=list(dict.fromkeys(plusVal))
            dictionary[DrugId]=removeDupl
        else:
            dictionary[DrugId]=relation_id
            print(DrugId)
            

import pandas as pd
path_ddi = 'H:\\Research Work\\TrinetXDataSet\\New Results for Jounral\\network_ndd_all_similarity_copy.csv'

df = pd.read_csv(path_ddi)
df = df[df['Predicted Prob'] >= .80]
# Count the rows where Predicted Prob is between 0.99 and 1.00
count_0_99_1_00 = df[(df['Predicted Prob'] >= 0.99) & (df['Predicted Prob'] <= 1.00)].shape[0]

# Count the rows where Predicted Prob is between 0.999 and 1.00
count_0_999_1_00 = df[(df['Predicted Prob'] >= 0.999) & (df['Predicted Prob'] <= 1.00)].shape[0]

# Total number of rows in the dataset
total_rows = df.shape[0]

# Calculate the percentage of rows in the specified ranges
percentage_0_99_1_00 = (count_0_99_1_00 / total_rows) * 100
percentage_0_999_1_00 = (count_0_999_1_00 / total_rows) * 100
count_1_00 = df[df['Predicted Prob'] == 1.0].shape[0]
percentage_1_00 = (count_1_00 / total_rows) * 100 if total_rows > 0 else 0
print(f"Rows with Predicted Prob equal to 1.00: {count_1_00} ({percentage_1_00:.2f}%)")
# Print the results62.
print(f"Rows with Predicted Prob between 0.99 and 1.00: {count_0_99_1_00} ({percentage_0_99_1_00:.2f}%)")
print(f"Rows with Predicted Prob between 0.999 and 1.00: {count_0_999_1_00} ({percentage_0_999_1_00:.2f}%)")

# Van diagram drawing 

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
path_ddi = 'H:\\Research Work\\TrinetXDataSet\\New Results for Jounral\\network_ndd_all_similarity_copy.csv'
data = pd.read_csv(path_ddi)

# Apply filtering criteria
prot_seq_filtered = set(
    data[(data['Predicted Prob'] >= 0.8) & (data['protein similarity'] >= 0.8)]
    .apply(lambda row: (row['Drug ID1'], row['Drug ID2']), axis=1)
)

prot_str_filtered = set(
    data[(data['Predicted Prob'] >= 0.8) & (data['pdb similarity'] >= 0.8)]
    .apply(lambda row: (row['Drug ID1'], row['Drug ID2']), axis=1)
)

predicted_filtered = set(
    data[data['Predicted Prob'] >= 0.8]
    .apply(lambda row: (row['Drug ID1'], row['Drug ID2']), axis=1)
)

# Plot Venn diagram
plt.figure(figsize=(8, 6))
venn3([prot_seq_filtered, prot_str_filtered, predicted_filtered],
      ('Protein Sequence Similarity', 'Protein Structure Similarity', 'Predicted Prob >= 0.8'))
plt.title("Venn Diagram of Filtered Drug-Drug Interactions")
plt.show()



from pyvis.network import Network

# Apply filtering criteria
filtered_data = data[(data['Predicted Prob'] >= 0.8) & 
                     (data['protein similarity'] >= 0.8) & 
                     (data['pdb similarity'] >= 0.8)]



# Create a Pyvis Network
drug_network = Network(height="750px", width="100%", notebook=True, cdn_resources='remote')

drug_network.barnes_hut()

# Add nodes and edges
for _, row in filtered_data.iterrows():
    # Add nodes
    drug_network.add_node(row['Drug ID1'], label=row['Drug ID1'])
    drug_network.add_node(row['Drug ID2'], label=row['Drug ID2'])
    
    # Scale the edge width based on `Predicted Prob` for better visualization
    scaled_width = max(0.1, min(10, row['Predicted Prob'] * 10))  # Adjust scaling as needed
    
    # Add edges with thickness and weight displayed as hover text
    drug_network.add_edge(
        row['Drug ID1'], 
        row['Drug ID2'], 
        width=scaled_width, 
        title=f"Predicted Prob: {row['Predicted Prob']:.2f}"  # Rounded to 2 decimals
    )

# Customize visualization options
drug_network.set_options('''var options = {
  "edges": {
    "color": {
      "inherit": true
    },
    "smooth": {
      "type": "continuous"
    },
    "font": {
      "align": "top",
      "size": 10
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -2000,
      "centralGravity": 0.1,
      "springLength": 100,
      "springConstant": 0.05,
      "avoidOverlap": 0.1
    },
    "minVelocity": 0.75,
    "solver": "barnesHut"
  }
}''')

# Display the network
drug_network.show("drug_network.html")

path = 'H:\\Study\\Research WVU\\PDB\\drug_similarity\\Generated Outputs\\matrices\\aligned\\network_final.csv'
import pandas as pd
import networkx as nx
import plotly.graph_objects as go

# Load CSV file
data = pd.read_csv(path)

# Apply filtering criteria
filtered_data = data

# Create a graph object
graph = nx.Graph()

# Add edges and nodes to the graph
for _, row in filtered_data.iterrows():
    graph.add_edge(row['Drug Id1'], row['Drug Id2'], weight=row['Similarity'])

# Extract positions for network layout
pos = nx.spring_layout(graph, seed=42, k=0.1)  # Adjusted k for tighter clustering

# Create edge traces
edge_trace = []
edge_label_x = []
edge_label_y = []
edge_labels = []

for edge in graph.edges(data=True):
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    weight = edge[2]['weight']
    edge_trace.append(go.Scatter(
        x=[x0, x1, None],
        y=[y0, y1, None],
        line=dict(width=0.5, color='blue'),  # Reduced edge thickness
        mode='lines'))
    # Calculate the midpoint for the edge label
    edge_label_x.append((x0 + x1) / 2)
    edge_label_y.append((y0 + y1) / 2)
    edge_labels.append(f"{weight:.2f}")

# Create edge label trace
edge_label_trace = go.Scatter(
    x=edge_label_x,
    y=edge_label_y,
    text=edge_labels,
    mode='text',
    textfont=dict(size=10, color='red'),
    hoverinfo='none'  # Disable hover info for labels
)

# Create node trace
node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode='markers+text',
    hoverinfo='text',
    marker=dict(
        size=10,
        color='lightblue',
        line=dict(width=2, color='black')
    ),
    textfont=dict(size=8)  # Smaller font for better visualization
)

# Add node positions to the trace
for node in graph.nodes():
    x, y = pos[node]
    node_trace['x'] += tuple([x])
    node_trace['y'] += tuple([y])
    node_trace['text'] += tuple([node])

# Create the figure
fig = go.Figure(data=edge_trace + [edge_label_trace] + [node_trace],
                layout=go.Layout(
                    title='Drug-Drug Interaction Network Diagram',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=0, l=0, r=0, t=40),
                    xaxis=dict(showgrid=False, zeroline=False),
                    yaxis=dict(showgrid=False, zeroline=False)))

# Show the plot
fig.show()


#..................................................................................................
#..............................................................................

import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from networkx.algorithms import community
from plotly.subplots import make_subplots

# Load CSV file
# data = pd.read_csv('your_file.csv')

# Apply filtering criteria
filtered_data = data[(data['Predicted Prob'] >= 0.8) & 
                     (data['protein similarity'] >= 0.8) & 
                     (data['pdb similarity'] >= 0.8)]

# Create a graph object
graph = nx.Graph()

# Add edges and nodes to the graph
for _, row in filtered_data.iterrows():
    graph.add_edge(row['Drug ID1'], row['Drug ID2'], weight=row['Predicted Prob'])

# Identify dense regions using the community detection algorithm
communities = list(community.greedy_modularity_communities(graph))
largest_community = communities[0]  # Select the densest region (largest community)

# Extract positions for network layout
pos = nx.spring_layout(graph, seed=42, k=0.1)  # Adjusted k for tighter clustering

# Full graph visualization
full_edge_trace = []
for edge in graph.edges(data=True):
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    full_edge_trace.append(go.Scatter(
        x=[x0, x1, None],
        y=[y0, y1, None],
        line=dict(width=0.5, color='blue'),
        mode='lines'))

full_node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode='markers+text',
    hoverinfo='text',
    marker=dict(size=10, color='lightblue', line=dict(width=2, color='black')),
    textfont=dict(size=8))

for node in graph.nodes():
    x, y = pos[node]
    full_node_trace['x'] += tuple([x])
    full_node_trace['y'] += tuple([y])
    full_node_trace['text'] += tuple([node])

# Highlight zoomed area in the full graph
bounding_box_x = [pos[node][0] for node in largest_community]
bounding_box_y = [pos[node][1] for node in largest_community]
zoom_box_trace = go.Scatter(
    x=[min(bounding_box_x), max(bounding_box_x), max(bounding_box_x), min(bounding_box_x), min(bounding_box_x)],
    y=[min(bounding_box_y), min(bounding_box_y), max(bounding_box_y), max(bounding_box_y), min(bounding_box_y)],
    mode='lines',
    line=dict(color='green', width=2, dash='dash'),
    showlegend=False)

# Zoomed-in graph (dense region) visualization
subgraph = graph.subgraph(largest_community)
sub_pos = {node: pos[node] for node in subgraph.nodes}

zoom_edge_trace = []
edge_label_x = []
edge_label_y = []
edge_labels = []

for edge in subgraph.edges(data=True):
    x0, y0 = sub_pos[edge[0]]
    x1, y1 = sub_pos[edge[1]]
    weight = edge[2]['weight']
    zoom_edge_trace.append(go.Scatter(
        x=[x0, x1, None],
        y=[y0, y1, None],
        line=dict(width=0.5, color='red'),
        mode='lines'))
    # Midpoint of the edge for labeling
    edge_label_x.append((x0 + x1) / 2)
    edge_label_y.append((y0 + y1) / 2)
    edge_labels.append(f"{weight:.2f}")

# Edge label trace
zoom_edge_label_trace = go.Scatter(
    x=edge_label_x,
    y=edge_label_y,
    text=edge_labels,
    mode='text',
    textfont=dict(size=8, color='darkred'),  # Smaller font for edge weights
    hoverinfo='none')

zoom_node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode='markers+text',
    hoverinfo='text',
    marker=dict(size=12, color='orange', line=dict(width=2, color='black')),
    textfont=dict(size=10))

for node in subgraph.nodes():
    x, y = sub_pos[node]
    zoom_node_trace['x'] += tuple([x])
    zoom_node_trace['y'] += tuple([y])
    zoom_node_trace['text'] += tuple([node])

# Create subplots for full and zoomed-in views
fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=("Full Visualization", "Zoomed-in Dense Region"),
    specs=[[{"type": "scatter"}, {"type": "scatter"}]],
    horizontal_spacing=0.005  
)

# Add traces to subplots
fig.add_traces(full_edge_trace + [full_node_trace, zoom_box_trace], rows=1, cols=1)
fig.add_traces(zoom_edge_trace + [zoom_edge_label_trace, zoom_node_trace], rows=1, cols=2)

# Set layout options
fig.update_layout(
    # title="Graph Visualization with Zoomed-in Dense Region",
    titlefont_size=16,
    showlegend=False,
    hovermode='closest',
    margin=dict(b=0, l=0, r=0, t=40),
    xaxis=dict(showgrid=False, zeroline=False),
    yaxis=dict(showgrid=False, zeroline=False))

fig.update_xaxes(showticklabels=False)
fig.update_yaxes(showticklabels=False)

# Show the plot
fig.show()


#.............................................


