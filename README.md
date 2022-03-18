# Proximity NetMed

Proximity is a package for the calculation of network-based distances like proximity and separation.

## Proximity
First introduced in the paper "Network-based in silico drug efficacy screening" by Guney et al. 2016, proximity is a distance measure based on shortest paths between two sets of nodes in the protein-protein interaction network.

## Separation
The separation between two set of nodes is a network-based measure introduced in the paper "Uncovering disease-disease relationships through the incomplete interactome" by Menche et al. 2015.


Implemented by: Rodrigo Dorantes Gilardi [rodgdor@gmail.com](mailto:rodgdor@gmail.com) for the Barabasi Lab.

## Usage
In order to use `proximity.py` you need to clone the repository and add the module to your path as follows:

```bash
git clone https://github.com/Barabasi-Lab/proximity ~/proximity
EXPORT PATH="$HOME/proximity:$PATH"
```

## Example
Here is an example to compute the proximity between a set of genes `T` and a set of disease genes `S`. The script assumes that the protein-protein interaction network is at `'../test_proximity/data/ppi.csv'` and that you have `networkx`, `graph-tool`, `numpy`, and `pandas` installed.

You can select the relevant network science module (`networkx` or `graph-tool`) and remove the code for the other module.

```python3
import proximity
import numpy as np
import pandas as pd
import networkx as nx
import graph_tool.all as gt

df = pd.read_csv('../test_proximity/data/ppi.csv').dropna()
df = df[['Symbol_A', 'Symbol_B']]

G = nx.from_pandas_edgelist(df, source='Symbol_A', target='Symbol_B')

np.random.seed(42)
T = np.random.choice(df['Symbol_A'].values, 10, replace=False)
S = np.random.choice(df['Symbol_B'].values, 10, replace=False)

g = gt.Graph(directed=False)
ids = g.add_edge_list(df.values, hashed=True)
g.vertex_properties['ids'] = ids

net = proximity.Network(g)
alt = proximity.Network(G)

print(net.get_proximity(T, S, n_iter=1000))
print(alt.get_proximity(T, S, n_iter=1000))
```

This will output something like this (result may vary):

```python
{'d_c': 1.8, 'z_score': -0.19576741426849495, 'mu': 1.8270333333333333, 'sigma': 0.13808903506411457}
{'d_c': 1.8, 'z_score': -0.5279439376809244, 'mu': 1.8639444444444444, 'sigma': 0.12111976268792904}
```