import pandas as pd 
import numpy as np 
import sys

homogenous_samples = sys.argv[1]
kinship_file = sys.argv[2]
output_file = sys.argv[3]

iid = pd.read_csv(homogenous_samples, sep=r'\s+', names=['IID'])['IID'].to_numpy()

kinship = pd.read_csv(kinship_file, sep=r'\s+')
kinship = kinship[
    kinship["ID1"].isin(iid) & kinship["ID2"].isin(iid)
]
kinship = kinship[kinship.Kinship >= 0.044][['ID1','ID2']]

unrel_list = set(iid) - set(kinship.ID1)
unrel_list = unrel_list - set(kinship.ID2)
unrel_list = list(unrel_list)

graph = {}
for i in range(len(kinship)):
    node1, node2 = kinship.iloc[i].values
    if node1 not in graph:
        graph[node1] = []
    if node2 not in graph:
        graph[node2] = []
    graph[node1].append(node2)
    graph[node2].append(node1)

label = {}
for node in graph:
    label[node] = False

for node in graph:
    all_false = True
    for neighbor in graph[node]:
        if label[neighbor]:
            all_false = False
            break
    if all_false:
        label[node] = True
mis = [node for node in graph if label[node]]

unrel_list = unrel_list + mis
unrel_list = np.unique(unrel_list)
print("Length of maximal independent set = " + str(len(unrel_list)))
df = pd.DataFrame(columns=['FID', 'IID'])
df['FID'] = unrel_list
df['IID'] = unrel_list
df.to_csv(output_file, sep='\t', index=None)
