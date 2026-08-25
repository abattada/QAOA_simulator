import pennylane as qml
import re
from pennylane import numpy as np
# print(qml.data.list_attributes(data_name="hamlib-maxcut"))
[dataset] = qml.data.load(data_name="hamlib-maxcut", attributes=["edges", "ids"])
# print(dataset.edges)
for i in range(102, 153):
    graph_data = dataset.ids[i]
    data_name = re.findall(r"b-(\d*)", graph_data)[0]
    file_name = f"./graph/qml_bi_{data_name}"
    nodes = int(data_name) *2
    edge_data = dataset.edges[i]
    edges = len(edge_data)
    with open(file_name, 'w') as file:
        file.write(str(nodes)+'\n')
        file.write(str(edges)+'\n')
        for (f,s) in edge_data:
            file.write(str(int(f)) + " " + str(int(s)) + '\n')
