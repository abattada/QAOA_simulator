import random

num_nodes = int(input("Number of nodes: "))
num_edges_wanted = int(input("Number of edges: "))

all_edges = []
for i in range(num_nodes):
    for j in range(i+1, num_nodes):
        all_edges.append((i, j))

# 檢查：如果想要的邊數 > 完全圖邊數，可能造成 sample() 報錯
max_possible_edges = len(all_edges)
if num_edges_wanted > max_possible_edges:
    raise ValueError(f"Cannot sample {num_edges_wanted} edges from only {max_possible_edges} possible edges.")

edges = random.sample(all_edges, num_edges_wanted)

print(f"Random {num_nodes}-node Graph with {num_edges_wanted} edges:")


with open(f"./graph/myGraph_{num_nodes}_{num_edges_wanted}", "w") as f:
    f.write(f"{num_nodes}\n{num_edges_wanted}\n")
    for e in edges:
        i, j = e
        f.write(f"{i} {j}\n")

print(f"./graph/myGraph_{num_nodes}_{num_edges_wanted}")

