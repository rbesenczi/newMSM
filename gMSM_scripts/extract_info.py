import csv
from os import environ
import graph_tool.all as gt

home = environ['HOME']
dataset = "HCP"
workdir = home + "/groupwise/" + dataset
clustering = workdir + "/frontal_subject_clusters_hcp.csv"
hierarchy = workdir + "/frontal_hierarchical_path_all.csv"
root = "NODE2218"

groups = {}
no_small_groups = []
paths = []
cg_path = {}

def proc_path(root):

	leaves = []
	for path in paths:
		try:
			idx = path.index(root)
		except:
			continue
		leaves.append(path[path.index(root)-1])
	leaves = list(set(leaves))
	cg_path[root] = leaves
	if leaves[0] not in no_small_groups:
		proc_path(leaves[0])
	if len(leaves) == 2 and leaves[1] not in no_small_groups:
		proc_path(leaves[1])


class Visitor(gt.DFSVisitor):

    def __init__(self, name, pred, time, paths, root):
        self.name = name
        self.pred = pred
        self.time = time
        self.last_time = 0
        self.path = []
        self.root = root

    def discover_vertex(self, u):
        self.path.append(self.name[u])
        self.time[u] = self.last_time
        self.last_time += 1
        if self.name[u] == self.root:
        	paths.append(self.path)

    def tree_edge(self, e):
        self.pred[e.target()] = int(e.source())

###############################################################################################

with open(clustering, "r", newline='') as csvfile:
	reader = csv.DictReader(csvfile, fieldnames=['line','subject','group'])
	for row in reader:
		if row['group'] in groups:
			groups[row['group']].append(row['subject'].split('\n')[0])
		else:
			groups[row['group']] = [ row['subject'].split('\n')[0] ]

	file = open(workdir + "/group_list.txt", 'w')

	for key in groups.keys():
		num_subs = len(groups[key])
		if num_subs > 4:
			no_small_groups.append(key)
			file.write(str(key) + ',' + str(num_subs) + '\n')
	
	file.close()

###############################################################################################

g = gt.Graph(directed=True)

with open(hierarchy, "r", newline='') as hierarchy_file:
	reader = csv.DictReader(hierarchy_file, fieldnames=['left','right','root'])
	src_trg_pair = []
	for row in reader:
		src_trg_pair.append((row['left'], row['root'], row['left']+row['root']))
		src_trg_pair.append((row['right'], row['root'], row['right']+row['root']))
	vertex_ids = g.add_edge_list(src_trg_pair, hashed=True, hash_type="string", eprops=[('name', 'string')])

key_to_id = { vertex_ids[i]: i for i in range(g.num_vertices()) }

name = vertex_ids
time = g.new_vertex_property("int")
pred = g.new_vertex_property("int64_t")

for key in no_small_groups:
	gt.dfs_search(g, g.vertex(key_to_id[key]), Visitor(name, pred, time, paths, root))

proc_path(root)

cg_path = dict(sorted(cg_path.items()))

###############################################################################################

lone_leaves = {}
iteration = True

while iteration == True:
	iteration = False
	for key in cg_path.keys():
		if len(cg_path[key]) == 1:
			iteration = True
			if cg_path[key][0] in lone_leaves.keys():
				cg_path[key][0] = lone_leaves[cg_path[key][0]]
			lone_leaves[key] = cg_path[key][0]

	for key in cg_path.keys():
		if cg_path[key][0] in lone_leaves.keys():
			cg_path[key][0] = lone_leaves[cg_path[key][0]]
		if len(cg_path[key]) == 2:
			if cg_path[key][1] in lone_leaves.keys():
				cg_path[key][1] = lone_leaves[cg_path[key][1]]

	for key in lone_leaves.keys():
		cg_path.pop(key)
	lone_leaves.clear()

file = open(workdir + "/frontal_hierarchical_path_study.csv", "w")
for key in cg_path.keys():
	file.write(cg_path[key][0] + ',' + cg_path[key][1] + ',' + key + '\n')

###############################################################################################

min_graph = gt.Graph(directed=True)
src_trg_pair = []

for key in cg_path.keys():
	src_trg_pair.append((cg_path[key][0], key, cg_path[key][0]+key))
	src_trg_pair.append((cg_path[key][1], key, cg_path[key][1]+key))

vids = min_graph.add_edge_list(src_trg_pair, hashed=True, hash_type="string", eprops=[('name', 'string')])

gt.graph_draw(min_graph, vertex_size=1.5, vertex_text=vids, output=workdir+"/frontal_hierarchical_graph_study.pdf")
gt.interactive_window(min_graph, vertex_size=3, vertex_text=vids)
