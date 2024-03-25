from os import environ
import csv
import copy

home = environ['HOME']
dataset = "HCP"
workdir = home + "/groupwise/" + dataset
clustering = home + "/groupwise/data/frontal_subject_clusters_HCP_base.csv"
path = workdir + "/frontal_hierarchical_path_study.csv"
order = workdir + "/order.txt"

clustering_file = open(clustering, "r", newline='')
path_file = open(path, "r", newline='')
order_file = open(order, "w")

cluster_reader = csv.DictReader(clustering_file, fieldnames=['linenum','subject','group'])
path_reader = csv.DictReader(path_file, fieldnames=['left','right','root'])

clusters_subjects = {}
gen_mean = []

for row in cluster_reader:
	if row['group'] in clusters_subjects.keys():
		clusters_subjects[row['group']].append(row['subject'])
	else:
		clusters_subjects[row['group']] = [ row['subject'] ]

have_mean = []
for key in clusters_subjects.keys():
	have_mean.append(key)

file_counter = 1
block_file = open(workdir + "/blocks/block_1.txt", "w")

for row in path_reader:
	if row['left'] not in have_mean or row['right'] not in have_mean:
		block_file.close()
		file_counter = file_counter + 1
		block_file = open(workdir + "/blocks/block_" + str(file_counter) + ".txt", "w")
		for mean in gen_mean:
			print(mean)
			order_file.write(mean + '\n')
			block_file.write(mean + '\n')
		gen_mean = []
		block_file.close()
		file_counter = file_counter + 1
		block_file = open(workdir + "/blocks/block_" + str(file_counter) + ".txt", "w")
	for subject in clusters_subjects[row['left']]:
		print('0,' + subject + ',' + row['left'] + ',' + row['right'] + ',' + row['root'])
		order_file.write('0,' + subject + ',' + row['left'] + ',' + row['right'] + ',' + row['root'] + '\n')
		block_file.write('0,' + subject + ',' + row['left'] + ',' + row['right'] + ',' + row['root'] + '\n')
		if row['root'] in clusters_subjects.keys():
			clusters_subjects[row['root']].append(subject)
		else:
			clusters_subjects[row['root']] = [ subject ]
	for subject in clusters_subjects[row['right']]:
		print('0,' + subject + ',' + row['right'] + ',' + row['left'] + ',' + row['root'])
		order_file.write('0,' + subject + ',' + row['right'] + ',' + row['left'] + ',' + row['root'] + '\n')
		block_file.write('0,' + subject + ',' + row['right'] + ',' + row['left'] + ',' + row['root'] + '\n')
		if row['root'] in clusters_subjects.keys():
			clusters_subjects[row['root']].append(subject)
		else:
			clusters_subjects[row['root']] = [ subject ]
	gen_mean.append('1,NA,' + row['left'] + ',' + row['right'] + ',' + row['root'])
