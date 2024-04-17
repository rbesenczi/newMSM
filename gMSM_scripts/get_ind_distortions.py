import numpy
import nibabel
import sys, getopt
import csv
from os import listdir, environ
from os.path import isfile, join, splitext
from scipy import stats
import matplotlib.pyplot as plt

dataset = ""

try:
	opts, args = getopt.getopt(sys.argv[1:], "h:d:", ["help","dataset="])
except getopt.GetoptError:
	print("extract_info.py -d {HCP,UKB}")
	sys.exit(2)
for o, a in opts:
	if o in ("-h", "--help"):
		print("extract_info.py -d {HCP,UKB}")
		sys.exit()
	if o in ("-d", "--dataset"):
		dataset = a

home = environ['HOME']
group_reg = home + "/groupwise/" + dataset
global_reg = home + "/" + dataset + "_to_template"
group_list = group_reg + "/group_list.txt"
group_subs_lists = home + "/groupwise/data/frontal_subject_clusters_" + dataset + "_base.csv"

mask_path = home + "/groupwise/NODE2218_frontal_mask.shape.gii"
percentile = 75

root_node_id = "NODE2218"

mask = nibabel.load(mask_path).darrays[0].data

with open(group_subs_lists, "r", newline='') as csvfile:
	file = open(group_reg+"/stats.csv", "w")
	filewriter = csv.writer(file)
	filewriter.writerow(['subject', 'group', 'global'])

	reader = csv.DictReader(csvfile, fieldnames=['line','subject','group'])
	
	subjects = []
	for row in reader:
		subjects.append(row['subject'])

	group_95 = []
	global_95 = []

	for subject in subjects:
		try:
			group_subject_file = nibabel.load(group_reg + "/output/" + root_node_id + "/groupwise." + root_node_id + ".sphere-" + subject + ".distortion.func.gii")
		except:
			continue
		group_subject_areal = group_subject_file.darrays[0].data * mask
		group_subject_shape = group_subject_file.darrays[1].data * mask

		global_subject_file = nibabel.load(global_reg + "/output/" + subject +  ".MSMSulc.ico6.sphere.distortion.func.gii")
		global_subject_areal = global_subject_file.darrays[0].data * mask
		global_subject_shape = global_subject_file.darrays[1].data * mask

		group_95.append(numpy.percentile(group_subject_areal, 95))
		global_95.append(numpy.percentile(global_subject_areal, 95))

		filewriter.writerow([subject, numpy.percentile(group_subject_areal, 95), numpy.percentile(global_subject_areal, 95)])

kstest_result = stats.kstest(group_95, global_95)
print(kstest_result)
