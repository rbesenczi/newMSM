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

	#subjects = list(set(subjects))
	group_95 = []
	global_95 = []

	for subject in subjects:
		try:
			group_subject_file = nibabel.load(group_reg + "/output/" + root_node_id + "/groupwise." + root_node_id + ".sphere-" + subject + ".distortion.func.gii")
		except:
			#print('file not found')
			continue
		group_subject_areal = group_subject_file.darrays[0].data * mask
		group_subject_shape = group_subject_file.darrays[1].data * mask

		global_subject_file = nibabel.load(global_reg + "/output/" + subject +  ".MSMSulc.ico6.sphere.distortion.func.gii")
		global_subject_areal = global_subject_file.darrays[0].data * mask
		global_subject_shape = global_subject_file.darrays[1].data * mask

		group_95.append(numpy.percentile(group_subject_areal, 95))
		global_95.append(numpy.percentile(global_subject_areal, 95))

		#print("Groupwise - subject {}: Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(subject, numpy.mean(group_subject_areal), numpy.max(group_subject_areal), numpy.percentile(group_subject_areal, 95), numpy.percentile(group_subject_areal, 98), numpy.mean(group_subject_shape), numpy.max(group_subject_shape)))
		#print("Globalreg - subject {}: Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(subject, numpy.mean(global_subject_areal), numpy.max(global_subject_areal), numpy.percentile(global_subject_areal, 95), numpy.percentile(global_subject_areal, 98), numpy.mean(global_subject_shape), numpy.max(global_subject_shape)))
		
		#subject_group = file.darrays[0].data * mask
		#subject_global = nibabel.load(global_reg + "/output/" + subject + ".MSMSulc.ico6.sphere.distortion.func.gii").darrays[0].data * mask

		#print(subject, numpy.mean(subject_group), numpy.mean(subject_global))

		filewriter.writerow([subject, numpy.percentile(group_subject_areal, 95), numpy.percentile(global_subject_areal, 95)])

kstest_result = stats.kstest(group_95, global_95)
print(kstest_result)

#plt.title("Line Graph")
#plt.xlabel("X Axis")
#plt.ylabel("Y Axis")

#plt.plot(list(range(group_95)), group_95, color="red")  # note a[0] instead of a
#plt.show()