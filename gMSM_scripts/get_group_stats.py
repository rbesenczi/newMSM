import numpy
import nibabel
import sys, getopt
import csv
from os import listdir, environ
from os.path import isfile, join, splitext

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
group_subs_lists = home + "/groupwise/data/frontal_subject_clusters_" + dataset + ".csv"

mask_path = home + "/groupwise/NODE2218_frontal_mask.shape.gii"
percentile = 75

groups = list(set([row['group'] for row in csv.DictReader(open(group_list, "r", newline=''), fieldnames=['group','size'])]))
#groups = ["NODE1750", "NODE1807"] #for testing

mask = nibabel.load(mask_path).darrays[0].data

def dice_overlap(subject1, subject2, perc=75):
	subject1perc = numpy.where((subject1 > numpy.percentile(subject1, perc)), 1, 0)
	subject2perc = numpy.where((subject2 > numpy.percentile(subject2, perc)), 1, 0)
	return (2 * numpy.sum(subject1perc * subject2perc)) / (numpy.sum(subject1perc) + numpy.sum(subject2perc))

for group_id in groups:
	with open(group_subs_lists, "r", newline='') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=['line','subject','group'])
		
		subjects = []
		for row in reader:
			if row['group'] == group_id:
				subjects.append(row['subject'])

		num_pairs = (len(subjects)*(len(subjects)-1))/2
		
		corrsum_group = 0.0
		dice_sum_group = 0.0
		corrsum_global = 0.0
		dice_sum_global = 0.0

		for i in range(len(subjects)):
			for j in range(i+1, len(subjects)):

				subject1_group = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[0].data * mask
				subject2_group = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[0].data * mask

				subject1_global = nibabel.load(global_reg + "/output/" + subjects[i] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data * mask
				subject2_global = nibabel.load(global_reg + "/output/" + subjects[j] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data * mask
				
				corrsum_group += numpy.corrcoef(subject1_group, subject2_group)[0,1]
				corrsum_global += numpy.corrcoef(subject1_global, subject2_global)[0,1]

				dice_sum_group += dice_overlap(subject1_group, subject2_group)
				dice_sum_global += dice_overlap(subject1_global, subject2_global)

		areal_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		areal_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		print("Groupwise - Group {}: CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(group_id.split('_')[0], corrsum_group/num_pairs, (dice_sum_group/num_pairs)*100, numpy.mean(areal_group), numpy.max(areal_group), numpy.percentile(areal_group, 95), numpy.percentile(areal_group, 98), numpy.mean(shape_group), numpy.max(shape_group)))
		print("HCP_templ - Group {}: CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(group_id.split('_')[0], corrsum_global/num_pairs, (dice_sum_global/num_pairs)*100, numpy.mean(areal_global), numpy.max(areal_global), numpy.percentile(areal_global, 95), numpy.percentile(areal_global, 98), numpy.mean(shape_global), numpy.max(shape_global)))
		print('\n')
