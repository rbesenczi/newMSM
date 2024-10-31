import numpy
import nibabel
import sys, getopt
import csv
from os import listdir, environ
from os.path import isfile, join, splitext
home = environ['HOME']

###########################################################
##################### SET BY THE USER #####################
group_reg = home + "/groupwise/HCP"
global_reg = home + "/HCP_to_template"
group_subs_lists = home + "/groupwise/HCP/frontal_subject_clusters_HCP.csv"

group_id = "NODE2078"

percentile = 75 #thresholding for DICE overlap
###########################################################

def dice_overlap(subject1, subject2, perc=75):
	subject1perc = numpy.where((subject1 > numpy.percentile(subject1, perc)), 1, 0)
	subject2perc = numpy.where((subject2 > numpy.percentile(subject2, perc)), 1, 0)
	return (2 * numpy.sum(subject1perc * subject2perc)) / (numpy.sum(subject1perc) + numpy.sum(subject2perc))

with open(group_subs_lists, "r", newline='') as csvfile:
	reader = csv.DictReader(csvfile, fieldnames=['line','subject','group'])
	
	subjects = []
	for row in reader:
		if row['group'] == group_id:
			subjects.append(row['subject'])

	num_subjects = len(subjects)
	num_pairs = (num_subjects*(num_subjects-1))/2
	
	corrsum_group_sulc = 0.0
	corrsum_group_curv = 0.0
	dice_sum_group_sulc = 0.0
	dice_sum_group_curv = 0.0
	corrsum_global_sulc = 0.0
	corrsum_global_curv = 0.0
	dice_sum_global_sulc = 0.0
	dice_sum_global_curv = 0.0

	for i in range(len(subjects)):
		for j in range(i+1, len(subjects)):

			subject1_group_sulc = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[0].data
			subject2_group_sulc = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[0].data

			subject1_global_sulc = nibabel.load(global_reg + "/output/" + subjects[i] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data
			subject2_global_sulc = nibabel.load(global_reg + "/output/" + subjects[j] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data
			
			corrsum_group_sulc += numpy.corrcoef(subject1_group_sulc, subject2_group_sulc)[0,1]
			corrsum_global_sulc += numpy.corrcoef(subject1_global_sulc, subject2_global_sulc)[0,1]

			subject1_group_curv = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[1].data
			subject2_group_curv = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[1].data

			subject1_global_curv = nibabel.load(global_reg + "/output/" + subjects[i] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[1].data
			subject2_global_curv = nibabel.load(global_reg + "/output/" + subjects[j] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[1].data
			
			corrsum_group_curv += numpy.corrcoef(subject1_group_curv, subject2_group_curv)[0,1]
			corrsum_global_curv += numpy.corrcoef(subject1_global_curv, subject2_global_curv)[0,1]

			dice_sum_group_sulc += dice_overlap(subject1_group_sulc, subject2_group_sulc)
			dice_sum_global_sulc += dice_overlap(subject1_global_sulc, subject2_global_sulc)
			dice_sum_group_curv += dice_overlap(subject1_group_curv, subject2_group_curv)
			dice_sum_global_curv += dice_overlap(subject1_global_curv, subject2_global_curv)

	areal_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data())).flatten())
	shape_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data())).flatten())

	areal_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".areal.distortion.merge.curv.affine.dedrifted.ico6.shape.gii").agg_data())).flatten())
	shape_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".shape.distortion.merge.curv.affine.dedrifted.ico6.shape.gii").agg_data())).flatten())

	cc_sim_group_sulc = corrsum_group_sulc/num_pairs
	cc_sim_global_sulc = corrsum_global_sulc/num_pairs

	d_sum_group_sulc = dice_sum_group_sulc/num_pairs
	d_sum_global_sulc = dice_sum_global_sulc/num_pairs

	cc_sim_group_curv = corrsum_group_curv/num_pairs
	cc_sim_global_curv = corrsum_global_curv/num_pairs

	d_sum_group_curv = dice_sum_group_curv/num_pairs
	d_sum_global_curv = dice_sum_global_curv/num_pairs

	areal_group_mean = numpy.mean(areal_group)
	areal_global_mean = numpy.mean(areal_global)

	areal_group_max = numpy.max(areal_group)
	areal_global_max = numpy.max(areal_global)

	areal_group_95 = numpy.percentile(areal_group, 95)
	areal_global_95 = numpy.percentile(areal_global, 95)

	areal_group_98 = numpy.percentile(areal_group, 98)
	areal_global_98 = numpy.percentile(areal_global, 98)

	shape_group_mean = numpy.mean(shape_group)
	shape_global_mean = numpy.mean(shape_global)

	shape_group_max = numpy.max(shape_group)
	shape_global_max = numpy.max(shape_global)

	print("\tStats for group " + group_id)
	print("\tSulc")
	print("\t\tCC similarity: {:.4}; Dice overlap: {:.4}".format(cc_sim_group_sulc, d_sum_group_sulc))
	print("\t\tCC similarity: {:.4}; Dice overlap: {:.4}".format(cc_sim_global_sulc, d_sum_global_sulc))
	print("\tCurv")
	print("\t\tCC similarity: {:.4}; Dice overlap: {:.4}".format(cc_sim_group_curv, d_sum_group_curv))
	print("\t\tCC similarity: {:.4}; Dice overlap: {:.4}".format(cc_sim_global_curv, d_sum_global_curv))
	print("\tDistortion")
	print("\t\tAreal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(areal_group_mean, areal_group_max, areal_group_95, areal_group_98, shape_group_mean, shape_group_max))
	print("\t\tAreal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(areal_global_mean, areal_global_max, areal_global_95, areal_global_98, shape_global_mean, shape_global_max))
	print('\n')
