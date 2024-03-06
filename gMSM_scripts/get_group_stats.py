import numpy
import nibabel
import sys
import csv
from os import listdir
from os.path import isfile, join, splitext

dataset = "HCP"
workdir = "/home/rb22/groupwise/" + dataset + "/"
global_reg = "/home/rb22/HCP_to_template/"
group_lists = "group_lists"
group_id_s = ["NODE2214"]

percentile = 75

list_path = workdir + group_lists

mask = nibabel.load("/home/rb22/groupwise/NODE2218_frontal_mask.shape.gii").darrays[0].data

for group_id in group_id_s:
	with open(list_path + "/" + group_id + "_HPC.csv", "r", newline='') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=['subject','group'])
		
		subjects = []
		for row in reader:
			subjects.append(row['subject'])

		num_pairs = (len(subjects)*(len(subjects)-1))/2
		corrsum = 0.0
		dice_sum = 0.0
		corrsum_global = 0.0
		dice_sum_global = 0.0

		for i in range(len(subjects)):
			for j in range(i+1, len(subjects)):

				sub1 = nibabel.load(workdir + "output/" + group_id + "_HPC/groupwise." + group_id + "_HPC.transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[0].data
				sub2 = nibabel.load(workdir + "output/" + group_id + "_HPC/groupwise." + group_id + "_HPC.transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[0].data

				sub1_global = nibabel.load(global_reg + "0/" + subjects[i] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data
				sub2_global = nibabel.load(global_reg + "0/" + subjects[j] + ".MSMSulc.ico6.transformed_and_reprojected.func.gii").darrays[0].data
				
				subject1 = sub1 * mask
				subject2 = sub2 * mask
				sub1_global = sub1_global * mask
				sub2_global = sub2_global * mask

				corrsum += numpy.corrcoef(subject1, subject2)[0,1]
				corrsum_global += numpy.corrcoef(sub1_global, sub2_global)[0,1]

				sub1perc = numpy.where((subject1 > numpy.percentile(subject1, percentile)), 1, 0)
				sub2perc = numpy.where((subject2 > numpy.percentile(subject2, percentile)), 1, 0)
				dice_sum += (2 * numpy.sum(sub1perc * sub2perc)) / (numpy.sum(sub1perc) + numpy.sum(sub2perc))

				sub1perc_global = numpy.where((sub1_global > numpy.percentile(sub1_global, percentile)), 1, 0)
				sub2perc_global = numpy.where((sub2_global > numpy.percentile(sub2_global, percentile)), 1, 0)
				dice_sum_global += (2 * numpy.sum(sub1perc_global * sub2perc_global)) / (numpy.sum(sub1perc_global) + numpy.sum(sub2perc_global))

		areal = numpy.abs((numpy.asarray(nibabel.load(workdir + "results/" + group_id + "_HPC/groupwise." + group_id + "_HPC.areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape = numpy.abs((numpy.asarray(nibabel.load(workdir + "results/" + group_id + "_HPC/groupwise." + group_id + "_HPC.shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		areal_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "distortions/" + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "distortions/" + group_id + ".shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		print("Groupwise - Group {}: CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(group_id.split('_')[0], corrsum/num_pairs, (dice_sum/num_pairs)*100, numpy.mean(areal), numpy.max(areal), numpy.percentile(areal, 95), numpy.percentile(areal, 98), numpy.mean(shape), numpy.max(shape)))
		print("HCP_templ - Group {}: CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(group_id.split('_')[0], corrsum_global/num_pairs, (dice_sum_global/num_pairs)*100, numpy.mean(areal_global), numpy.max(areal_global), numpy.percentile(areal_global, 95), numpy.percentile(areal_global, 98), numpy.mean(shape_global), numpy.max(shape_global)))
