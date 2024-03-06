import numpy
import nibabel
import sys
import csv
from os import listdir
from os.path import isfile, join, splitext

dataset = "HCP"
workdir = "/home/rb22/groupwise/" + dataset + "/"
group_lists = "group_lists"

percentile = 75

list_path = workdir + group_lists
onlyfiles = [f for f in listdir(list_path) if isfile(join(list_path, f))]

results = open(workdir + "results.csv", "w")
results.write("Group ID (numsubjects),CC similarity,Dice overlap,Areal mean,Areal Max,Areal 95%,Areal 98%,Shape mean,Shape Max\n")

mask = nibabel.load("/home/rb22/groupwise/NODE2218_frontal_mask.shape.gii").darrays[0].data

for file in onlyfiles:
	with open(list_path + "/" + file, "r", newline='') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=['subject','group'])
		
		group_id = splitext(file)[0]
		subjects = []
		for row in reader:
			subjects.append(row['subject'])

		num_pairs = (len(subjects)*(len(subjects)-1))/2
		corrsum = 0.0
		dice_sum = 0.0

		for i in range(len(subjects)):
			for j in range(i+1, len(subjects)):

				sub1 = nibabel.load(workdir + "output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[0].data
				sub2 = nibabel.load(workdir + "output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[0].data

				subject1 = sub1 * mask
				subject2 = sub2 * mask

				corrsum += numpy.corrcoef(subject1, subject2)[0,1]
				sub1perc = numpy.where((subject1 > numpy.percentile(subject1, percentile)), 1, 0)
				sub2perc = numpy.where((subject2 > numpy.percentile(subject2, percentile)), 1, 0)
				dice_sum += (2 * numpy.sum(sub1perc * sub2perc)) / (numpy.sum(sub1perc) + numpy.sum(sub2perc))

		areal = numpy.abs((numpy.asarray(nibabel.load(workdir + "results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape = numpy.abs((numpy.asarray(nibabel.load(workdir + "results/" + group_id + "/groupwise." + group_id + ".shape.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		print("Group {}: CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(group_id.split('_')[0], corrsum/num_pairs, (dice_sum/num_pairs)*100, numpy.mean(areal), numpy.max(areal), numpy.percentile(areal, 95), numpy.percentile(areal, 98), numpy.mean(shape), numpy.max(shape)))
		results.write("{} ({}),{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}\n".format(group_id.split('_')[0], len(subjects), corrsum/num_pairs, dice_sum/num_pairs, numpy.mean(areal), numpy.max(areal), numpy.percentile(areal, 95), numpy.percentile(areal, 98), numpy.mean(shape), numpy.max(shape)))
