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
global_reg = home + "/" + dataset + "Sulc_Curv_to_template"
group_list = group_reg + "/group_list.txt"
group_subs_lists = home + "/groupwise/" + dataset + "/frontal_subject_clusters_" + dataset + ".csv"

stats_file = open(group_reg + "/dist.csv", "a")
stats_writer = csv.writer(stats_file, delimiter=',')
stats_writer.writerow(["Group","Size","Lambda","CC similarity","Dice","Areal 95%"])

mask_path = home + "/groupwise/NODE2218_frontal_mask.shape.gii"
percentile = 75

#groups = [row['group'] for row in csv.DictReader(open(group_list, "r", newline=''), fieldnames=['group','size'])]
groups = ["NODE2078","NODE2152","NODE2162","NODE2139","NODE2158","NODE2146"]

mask = nibabel.load(mask_path).darrays[0].data

def dice_overlap(subject1, subject2, perc=75):
	subject1perc = numpy.where((subject1 > numpy.percentile(subject1, perc)), 1, 0)
	subject2perc = numpy.where((subject2 > numpy.percentile(subject2, perc)), 1, 0)
	return (2 * numpy.sum(subject1perc * subject2perc)) / (numpy.sum(subject1perc) + numpy.sum(subject2perc))

for group_id in groups:
	with open(group_subs_lists, "r", newline='') as csvfile:
		reader = csv.DictReader(csvfile, fieldnames=['line','subject','group'])
		
		subjects = []
		size = 0
		for row in reader:
			if row['group'] == group_id:
				subjects.append(row['subject'])
				size += 1

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

				subject1_group_sulc = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[0].data * mask
				subject2_group_sulc = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[0].data * mask

				subject1_global_sulc = nibabel.load(global_reg + "/output/" + subjects[i] + ".MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii").darrays[0].data * mask
				subject2_global_sulc = nibabel.load(global_reg + "/output/" + subjects[j] + ".MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii").darrays[0].data * mask
				
				corrsum_group_sulc += numpy.corrcoef(subject1_group_sulc, subject2_group_sulc)[0,1]
				corrsum_global_sulc += numpy.corrcoef(subject1_global_sulc, subject2_global_sulc)[0,1]

				subject1_group_curv = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[i] + ".func.gii").darrays[1].data * mask
				subject2_group_curv = nibabel.load(group_reg + "/output/" + group_id + "/groupwise." + group_id + ".transformed_and_reprojected.dedrift-" + subjects[j] + ".func.gii").darrays[1].data * mask

				subject1_global_curv = nibabel.load(global_reg + "/output/" + subjects[i] + ".MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii").darrays[1].data * mask
				subject2_global_curv = nibabel.load(global_reg + "/output/" + subjects[j] + ".MSMSulc.Curv.ico6.transformed_and_reprojected.func.gii").darrays[1].data * mask
				
				corrsum_group_curv += numpy.corrcoef(subject1_group_curv, subject2_group_curv)[0,1]
				corrsum_global_curv += numpy.corrcoef(subject1_global_curv, subject2_global_curv)[0,1]

				dice_sum_group_sulc += dice_overlap(subject1_group_sulc, subject2_group_sulc)
				dice_sum_global_sulc += dice_overlap(subject1_global_sulc, subject2_global_sulc)
				dice_sum_group_curv += dice_overlap(subject1_group_curv, subject2_group_curv)
				dice_sum_global_curv += dice_overlap(subject1_global_curv, subject2_global_curv)

		areal_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape_group = numpy.abs((numpy.asarray(nibabel.load(group_reg + "/results/" + group_id + "/groupwise." + group_id + ".areal.distortion.merge.sulc.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		areal_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".areal.distortion.merge.curv.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())
		shape_global = numpy.abs((numpy.asarray(nibabel.load(global_reg + "/results/" + group_id + ".shape.distortion.merge.curv.affine.dedrifted.ico6.shape.gii").agg_data()) * mask).flatten())

		print("Group {} with {} subjects".format(group_id.split('_')[0], num_subjects))
		print("\tSulc")
		print("\t\tCC similarity: {:.4}; Dice overlap: {:0,.2f}%".format(corrsum_group_sulc/num_pairs, (dice_sum_group_sulc/num_pairs)*100))
		print("\t\tCC similarity: {:.4}; Dice overlap: {:0,.2f}%".format(corrsum_global_sulc/num_pairs, (dice_sum_global_sulc/num_pairs)*100))
		print("\tCurv")
		print("\t\tCC similarity: {:.4}; Dice overlap: {:0,.2f}%".format(corrsum_group_curv/num_pairs, (dice_sum_group_curv/num_pairs)*100))
		print("\t\tCC similarity: {:.4}; Dice overlap: {:0,.2f}%".format(corrsum_global_curv/num_pairs, (dice_sum_global_curv/num_pairs)*100))
		print("\tDistortion")
		print("\t\tAreal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(numpy.mean(areal_group), numpy.max(areal_group), numpy.percentile(areal_group, 95), numpy.percentile(areal_group, 98), numpy.mean(shape_group), numpy.max(shape_group)))
		print("\t\tAreal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(numpy.mean(areal_global), numpy.max(areal_global), numpy.percentile(areal_global, 95), numpy.percentile(areal_global, 98), numpy.mean(shape_global), numpy.max(shape_global)))
		print('\n')

		stats_writer.writerow([group_id,size,"0.6",corrsum_group_sulc/num_pairs,(dice_sum_group_sulc/num_pairs)*100,numpy.percentile(areal_group, 95)])

stats_file.close()
