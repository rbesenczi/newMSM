import numpy
import nibabel
import sys

###########################################################
## Set the following lines according to your settings and cohort.
workdir = "/home/rb22/groupwise/"

subjects = [ 117122, 117123, 128026, 139637, 561444, 660951, 677968 ]
group_id = 2133

percentile = 75
###########################################################

num_pairs = (len(subjects)*(len(subjects)-1))/2
corrsum = 0.0
dice_sum = 0.0

for i in range(len(subjects)):
	for j in range(i+1, len(subjects)):
		subject1 = nibabel.load(workdir + "output/" + str(group_id) + "/groupwise." + str(group_id) + ".transformed_and_reprojected.dedrift-" + str(subjects[i]) + ".func.gii").darrays[0].data
		subject2 = nibabel.load(workdir + "output/" + str(group_id) + "/groupwise." + str(group_id) + ".transformed_and_reprojected.dedrift-" + str(subjects[j]) + ".func.gii").darrays[0].data
		
		corrsum += numpy.corrcoef(subject1, subject2)[0,1]
		
		sub1perc = numpy.where((subject1 > numpy.percentile(subject1, percentile)), 1, 0)
		sub2perc = numpy.where((subject2 > numpy.percentile(subject2, percentile)), 1, 0)
		dice_sum += (2 * numpy.sum(sub1perc * sub2perc)) / (numpy.sum(sub1perc) + numpy.sum(sub2perc))

areal = numpy.abs(numpy.asarray(nibabel.load(workdir + "results/" + str(group_id) + "/groupwise." + str(group_id) + ".areal.distortion.merge.L.sulc.affine.dedrifted.ico6.shape.gii").agg_data()).flatten())
shape = numpy.abs(numpy.asarray(nibabel.load(workdir + "results/" + str(group_id) + "/groupwise." + str(group_id) + ".shape.distortion.merge.L.sulc.affine.dedrifted.ico6.shape.gii").agg_data()).flatten())

print("CC similarity: {:.4}; Dice overlap: {:0,.2f}%; Areal mean: {:.4}; Areal Max: {:.4}; Areal 95%: {:.4}; Areal 98%: {:.4}; Shape mean: {:.4}; Shape Max: {:.4}".format(corrsum/num_pairs, (dice_sum/num_pairs)*100, numpy.mean(areal), numpy.max(areal), numpy.percentile(areal, 95), numpy.percentile(areal, 98), numpy.mean(shape), numpy.max(shape)))

#This is for a CSV like output
#print("{} ({}),{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}".format(group_id, len(subjects), corrsum/num_pairs, dice_sum/num_pairs, numpy.mean(areal), numpy.max(areal), numpy.percentile(areal, 95), numpy.percentile(areal, 98), numpy.mean(shape), numpy.max(shape)))
