import os
import numpy as np

def mass_hit(mz1,mz2,tol):
    if 1e6*np.abs(mz1-mz2)/mz1 < tol:
        return True
    else:
        return False

def rt_hit(rt1,rt2,tol):
    if np.abs(rt1-rt2) < tol:
        return True
    else:
        return False
    
                
    

if __name__ == '__main__':
	import django
	os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'pimp.settings_dev')
	django.setup()
	from frank.models import Experiment,Peak,SampleFile,FragmentationSet,AnnotationQuery,CandidateAnnotation,Compound
	fs = FragmentationSet.objects.filter(slug='beer-3')
	fs = fs[0]
	peaks = Peak.objects.filter(fragmentation_set = fs,msn_level=1)
	print "Found {} peaks".format(len(peaks))

	to_find = [(314.14,269.69,'both'),
				(540.2707,263.085,'ferulic'),
				(307.1767,547.239,'ferulic'),
				(498.2599,615.862,'ferulic'),
				(265.1545,1100.580,'ferulic'),
				(194.0812,364.448,'ferulic'), 
				(538.2802,240.354,'ethyl'),
				(166.1226,692.788,'ethyl'),
				(232.0638,549.717,'ethyl'),
				(246.1334,410.326,'ethyl'),
				(328.1755,652.517,'ethyl'),]
 
	for mz,rt,prefix in to_find:
		hits = []
		for peak in peaks:
			if mass_hit(mz,float(peak.mass),100) and rt_hit(rt,float(peak.retention_time),20):
				hits.append(peak)
		print hits
		if len(hits) == 1:
			peak = hits[0]
			fname = "simon_data/" + prefix + "_" + str(peak.mass) + "_" + str(peak.retention_time) + ".txt"
			with open(fname,'w') as f:
				child_peaks = Peak.objects.filter(parent_peak = peak)
				for child in child_peaks:
					f.write("" + str(child.mass) + "," + str(child.intensity) + "\n")
				f.write("\n\n")
				ca = CandidateAnnotation.objects.filter(peak=peak)
				if ca:
					for i,c in enumerate(ca):
						f.write("Hit " + str(i) + ":" + c.compound.name + "," + c.compound.formula + "\n")




