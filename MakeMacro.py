



import os 
from os import listdir
from os.path import isfile, join

mypath='PT_S84_m5_mut/Figure_tif/'


onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]


s=mypath.split('/')
print(s[0])


myoutput='/home/student/ankit_code_paul_maps_registration/data/ColumnRelated/only_PZ_column4/Videos/'+s[0]
answer=os.path.isdir(myoutput)
if answer==True:
	pass
else:
	os.mkdir(myoutput)





f=open('MakeVideo'+s[0]+'.ijm','w')

print('# of files ', len(onlyfiles))

for i in range(len(onlyfiles)):

	name=onlyfiles[i]
	t1=name.split('Output_cluster_')
	t=t1[1].split('.tif')

	x=t[0].split('_')

	#if int(x[1])==0:
	if True:

		f.write('open("/home/student/ankit_code_paul_maps_registration/data/ColumnRelated/only_PZ_column4/Cluster_Structure_Prediction/'+mypath+name+'");\n')
		f.write('name=getTitle();\n')
		f.write('run("3D Viewer");\n')
		f.write('call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");\n')
		f.write('call("ij3d.ImageJ3DViewer.add", name, "None", name, "0", "true", "true", "true", "2", "0");\n')
		f.write('call("ij3d.ImageJ3DViewer.record360");\n')
		f.write('selectWindow("Movie");\n')
		f.write('run("AVI... ", "compression=JPEG frame=7 save=/home/student/ankit_code_paul_maps_registration/data/ColumnRelated/only_PZ_column4/Videos/'+s[0]+'/Movie_'+t[0]+'.avi");\n')
		f.write('close();\n')
		f.write('selectWindow(name);\n')
		f.write('close();\n')
		f.write('call("ij3d.ImageJ3DViewer.close");\n')

		f.write('\n\n')


