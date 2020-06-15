
import math 
import numpy as np 
import networkx as nx 
import scipy.io as sio 
import matplotlib.pyplot as plt
import os

'''
maindir='Columnar_Structure_Prediction_8_100'

answer=os.path.isdir(maindir)
if answer==True:
	pass
else:
	os.mkdir(maindir)
'''



dt_path_wt=[ '../data/Nuclei_and_Cells_DT_S18_m6_wt/', '../data/Nuclei_and_Cells_DT_S17_m2_wt/',\
			 '../data/Nuclei_and_Cells_DT_S84_m3_wt/', '../data/Nuclei_and_Cells_DT_S51_m2_wt/',\
			 '../data/Nuclei_and_Cells_DT_S84_m4_wt/'];

pt_path_wt = [  '../data/Nuclei_and_Cells_PT_S18_m6_wt/','../data/Nuclei_and_Cells_PT_S17_m2_wt/',\
				'../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/',\
				'../data/Nuclei_and_Cells_PT_S84_m4_wt/'];

dt_path_mut= ['../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,\
			  '../data/Nuclei_and_Cells_DT_S84_m1_mut/', '../data/Nuclei_and_Cells_DT_S84_m5_mut/'];
	
pt_path_mut = ['../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',\
			   '../data/Nuclei_and_Cells_PT_S84_m1_mut/', '../data/Nuclei_and_Cells_PT_S84_m5_mut/', ];
			   
du_path_wt=['../data/Nuclei_and_Cells_DU_S51_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m3_wt/'];

	   



def read_matlab_files(filedata):
	d={}
	d2={}
	d3={}
	count=0
	for line in filedata:
		count+=1
		l=line.split()
		a=[]
		for j in range(len(l)):
			a.append(int(l[j]))

		a=sorted(a)
		name=''
		for i in a:
			name+=str(i)+'-'

		if name in d:
			d[name]+=1
			d3[name].append(count)
		else:
			d[name]=1
			d3[name]=[]
			d3[name].append(count)

	
		d2[name]=count



	return [d,d2,d3] 



allpath=[dt_path_wt, pt_path_wt, dt_path_mut, pt_path_mut, du_path_wt]; 


for i in range(0,len(allpath)):

	growthPlate={}
	for j in range(len(allpath[i])):
		path=allpath[i][j]
		mypath=path[0:-1].split('Nuclei_and_Cells_')	
		#d=read_matlab_files('MakeListColumnarStructurePrediction/'+mypath[1]+'/')
		mypath=path[0:-1].split('Nuclei_and_Cells_')

		f=open('degree_sequence/degree_'+mypath[1]+'.dat')
		[d,d2,d3]=read_matlab_files(f)
	

		#print(mypath[1])
		a=[]
		for key in d:
			#print(key, '\t',d[key],'\t\t', d2[key])
			if key in growthPlate:	
				growthPlate[key]+=d[key]
			else:
				growthPlate[key]=d[key]
			a.append(d2[key])

		#print(a)
		#print('\n')


		try:
			a=[]
			a=a+d3['1-1-2-']
			a=a+d3['1-1-2-2-']
			a=a+d3['1-1-2-2-2-']
			a=a+d3['1-1-2-2-2-2-']
			a=a+d3['1-1-2-2-2-2-2-']
			a=a+d3['1-1-2-2-2-2-2-2-']
		except KeyError:
			pass
		


		#newname=maindir+'/'+mypath[1]
		newname='./StraightLine/'
		answer=os.path.isdir(newname)
		if answer==True:
			pass
		else:
			os.mkdir(newname)
		

		f=open(newname+'pureColum_'+mypath[1]+'.dat','w')
		for k in a:
			f.write(str(k)+'\n')
		
	
	print('\nGP', i+1,'\n')
	for key in growthPlate:
		print(key,'\t',growthPlate[key])
		


