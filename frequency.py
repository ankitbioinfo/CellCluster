

f=open('GPMotif')

cont=f.readlines()

d={}
pos=[]
for i in range(len(cont)):
	l=cont[i].split()
	if len(l)==2:
		if l[0]=='GP':
			pos.append(i)
		else:
			d[l[0]]=1
	
pos.append(len(cont))

total=[]
for i in range(len(pos)-1):
	ld={}

	for key in d:
		ld[key]=0

	for j in range(pos[i],pos[i+1]):
		l=cont[j].split()
		if len(l)==2:
			ld[l[0]]=int(l[1])
	total.append(ld)


for key in d:
	

	a=[]
	for i in range(len(total)):
		a.append(total[i][key])

	if (sum(a)==3):
		print(key, a)



	
