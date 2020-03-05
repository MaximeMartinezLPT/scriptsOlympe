import sys
import os
sys.path.insert(0, '..')
sys.path.insert(0, '../..')

from inputparams import *
from dynamics1D.potential import *
from dynamics1D.classical import *

mode=sys.argv[1]
wdir=sys.argv[2]


if mode=="initialize":
	os.mkdir(wdir)
	os.mkdir(wdir+"dataruns")
	
	inputfile=sys.argv[3]+".txt"
	os.system("mv "+inputfile+" "+wdir+"params.txt")
	nruns=int(sys.argv[4])
	addParams(wdir+"params.txt",{'nruns':nruns})
	
if mode=="compute":
	runid=int(sys.argv[3])-1
	
	params=readInput(wdir+"params.txt")
	gamma=float(params['gamma'])
	e=float(params['epsilon'])
	tmax=int(params['tmax'])
	nruns=int(params['nruns'])
	
	icg=InitialConditionGenerator(nruns)
	pot=DoubleWell(e,gamma)
	
	tp=ContinuousTimePropagator(pot,T0=0.01,ndt=10) #Propagateur entre deux points du portrait de phase
	pp=PhasePortrait(tmax,nruns,tp,icg.generateXP)
	x,p=pp.computeOrbit(runid)
	c=pp.getChaoticity(x,p)

	np.savez(wdir+"dataruns/"+str(runid),"w", x=modc(x,2*np.pi),p=p,c=c)
	
if mode=="gather":
	params=readInput(wdir+"params.txt")
	nruns=int(params['nruns'])
	tmax=int(params['tmax'])
	
	c=np.zeros((nruns,tmax))
	x=np.zeros((nruns,tmax))
	p=np.zeros((nruns,tmax))

	for i in range(nruns):
		data=np.load(wdir+"dataruns/"+str(i)+".npz")
		x[i]=data['x']
		p[i] = data['p']
		c[i]= data['c']*np.ones(tmax)
		data.close()
		
	np.savez(wdir+"data","w", x=x, p=p,c=c/np.max(c))
	os.system("rm -r "+wdir+"dataruns/")

	
if mode=="plothusimi":
	# Loading inpute file
	fig, ax = plt.subplots(figsize=(np.pi*2,4),frameon=False)
	ax.set_xlim(-np.pi,np.pi)
	ax.set_ylim(-2.0,2.0)
	ax.set_xticks([])
	ax.set_yticks([])
	
	data=np.load(wdir+"data.npz")
	x=data["x"]
	p=data["p"]
	c=data["c"]
	data.close()
	
	
	plt.scatter(x,p,s=0.1**2,c="black")
	
	fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
	plt.savefig(wdir+"PP-husimi.png",dpi=500)
	
if mode=="plot":
	# Loading inpute file
	fig, ax = plt.subplots(figsize=(np.pi*2,4))
	
	ax.set_xlim(-np.pi,np.pi)

	
	ax.set_xlabel(r"Position")
	ax.set_ylabel(r"Vitesse")

	data=np.load(wdir+"data.npz")
	x=data["x"]
	p=data["p"]
	c=data["c"]
	data.close()
	
	condreg=c<0.5	
	condchaotic=c>0.5
	plt.scatter(x,p,s=0.1**2,c="black")
	
	
	plt.savefig(wdir+"PP.png",dpi=500)
	

	
	
	
