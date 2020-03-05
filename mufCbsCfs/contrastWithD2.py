import sys
import os
sys.path.insert(0, '..')
sys.path.insert(0, '../..')

from inputparams import *
from dynamics1D.potential import *
from dynamics1D.quantum import *

mode=sys.argv[1]
wdir=sys.argv[2]

def naverage(N):
		p=np.log(N)/np.log(2)
		dp=p-pNmin
		return int(2**(10-dp))

if mode=="initialize":
	os.mkdir(wdir)
	os.mkdir(wdir+"dataruns")
	os.mkdir(wdir+"d2control")
	
	inputfile=sys.argv[3]+".txt"
	os.system("cp "+inputfile+" "+wdir+"params.txt")
	nruns=int(sys.argv[4])
	addParams(wdir+"params.txt",{'nruns':nruns})
	addParams(wdir+"params.txt",{'nN':5})
	addParams(wdir+"params.txt",{'nalpha':int(nruns/5)})
	addParams(wdir+"params.txt",{'pNmin':7})

if mode=="compute":
	# Loading input file
	
	runid=int(sys.argv[3])-1
	
	params=readInput(wdir+"params.txt")
	potential=params['potential']
	potential=params['potential']
	nalpha=int(params['nalpha'])
	pNmin=int(params['pNmin'])
	nN=int(params['nN'])
	
	N=int(np.tile(2**np.arange(pNmin,pNmin+nN),nalpha)[runid])

	if potential=="Rectangle":
		alpha=np.repeat(np.linspace(0,0.5,nalpha),nN)[runid]
		pot=Rectangle(alpha)
		i0=int(N/8)
	
	if potential=="SawTooth":
		alpha=np.repeat(np.linspace(0,1.0,nalpha),nN)[runid]
		pot=SawTooth(alpha)
		i0=int(N/4)	
	
	momenta=0
	phi2=np.zeros(N)
	phi4=np.zeros(N)

	for i in range(naverage(N)):
		grid=Grid(N,h=1,xmax=2*np.pi)
	
		fo=FloquetRandomPhasePropagator(grid,pot,T0=1,idtmax=1)
		fo.diagonalize()
		
		momenta+=np.sum(np.array([fo.eigenvec[i].getMomentum("p",2) for i in range(N)]))
		
		phi2+=np.array([np.abs(fo.eigenvec[j].x[i0])**2 for j in range(N)])/naverage(N) # axis 1 otherwise
		phi4+=np.array([np.abs(fo.eigenvec[j].x[i0])**4 for j in range(N)])/naverage(N)

	momenta/=naverage(N)*N
	Cinf=np.sum(phi4)/np.sum(phi2**2)-1
		
	np.savez(wdir+"dataruns/"+str(runid),momenta=momenta,alpha=alpha,N=N,Cinf=Cinf)
	
if mode=="gather":
	params=readInput(wdir+"params.txt")
	nalpha=int(params['nalpha'])
	pNmin=int(params['pNmin'])
	nN=int(params['nN'])
	
	D2=np.zeros(nalpha)
	alpha=np.zeros(nalpha)
	Cinf=np.zeros((nN,nalpha))
	

	for ialpha in range(nalpha):
		N=np.zeros(nN)
		momenta=np.zeros(nN)
		for iN in range(nN):
			runid=ialpha*nN+iN
			data=np.load(wdir+"dataruns/"+str(runid)+".npz")
			N[iN]=data['N']
			momenta[iN]=data['momenta']
			alpha[ialpha]=data['alpha']
			Cinf[iN,ialpha]=data['Cinf']
			data.close()
			
		fit = np.polyfit(np.log(N),np.log(momenta), 1)
		D2[ialpha]=-fit[0]
		print(alpha[ialpha],D2[ialpha])
			
		ax=plt.gca()
		ax.set_title(r"$alpha={:.2f}$".format(alpha[ialpha]))
			
		plt.scatter(np.log(N),np.log(momenta))
		plt.plot(np.log(N),fit[0]*np.log(N)+fit[1],c="red")
		plt.savefig(wdir+"d2control/"+str(ialpha)+".png", bbox_inches = 'tight',format="png")
		plt.clf()
		
	np.savez(wdir+"data",alpha=alpha,D2=D2,Cinf=Cinf)
	os.system("rm -r "+wdir+"dataruns/")
		
if mode=="plot":
	data=np.load(wdir+"data.npz")
	alpha=data['alpha']
	D2=data['D2']
	Cinf=data['Cinf']
	data.close()
	
	ax=plt.subplot(1,2,1)
	ax.set_xlim(np.min(alpha),np.max(alpha))
	ax.set_ylim(0.0,1.0)
	ax.set_xlabel(r"$\alpha$")
	ax.set_ylabel(r"$D_2$")
	
	# ~ D2th=(2*alpha*sc.gamma(2-0.5*np.ones(alpha.size)))/(np.sqrt(np.pi)*sc.gamma(2))
	
	ax.grid()
	
	plt.scatter(alpha,D2,c="red")
	
	ax=plt.subplot(1,2,2)
	for i in range(5):
		ax.scatter(D2,Cinf[i])
		
	# ~ ax.set_ylim(0.0,1.0)
	ax.set_xlim(0.0,np.max(D2))
	ax.set_xlabel(r"$D_2$")
	ax.set_ylabel(r"$\Lambda_\infty=\frac{\sum_n \langle |\phi_n(x_0)|^4 \rangle}{\sum_n \langle |\phi_n(x_0)|^2 \rangle^2}$")
	ax.grid()
	
	# ~ plt.show()

	plt.savefig(wdir+"d2.png", bbox_inches = 'tight',format="png")



				
	


