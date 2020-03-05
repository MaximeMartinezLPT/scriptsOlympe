import sys
import os
sys.path.insert(0, '..')
sys.path.insert(0, '../..')

from inputparams import *
from dynamics1D.potential import *
from dynamics1D.quantum import *

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
	h=float(params['h'])
	phi=float(params['phi'])
	Npcell=int(params['Npcell'])
	nruns=int(params['nruns'])
	
	beta0=np.linspace(0,1,nruns,endpoint=False)
	beta0=np.sort(beta0*(beta0<=0.5)+(beta0-1)*(beta0>0.5))
	beta=beta0[runid]

	grid=Grid(Npcell,h)
	pot=ModulatedPendulum(e,gamma,phi=phi)
	floquet=FloquetPropagator(grid,pot,beta=beta,T0=4*np.pi,idtmax=1000)

	wf=WaveFunction(grid)
	wf.setState("coherent",2.0)
	
	floquet.diagonalize()
	
	overlaps=floquet.orderEigenstatesWithOverlapOn(wf)
	quasienergies=floquet.quasienergy
	
	evec0p=floquet.eigenvec[0].p
	
	np.savez(wdir+"dataruns/"+str(runid),"w", quasienergies=quasienergies, overlaps=overlaps,evec0p=evec0p)
	
if mode=="gather":
	params=readInput(wdir+"params.txt")
	Npcell=int(params['Npcell'])
	nruns=int(params['nruns'])
	h=float(params['h'])
	
	quasienergies=np.zeros((nruns,Npcell))
	overlaps=np.zeros((nruns,Npcell),dtype=complex)
	beta=np.zeros((nruns,Npcell))
	
	beta0=np.linspace(0,1,nruns,endpoint=False)
	beta0=np.sort(beta0*(beta0<=0.5)+(beta0-1)*(beta0>0.5))
	
	Ncell=nruns
	
	grid=Grid(Ncell*Npcell,h,xmax=Ncell*2*np.pi)
	wannier=WaveFunction(grid)
	ind=np.flipud(np.arange(0,Ncell*Npcell,Ncell))
	
	for irun in range(nruns):
		beta[irun]=np.ones(Npcell)*np.sort(beta0*(beta0<=0.5)+(beta0-1)*(beta0>0.5))[irun]
		
		data=np.load(wdir+"dataruns/"+str(irun)+".npz")
		quasienergies[irun]=data['quasienergies']
		overlaps[irun]=data['overlaps']
		wannier.p[ind]=np.abs(data['evec0p'])
		data.close()

		ind=ind+1

	wannier.p=np.roll(wannier.p,int(Ncell/2+1))
	wannier.p=np.fft.fftshift(wannier.p*grid.phaseshift)

	wannier.normalize("p")


	np.savez(wdir+"data","w", quasienergies=quasienergies, overlaps=overlaps,beta=beta,wannierx=wannier.x,gridx=grid.x)	
	
	os.system("rm -r "+wdir+"dataruns/")
	
if mode=="plot":
	
	data=np.load(wdir+"data.npz")
	quasienergies=data['quasienergies']
	overlaps=data['overlaps']
	beta=data['beta']
	gridx=data['gridx']
	wannierx=data['wannierx']
	
	ax=plt.subplot(1,2,1)
	
	
	# ~ ax.set_title(r"$\varepsilon={:.3f} \quad \gamma={:.3f} \quad 1/h={:.3f} \quad N_p={:d}$".format(float(e),float(gamma),float(1/h),int(beta.size)))
	ax.set_xlabel(r"$\beta$")
	ax.set_ylabel(r"$qE/h$")
	ax.set_xlim(np.min(beta),np.max(beta))
	
	# ~ print(np.min(np.mean(quasienergies[:,0])),1.2*np.max(np.mean(quasienergies[:,0])))
	# ~ ax.set_ylim(3*np.min(np.mean(quasienergies[:,0])),3*np.max(np.mean(quasienergies[:,0])))		

	ax.scatter(beta[:,0],quasienergies[:,0],c="red",s=5**2,zorder=2)
	
	ax=plt.subplot(1,2,2)
	
	ax.plot(gridx,np.abs(wannierx))
	
	ax.set_xlim(-5*np.pi,5*np.pi)
	
	# ~ plt.show()

	plt.savefig(wdir+"spectrum.png", bbox_inches='tight',dpi=250)
	
	
	
