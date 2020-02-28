import sys

def readInput(dfile):
	paramsinput={}
	with open(dfile) as f:
		for line in f:
			(key, val) = line.split()
			paramsinput[key] = val
	return paramsinput
		
			
def writeInput(dfile,params):
	f = open(dfile,"w")
	for p in params:
		f.write(str(p)+" "+str(params[p])+"\n")
	f.close()
	
def addParams(dfile,params):
	f = open(dfile,"a+")
	for p in params:
		f.write(str(p)+" "+str(params[p])+"\n")
	f.close()
