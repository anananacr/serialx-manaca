import re
import math
import subprocess as sub
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
import csv
import subprocess as sub
from datetime import datetime
from itertools import product
import crystplots
import peakopt
import random

"""
This function generates the geometry file for a detector to sample distance of det_dist added a step of coff
"""
def gen_geom(coff, det_dist):
	f=open("pilatus2mpanel.geom","w+")
	up=";Pilatus 2M\nphoton_energy = 12688\nadu_per_eV = 0.0001\nclen = "
	down=" \ncoffset=0 \nres = 5814.0  ; 172 micron pixel size\n\n; Define rigid group quadrant for a single panel, asic group and collections for geoptimiser\nrigid_group_q0 = 0\nrigid_group_a0 = 0\nrigid_group_collection_quadrants = q0\nrigid_group_collection_asics = a0\n\n\n; corner_{x,y} set the position of the corner of the detector (in pixels)\n; relative to the beam\n\n0/min_fs = 0\n0/max_fs = 1474\n0/min_ss = 0\n0/max_ss = 1678\n0/corner_x = -736\n0/corner_y = -858\n0/fs = x\n0/ss = y\n\nbad_beamstop/min_x = -736\nbad_beamstop/max_x = 22\nbad_beamstop/min_y = -66\nbad_beamstop/max_y = 24"
	f=open("pilatus2mpanel.geom","a+")
	f.write(up+str(coff+det_dist)+down)
	#print(up+str(coff)+down)
	f.close()
"""
This function calls cell_explorer CrystFEL script for unit cell parameters fitting
"""
def cellfit(curves,out,cell):
	count=0
	for i in curves:
		print("Save as : "+str(count)+"_"+cell)
		sub.call("cell_explorer "+out[count], shell= True)
		count+=1

"""
This function calls indexamajig CrystFEL script with a list of peak search parameters (curves), input file *.lst (inp), output files *.stream (out), geometry file *geom (geom), unit cell file *.cell (cell), string of indexing methods (idx), list of integration parameters (int), number of processors (proc)
"""
def runindexamajig(curves,inp, out, geom, cell, idx, int, proc):
        #curves rebe parâmetros de uma rodada do indexamjig
	count=1
	f=open(out+".tab","w+")
	i=0
	meu_string=" "
	for i in curves:
		if i[0]=='zaef':
			meu_string="indexamajig -i "+inp+" -g "+geom+" --peaks=zaef --indexing="+idx+" --min-peaks=4 --integration="+int[0]+"-"+int[1]+" --int-rad="+int[2]+" -o "+out+".stream -j "+str(proc)+" --threshold="+str(i[1])+" --min-squared-gradient="+str(i[2])+" --min-snr="+str(i[3])+" --peak-radius="+str(i[4])+" --profile -p "+cell
			if len(i)>5:
				if i[5]==0: meu_string=meu_string+" --filter-noise"
				if i[5]==' ': meu_string=meu_string
				else: meu_string=meu_string+" --median-filter="+str(i[5])
		if i[0]=='peakfinder8':
			meu_string="indexamajig -i "+inp+" -g "+geom+" --peaks=peakfinder8 --indexing="+idx+" --min-peaks=4 --integration="+int[0]+"-"+int[1]+" --int-rad="+int[2]+" -o "+out+".stream -j "+str(proc)+" --threshold="+str(i[1])+" --min-snr="+str(i[2])+" --min-pix-count="+str(i[3])+" --max-pix-count="+str(i[4])+" --local-bg-radius="+str(i[5])+" --profile -p "+cell
		#print(meu_string)
		peakopt.grepindexamajig(meu_string, out+".stream", count)
		sub.call("mv output.tab "+out+".tab",shell=True)
	sub.call("mkdir "+out+"/; mv "+out+".tab "+out+"/; mv "+out+".stream "+out+"/ ;mv *.tab "+out+"/", shell=True)
"""
This function calls indexamajig changing detector distance to sample from vmin to vmax in steps of step.
"""
def clen_opt(vmin,vmax,step, inp, geom, curves,cellf):
	n=int(((vmax-vmin)/step)+1)
	offset=np.linspace(vmin,vmax,n)
	np.round(offset,6)
	#print(offset)
	run=0
	cell=[cellf]
	det=0.125

	for i in offset:
		gen_geom(i, det)
		count=0
		#print(offset)
		for j in curves:
			runindexamajig([j],inp,"clen_"+str(count)+"_"+str(run), geom, cell[0], 'mosflm-latt-nocell',['rings','nocen','4,5,7'],400)
			methods,total=peakopt.filesearch_crystal("clen_"+str(count)+"_"+str(run)+"/clen_"+str(count)+"_"+str(run)+".tab",1)
			crystplots.plot_hist(methods, total, 0)
			#crystplots.plot_crystals_cell(methods, total, 0)
			sub.call("mv hist_0.png hist_"+str(run)+".png", shell=True)
			#sub.call("mv cell_id_0.png cell_id_"+str(run)+".png", shell=True)
			sub.call("mv *.png clen_"+str(count)+"_"+str(run)+"/", shell=True)
			count+=1
		run+=1
"""
This function calls detector-shift CrystFEL script 3 times in arrow, user should select beam shift clusters to correct beam center position.
"""
def detshift(inpf,curves, geomf,cellf):
	#curves=[['zaef',70,5000,10,[4,5,7],1],['peakfinder8',40,10,2,200,3],['zaef',70,5000,8,[2,3,4],0],['peakfinder8',70,10,1,200,2]]
	runs=[1,2,3]
	tests=range(len(curves))
	geom=[geomf]
	inp=[inpf]
	for i in runs:
		for j in tests:
			meu_string="pilatus2mpanel-"+str(j)+"-"+str(i)+".geom"
			geom.append(meu_string)
	print(geom)
	stream=[]
	for i in tests:
		stream.append("calibration_"+str(i)+".stream")

	for i in runs[:-1]:
		for j in tests:
			meu_string="detshift_"+str(j)+"_run_"+str(i)+"/detshift_"+str(j)+"_run_"+str(i)+".stream"
			stream.append(meu_string)
	print(stream)
	out=[]
	cell=[cellf]
	for i in tests:
		out.append("detshift_"+str(i))
	#print(cell)
	#print(out)
	#ex
	#geom=["pilatus2mpanel.geom", "pilatus2mpanel-0-1.geom", "pilatus2mpanel-1-1.geom","pilatus2mpanel-0-2.geom", "pilatus2mpanel-1-2.geom", "pilatus2mpanel-0-3.geom", "pilatus2mpanel-1-3.geom" ]
	#stream=["calibration_0.stream","calibration_1.stream","detshift_0_run_1/detshift_0_run_1.stream","detshift_1_run_1/detshift_1_run_1.stream", "detshift_0_run_2/detshift_0_run_2.stream","detshift_1_run_2/detshift_1_run_2.stream" ]
	#inp=["filescbf.lst"]
	#out=["detshift_0","detshift_1"]

	#refinement_run1
	
	sub.call("mkdir norefine/",shell=True)
	for j in tests:
		
		sub.call("./detector-shift "+stream[j]+" "+geom[0]+">> "+out[j]+".tab", shell=True)
		sub.call("mv pilatus2mpanel-predrefine.geom "+geom[j+1],shell=True)
		runindexamajig([curves[j]],inp[0],out[j]+"_run_1",geom[j+1], cell[0], 'mosflm-latt-nocell', ['rings','nocen','4,5,7'],400)
		methods,total=peakopt.filesearch_crystal(out[j]+"_run_1/"+out[j]+"_run_1.tab",1)
		crystplots.shift_map(methods, total, 0)
		sub.call("mv "+stream[j]+" norefine/",shell=True)
		sub.call("mv *.png "+out[j]+"_run_1/",shell=True)
	sub.call("mv "+geom[0]+" norefine/",shell=True)

	#refinement_run2 and run3
	cont=len(curves)
	beg=len(curves)-1

	for i in runs[1:]:
		for j in tests:
			sub.call("./detector-shift "+stream[cont]+" "+geom[cont-beg]+">> "+out[j]+".tab", shell=True)
			sub.call("mv pilatus2mpanel-"+str(j)+"-"+str(i-1)+".geom "+out[j]+"_run_"+str(i-1),shell=True)
			sub.call("mv pilatus2mpanel-"+str(j)+"-"+str(i-1)+"-predrefine.geom "+geom[cont+1],shell=True)
			runindexamajig([curves[j]],inp[0],out[j]+"_run_"+str(i),geom[cont+1], cell[0], 'mosflm-latt-nocell',['rings','nocen','4,5,7'],400)
			methods,total=peakopt.filesearch_crystal(out[j]+"_run_"+str(i)+"/"+out[j]+"_run_"+str(i)+".tab",1)
			crystplots.shift_map(methods, total, 0)
			sub.call("mv *.png "+out[j]+"_run_"+str(i)+"/",shell=True)
			cont+=1

"""
This function randonmly shuffle indexing methods order.
"""
def rand_tests(init):
	save=[]
	tests=[]
	for i in range(len(init)):
		save.append(init[i])
	tests.append(save)
	save=[]
	k=0
	iter=0
	while k<3:
		random.shuffle(init)
		for iter in range(len(tests)):
			same= [i for i, j in zip(init, tests[iter]) if i == j]
			print(same)
			if len(same)==len(tests[iter]):
				random.shuffle(init)
				iter=0
		for i in range(len(init)):
			save.append(init[i])
		tests.append(save)
		save=[]
		k+=1
	return(tests)
"""
This function perfoms all possible combination of r indexing methods from opt list of indexing methods to be tested. 
"""
def select_idx(opt, r):
	n=len(opt)
	same=[]
	tests=[]
	total=math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
	intersect=[]
	while len(tests)<total:
		flag=0
		update= random.sample(opt, r)
		for k in tests:
			same= set(k) & set(update)
			intersect.append(same)
		#print(intersect)
		for k in intersect:
			if len(k)==r:
				flag+=1
		intersect=[]
		if flag==0:
			tests.append(update)
	return(tests)

"""
This function tests different indexing methods available in CrystFEL and plot their performance comparatevely.
"""
def indexingopt(inp, out, peakpar, geom, cell,int, r, proc):
	#idx=['mosflm-nolatt-nocell','mosflm-latt-nocell','mosflm-nolatt-cell','mosflm-latt-cell','dirax-nolatt-nocell','asdf-nolatt-nocell','asdf-nolatt-cell','xds-nolatt-nocell','xds-latt-cell','xgandalf-nolatt-nocell','xgandalf-nolatt-cell','taketwo-latt-cell']
	#idx=['mosflm-nolatt-nocell','mosflm-latt-nocell','mosflm-nolatt-cell','mosflm-latt-cell','dirax-nolatt-nocell','asdf-nolatt-nocell','asdf-nolatt-cell','xds-nolatt-nocell','xds-latt-cell','xgandalf-nolatt-nocell','xgandalf-nolatt-cell']
	#at least 3 options should be passed to rand_tests it will return 4 tests with the methods passed but in random positions
	#idx=['xds-latt-cell','xgandalf-nolatt-cell','mosflm-nolatt-nocell']
	#tests=rand_tests(idx)
	#string combination of methods in 4 random orders
	count=0
	a=0
	b=0
	c=0
	d=0
	e=0
	xg=0
	
	labels=[]
	tick=""
	#tests=select_idx(idx,r)
	#tests=[['asdf-nolatt-cell'], ['mosflm-nolatt-nocell'], ['mosflm-latt-cell'], ['mosflm-latt-nocell'], ['dirax-nolatt-nocell'], ['xds-nolatt-nocell'], ['taketwo-latt-cell'], ['xgandalf-nolatt-cell'], ['xds-latt-cell'], ['xgandalf-nolatt-nocell'], ['mosflm-nolatt-cell'], ['asdf-nolatt-nocell']]
	tests=[['xds-latt-cell,mosflm-latt-cell,xgandalf-nolatt-cell,mosflm-latt-nocell,mosflm-nolatt-cell,mosflm-nolatt-nocell,dirax-nolatt-nocell,xgandalf-nolatt-nocell,asdf-nolatt-nocell,xds-nolatt-nocell']]
	#print(tests)
	f=open("xlabel.tab","w+")
	f.write(str(tests))
	f.close()

	for i in tests:
		idx=""
		for j in i:
			idx=idx+j+","
			if j[0]=='a':
				a+=1
				if j[12]=='n':
					tick+='ANLNC'
				if j[12]=='c':
					tick+='ANLC'
			if j[0]=='d':
				b+=1
				tick+='DNLNC'
			if j[0]=='m':
				c+=1
				if j[7]=='l':
					if j[12]=='n':
						tick+='MLNC'
					if j[12]=='c':
						tick+='MLC'
				if j[7]=='n':
					if j[14]=='c':
						tick+='MNLC'
					if j[14]=='n':
						tick+='MNLNC'
			if j[0]=='t':
				d+=1
				tick+='TLC'
			if j[0]=='x' and j[1]=='d':
				e+=1
				if j[4]=='n':
					tick+='XDNLNC'
				if j[4]=='l':
					tick+='XDLC'
			if j[0]=='x' and j[1]=='g':
				xg+=1
				if j[16]=='n':
					tick+='XGNLNC'
				if j[16]=='c':
					tick+='XGNLC'
			tick+="+"
		labels.append(tick[:-1])
		tick=""
		runindexamajig(peakpar,inp,out+'_'+str(r)+'_'+str(count),geom,cell,idx[:-1],int,proc)
		peakopt.fileformat(out+'_'+str(r)+'_'+str(count)+"/"+out+'_'+str(r)+'_'+str(count)+'.tab',"peakopt.tab",[[0,count]], [0])
		sub.call('mv peakopt.tab '+out+'_'+str(r)+'_'+str(count)+"/", shell=True)
		sub.call('mkdir idx_'+str(r)+'/; mv '+out+'_'+str(r)+'_'+str(count)+'/ idx_'+str(r)+'/', shell=True)
		count+=1
	print('asdf',a,'dirax',b,'mosflm',c,'take',d,'xds',e,'xgandalf',xg)
	print(labels)
	crystplots.plot_idx(labels,r)
	return(labels)

"""
This function optimizes integration parameters.
"""

def intopt(inp,out, peakpar,geom,cell):
	#int=[['rings','nocen'],['rings','cen'],['prof2d','cen'],['prof2d','nocen']]
	int=[['prof2d','cen']]
	k=0
	rad=["2,3,4","2,4,5","3,4,5","3,5,7","4,5,7"]
	for i in int:
		f=open("output.tab","w+")
		f.write("Integration optmization method: "+str(i))
		for j in rad:
			runindexamajig(peakpar,inp,out,geom,cell,'mosflm',i+[j])
			sub.call('cat '+out+'/'+out+'.tab>>summary.txt', shell=True)
		sub.call('cat summary.txt >output.tab', shell=True)
		sub.call('rm summary.txt', shell=True)
		f.close()
		methods,total=peakopt.filesearch_crystal('output.tab')
		crystplots.plot_crystals_vol(methods,total,rad)
		sub.call("mv *.png "+out+"/"+";mv output.tab output_"+str(k)+".tab", shell=True)
"""
This function runs merging process for partialator or process_hkl methods. Here, scaling, partialities and/or post-refinement can be turned on/off. 
"""

def runmerge(dir, y, adu, method):
	meu_string=""
	rx = re.compile(r'\.(stream)')
	r = []
	for path, dnames, fnames in os.walk(dir):
		r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])
	beg=len(dir)
	k=beg+1
	name=''
	out=[]
	for i in r:
		while i[k]!='/':
			name=name+i[k]
			k+=1
		out.append(name+".hkl")
		name=''
		k=beg+1

	count=0
	"""
	method=int(input("Merge with 1-process_hkl or 2-partialator:"))
	err=input("You chose to merge with "+str(method)+".Do you wanna edit your option y/n")
	while err!='n':
		method=int(input("Merge with 1-process_hkl or 2-partialator:"))
		err=input("You chose to merge with "+str(method)+".Do you wanna edit your option y/n")
	"""
	if method==2:
		#partialator
		shelf=["partialator -j 200","-i "," -o "," -y "," --max-adu="," --model="," --iterations="," --no-pr"," --no-scale"]
		sca=1
		part=1
		pr=1
		"""
		sca=int(input("Merge with scaling 1-yes 0-no:"))
		part=int(input("Merge with partialities 1-yes 0-no:"))
		pr=int(input("Merge with post-refinement 1-yes 0-no:"))
		err=(input("You chose (1) yes (0) no: scaling "+str(sca)+" partialities "+str(part)+" and post-refine "+str(pr)+". Do you wanna edit your options y/n"))
		while err!='n':
			sca=int(input("Merge with scaling 1-yes 0-no:"))
			part=int(input("Merge with partialities 1-yes 0-no:"))
			pr=int(input("Merge with post-refinement 1-yes 0-no:"))
			err=(input("You chose (1) yes (0) no: scaling "+str(sca)+" partialities "+str(part)+" and post-refine "+str(pr)+". Do you wanna edit your options y/n"))
		"""
		for i in r:
			if sca==0 and part==0 and pr==0:
                        	curves=[[" ",i, out[count], y, adu, 'unity', '0']]
			if sca==1 and part==0 and pr==0:
				curves=[[" ",i, out[count], y, adu, 'unity', '1']]
			if sca==0 and part==1 and pr==0:
				curves=[[" ",i, out[count], y, adu, 'xsphere', '0']]
			if sca==1 and part==1 and pr==0:
				curves=[[" ",i, out[count], y, adu, 'xsphere', '1', ' ']]
			if sca==1 and part==1 and pr==1:
				curves=[[" ",i, out[count][:-4]+"_0"+out[count][-4:], y,adu, 'xsphere', '1']]
			if sca==0 and part==1 and pr==1:
				curves=[[" ",i, out[count], y, adu, 'xsphere', '1', " "]]
				shelf.pop(7)
			flag=0
			for k in curves:
				for j in range(len(k)):
					meu_string=meu_string+shelf[j]+k[j]
				#print(meu_string)
				sub.call(meu_string,shell=True)
				meu_string=""
				sub.call("mv pr-logs/ pr-logs_"+out[count][:-4]+"_"+str(flag)+"/", shell=True)
				flag+=1
			count+=1
		sub.call("mkdir partialator/; mv pr-logs*/ partialator/; mv *hkl* partialator/", shell=True)
	if method==1:
                #process_hkl
		shelf=["process_hkl","-i "," -o "," -y "," --max-adu="," --scale"]
		for i in r:
			curves=[[" ", i, out[count], y, adu, " "],[" ", i, out[count]+"1 --even-only", y, adu, " "],[" ", i, out[count]+"2 --odd-only", y, adu, " "]]
			for k in curves:
				for j in range(len(shelf)):
					meu_string=meu_string+shelf[j]+k[j]
				sub.call(meu_string, shell=True)
				meu_string=""
			count+=1
		sub.call("mkdir process_hkl/; mv *hkl* process_hkl/", shell=True)

"""
This function calls compare_hkl and check_hkl for figures of merit calculation.
"""
def calcfig(dir, cell, y,method):

	meu_string=""
	count=0
	rx = re.compile(r'\.(hkl$)')
	r = []
	for path, dnames, fnames in os.walk(dir):
		r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])

	beg=len(dir)
	k=beg+1
	name=''
	out=[]

	for i in r:
		while i[k]!='.':
			name=name+i[k]
			k+=1
		out.append(dir[:4]+"_"+name)
		name=''
		k=beg+1
	print(out)

	f=open("overall.dat","w+")

	if method==1:
		#check_hkl
		shelf=["check_hkl "," -p "," --shell-file="," -y "," --ignore-negs"," --zero-negs", " --ltest", " --wilson"]
		for i in r:
			curves=[[i, cell, out[count]+".dat", y]]
			for k in curves:
				for j in range(len(k)):
					meu_string=meu_string+shelf[j]+k[j]
				resumo=sub.run(meu_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
				resumoout=resumo.stdout.decode('utf-8')
				resumoerr=resumo.stderr.decode('utf-8')
				f=open("overall.dat","a+")
				f.write(resumoerr)
				f.close()
				print(meu_string)
				meu_string=""
			count+=1
		sub.call("mkdir check_hkl/; mv *.dat check_hkl/", shell=True)

	if method==2:
	#compare_hkl
		shelf=["compare_hkl "," -p "," --shell-file="," -y "," --fom="," --zero-negs", " --ignore-negs", " -u"]
		fom=["Rsplit","CC", "CCstar"]
		for i in r:
			for k in fom:
				curves=[i+'1 '+i+'2', cell,out[count]+"_"+k+".dat", y]
				print(curves)
				for j in range(len(curves)+1):
					if j==4:
						meu_string=meu_string+shelf[j]+k
					else:
						meu_string=meu_string+shelf[j]+curves[j]
				#print(meu_string)
				resumo=sub.run(meu_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
				resumoout=resumo.stdout.decode('utf-8')
				resumoerr=resumo.stderr.decode('utf-8')
				f=open("overall.dat","a+")
				f.write(resumoerr)
				f.close()
				print(meu_string)
				meu_string=""
			count+=1
		sub.call("mkdir compare_hkl/; mv *.dat compare_hkl/", shell=True)
"""
This function calls create-mtz CrystFEL script with user's specific  space group and unit cell parameters.
"""

def convertmtz(inp):
	group='P43212'
	cell="79 79 37 90 90 90"
	f=open("create-mtz","r")
	g=open("create-mtz-mod","+w")
	g=open("create-mtz-mod","+w")
	count=0
	for x in f:
		if count==30:
			g.write("CELL "+cell+"\n")
		elif count==31:
			g.write("SYMM "+group+"\n")
		else:
			g.write(x)
		count+=1
	g.close()
	f.close()
	sub.call("chmod +x create-mtz-mod;./create-mtz-mod "+inp, shell=True)

"""
This function calls create-xscale CrystFEL script with user's specific space group and unit cell parameters. 
"""

def convertxscale(inp):
        group='96'
        cell="79 79 38 90 90 90"
        f=open("create-xscale","r")
        g=open("create-xscale-mod","+w")
        g=open("create-xscale-mod","+w")
        count=0
        for x in f:
                if count==14:
                        g.write(x[0:31]+cell+x[-6:])
                elif count==13:
                        g.write(x[0:28]+group+x[-6:])
                else:
                        g.write(x)
                count+=1
        g.close()
        f.close()
        sub.call("chmod +x create-xscale-mod;./create-xscale-mod "+inp+">"+inp[:-4]+"_xds.hkl", shell=True)

def main():
	#### Start detector corrections #######

	#curves=[['zaef',70,5000,5,[4,5,7],1],['peakfinder8',40,5,2,200,3],['zaef',200,25000,5,[3,4,5]],['peakfinder8',200,6,2,200,3]]
	#curves=[['zaef',80,5000,4,'2,3,4',0],['zaef',80,5000,4,'2,3,4',6],['zaef',100,5000,4,'3,4,5',6],['peakfinder8',80,5,2,50,6],['peakfinder8',100,5,2,50,6]]
	#curves=[['zaef',10,5000,4,'4,5,7',' '],['zaef',100,5000,4,'4,5,7',' '],['peakfinder8',10,3,2,20,6],['peakfinder8',100,3,2,20,6]]
	curves=[['zaef',10,5000,4,'4,5,7',' ']]
	stream=["calibration_0.stream"]
	inp='files_500.lst'
	geom="pilatus2mpanel.geom"
	
	"""
	Detector distance to sample optimization.
	"""

	#cell="calib.cell"
	#clen_opt(-5*1e-4,5*1e-4,1e-4, inp, geom, curves,cell)
	#clen_opt(4*1e-4,6*1e-4,2*1e-5, inp, geom, curves,cell)
	det_dist=0.1255
	#update clen in the geom file call finalpeakopt ame.cell to see if the indexing is better

	"""
	Unit cell file, if no unit cell file should be included cell=0
	"""
	#cell=0
	#cell="liso.cell"

	#peakopt.finalpeakopt(inp,curves,geom,cell,400)
	#peakopt.fileformat("liso_1_0.tab","peakopt.tab",curves, [0])
	#methods, total=peakopt.filesearch_crystal('liso_1_0.tab',1)

	"""
	Plot of unit cell parameters histograms.
	"""

	'''
	for i in range(len(curves)):
		print(i)
		crystplots.plot_hist(methods,total, i)
		crystplots.plot_crystals_cell(methods,total, i)
	'''

	#crystplots.compare_hist(methods,total)
	
	"""
	Unit cell fitting with cell_explorer from CrystFEL
	"""

	#cellfit(curves, stream, cell)


	"""
	Beam shift position corrections
	"""
	#detshift(inp,curves, geom, cell)

	#find det_0_run_3.tab to see if scores are better than before det_shift

	####### End detector corrections ########

	

	############## Start indexing optimization ###################
	#Choose one peakparam option from now
	#Indexing long file.lst (few thousands of images) and unit cell parameters file fitted in last step.

	inp='files_16.lst'
	geom="pilatus2mpanel.geom"
	cell="0_liso_after.cell"
	out='intopt'
	"""
	Integration parameters optimization
	"""
	#intopt(inp,out, curves,geom,cell)

	"""
	Indexing methods optimization
	"""

	out="liso"
	int=['rings','nocen','4,5,7']
	r=1

	#labels=indexingopt(inp,out,curves,geom,cell,int,r, proc)

	xlabel=['0','1','2','3','4','5','6','7','8','9','10','11']
	#crystplots.plot_idx(xlabel,1)
	

	###############End indexing optimization##################
	
	############### Start merging and FoM calculation ##########################

	results='./results/'
	stream=[results+'idx_all_64/idx_1/',results+'idx_all_48/idx_1/',results+'idx_all_32/idx_1/',results+'idx_all_16/idx_1/']
	methods=[1,2]
	dir_compare=[]
	dir_check=[]

	"""
	Merging and FoM calc
	"""

	'''
	for j in methods:
		for i in stream:
			#runmerge(i,'4/mmm','100000',j)
			#sub.call('mv p*/ '+i, shell=True)
			calcfig(i+'partialator',cell,'4/mmm',j)
			sub.call('mv c*/ '+i+'partialator', shell=True)
			calcfig(i+'process_hkl',cell,'4/mmm',j)
			sub.call('mv c*/ '+i+'process_hkl', shell=True)
			sub.call('mv p*/ '+i, shell=True)
	'''
	
	"""
	FoM comparative plots
	"""

	for i in stream:
		dir_check.append(i+'partialator/check_hkl/')
		dir_compare.append(i+'partialator/compare_hkl/')

	label=[['64'],['48'],['32'],['16']]
	crystplots.plot_check(label,'part','0', dir_check)
	crystplots.plot_compare(label,'part','0', dir_compare)
        

	#ame y=2_uab
	#liso y=4/mmm

	################### End merging and FoM calculation #########################

	################## Start conversions ########################################

	#user decide final data to convert
	#inp="process_hkl_1/_3.hkl"
	#convertmtz(inp)
	#convertxscale(inp)

	################# End conversions ###########################################
