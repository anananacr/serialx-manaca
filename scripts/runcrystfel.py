import re
import math
import subprocess as sub
import numpy as np
import os
import sys
import argparse
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
from proc_config import proc_param

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
def cellfit(curves,out):
	count=0
	for i in curves:
		print("Save as : "+str(count)+"_cell.cell")
		sub.call("cell_explorer "+out[count], shell= True)
		count+=1

"""
This function calls indexamajig CrystFEL script with a list of peak search parameters (curves), input file *.lst (inp), output files *.stream (out), geometry file *geom (geom), unit cell file *.cell (cell), string of indexing methods (idx), list of integration parameters (int), number of processors (proc)
"""
def runindexamajig(curves,inp, out, geom, cell, idx, int, proc):
	#set up your environment in cmd	
	cmd='source /etc/profile.d/modules.sh; module load xray; module load hdf5-openmpi/1.10.5; module load crystfel; module load hdf5/1.10.5'
	sub.call(cmd,shell=True)
	#curves rebe parâmetros de uma rodada do indexamjig
	count=1
	f=open(out+".tab","w+")
	i=0
	meu_string=" "
	for i in curves:
		if i[0]=='zaef':
			if idx==0:
				meu_string="/software/crystfel/0.10.1/bin/indexamajig -i {inp} -g {geom} --peaks=zaef --min-peaks=14 --integration={int[0]}-{int[1]} --int-rad={int[2]} -o {out}.stream -j {proc} --threshold={i[1]} --min-squared-gradient={i[2]} --min-snr={i[3]} --peak-radius={i[4]} --profile -p {cell}"
			else:
				meu_string="/software/crystfel/0.10.1/bin/indexamajig -i {inp} -g {geom} --peaks=zaef --indexing={idx} --min-peaks=14 --integration={int[0]}-{int[1]} --int-rad={int[2]} -o {out}.stream -j {proc} --threshold={i[1]} --min-squared-gradient={i[2]} --min-snr={i[3]} --peak-radius={i[4]} --profile -p {cell}"
				if len(i)>5:
					if i[5]==0: meu_string=meu_string+" --filter-noise"
					if i[5]==' ': meu_string=meu_string
					else: meu_string=meu_string+" --median-filter="+str(i[5])
		if i[0]=='peakfinder8':
			if idx==0:
				meu_string=f"/software/crystfel/0.10.1/bin/indexamajig -i {inp} -g {geom} -o {out}.stream -j {proc} -p {cell} --peaks=peakfinder8 --threshold={i[1]} --min-snr={i[2]} --min-pix-count={i[3]} --max-pix-count={i[4]} --local-bg-radius={i[5]} --min-res=10 --max-res=1000 --min-peaks=14 --peak-radius=4.0,5.0,7.0 --integration={int[0]}-{int[1]} --int-rad={int[2]} --multi --profile "
			else:
				meu_string=f"/software/crystfel/0.10.1/bin/indexamajig -i {inp} -g {geom} -o {out}.stream -j {proc} -p {cell} --peaks=peakfinder8 --indexing={idx} --threshold={i[1]} --min-snr={i[2]} --min-pix-count={i[3]} --max-pix-count={i[4]} --local-bg-radius={i[5]} --min-res=10 --max-res=1000 --min-peaks=14 --peak-radius=4.0,5.0,7.0 --integration={int[0]}-{int[1]} --int-rad={int[2]} --multi --profile "
		print(meu_string)
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
		#runindexamajig(peakpar,inp,out+'_'+str(r)+'_'+str(count),geom,cell,idx[:-1],int,proc)
		runindexamajig(peakpar,inp,out+'_'+str(r)+'_'+str(count),geom,cell,0,int,proc)
		peakopt.fileformat(out+'_'+str(r)+'_'+str(count)+"/"+out+'_'+str(r)+'_'+str(count)+'.tab',"peakopt.tab",[[0,count]], [0])
		sub.call('mv peakopt.tab '+out+'_'+str(r)+'_'+str(count)+"/", shell=True)
		sub.call('mkdir idx_'+str(r)+'/; mv '+out+'_'+str(r)+'_'+str(count)+'/ idx_'+str(r)+'/', shell=True)
		count+=1
	print('asdf',a,'dirax',b,'mosflm',c,'take',d,'xds',e,'xgandalf',xg)
	print(labels)
	#crystplots.plot_idx(labels,r)
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

def runmerge(dir, y, adu,label,method):
	meu_string=""
	rx = re.compile(r'\.(stream)')
	r = []
	for path, dnames, fnames in os.walk(dir):
		r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])
	beg=len(dir)
	k=beg
	name=''
	print(r)
	out=[]
	count=0
	for i in r:
		print(name)
		out.append(f"{label}_{count}.hkl")
		count+=1

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
		shelf=["/software/crystfel/0.10.1/bin/partialator -j 40 --polarisation=horiz --min-measurements=2 --min-res=inf --push-res=inf --no-logs","-i "," -o "," -y "," --max-adu="," --model="," --iterations="," --no-pr"]

		for i in r:
			curves=[[" ",i, out[count], y, adu, 'unity', '3', ' ']]
			
			flag=0
			for k in curves:
				for j in range(len(k)):
					meu_string=meu_string+shelf[j]+k[j]
				print(meu_string)
				sub.call(meu_string,shell=True)
				meu_string=""
				sub.call("mv pr-logs/ pr-logs_"+out[count][:-4]+"_"+str(flag)+"/", shell=True)
				flag+=1
			count+=1
		sub.call("mkdir partialator/; mv pr-logs*/ partialator/; mv *hkl* partialator/", shell=True)
	if method==1:
		#process_hkl
		shelf=["/software/crystfel/0.10.1/bin/process_hkl","-i "," -o "," -y "," --max-adu="," --scale --polarisation=horiz --min-measurements=2 --min-res=inf --push-res=inf"]
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
		out.append(name)
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
				print(meu_string)
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

def export_file(inp,cell,out,file_format):
	cmd=f"get_hkl -i {inp} -p {cell} -o {out} --output-format={file_format}"
	sub.call(cmd, shell=True)

def index_no_cell(var):
	peakopt.finalpeakopt(var.inp,var.curves,var.label,var.geom,var.cell,var.n_proc)
	peakopt.fileformat("output.tab","peakopt.tab",var.curves, [0])

	#methods, total=peakopt.filesearch_crystal('output.tab',0)
	
	#for i in range(len(var.curves)):
	#	print(i)
	#	crystplots.plot_hist(methods,total, i)
	#	crystplots.plot_crystals_cell(methods,total, i)
	#crystplots.compare_hist(methods,total)
	cellfit(var.curves, var.stream)
def index_all(var):
	labels=indexingopt(var.inp,var.label,var.curves,var.geom,var.cell,var.integ,var.r, var.n_proc)

def merge_all(var):

	cmd=f'mkdir ../results;rm -r ../results/{var.label};mkdir ../results/{var.label};rm -r mv ./idx_1/p*/; mv ./idx_1/ ../results/{var.label}/'
	sub.call(cmd,shell=True)
	stream=[f'../results/{var.label}/idx_1/']
	methods=[1,2]


	"""
	Merging and FoM calc
	"""

	for j in methods:
		for i in stream:
			runmerge(i,var.sym,var.sat,var.label,j)
			sub.call('mv p*/ '+i, shell=True)
	for j in methods:
		for i in stream:
			calcfig(i+'partialator',var.cell,var.sym,j)
			sub.call('mv c*/ '+i+'partialator/', shell=True)
	for j in methods:
		for i in stream:
			calcfig(i+'process_hkl',var.cell,var.sym,j)
			sub.call('mv c*/ '+i+'process_hkl', shell=True)
			sub.call('mv pr-*/ '+i, shell=True)

def show_plots(var):
	"""
	FoM comparative plots
	"""
	results='../results/'
	stream=[f'{results}/{var.label}/idx_1/']
	dir_check=[]
	dir_compare=[]
	for i in stream:
		dir_check.append(i+f'{var.merge_method}/check_hkl/')
		dir_compare.append(i+f'{var.merge_method}/compare_hkl/')
	label=[['final']]
	crystplots.plot_check(label,f'{var.merge_method}'[:4],'0', dir_check)
	crystplots.plot_compare(label,f'{var.merge_method}'[:4],'0', dir_compare)
	
	
def main(raw_args=None):
	#### Straightforward #####
	
	parser = argparse.ArgumentParser(
	description="Decode and sum up Jungfrau images.")
	parser.add_argument("-m", "--mode", type=str, action="store",
	help="mode option in CrystFEL: index_no_cell, index_cell, index_all, merge_all, fom_plots, export_mtz ")
	args = parser.parse_args(raw_args)

	var=proc_param()

	if args.mode=='index_no_cell':
		index_no_cell(var)
	if args.mode=='index_cell':
		var.cell='../0_cell.cell'
		index_no_cell(var)
	if args.mode=='index_all':
		var.cell='../0_cell.cell'
		index_all(var)
	if args.mode=='merge_all':
		var.cell='../0_cell.cell'
		merge_all(var)
	if args.mode=='fom_plots':
		show_plots(var)
	if args.mode=='export_mtz':
		var.cell='../0_cell.cell'
		inp=f"../results/{var.label}/idx_1/{var.merge_method}/{var.label}_0.hkl"
		out=f"../results/{var.label}/idx_1/{var.merge_method}/{var.label}.mtz"
		export_file(inp,var.cell,out,'mtz')
		print(f'Final .mtz file: {out}')

	
        #### Advanced adjustments #####

        #### Start detector corrections #######


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
	cell="0_cell.cell"

	#peakopt.finalpeakopt(inp,curves,geom,cell,40)
	#peakopt.fileformat("output.tab","peakopt.tab",curves, [0])
	#methods, total=peakopt.filesearch_crystal('output.tab',0)


	"""
	Plot of unit cell parameters histograms.
	"""
	
	#for i in range(len(curves)):
	#	print(i)
	#	crystplots.plot_hist(methods,total, i)
	#	crystplots.plot_crystals_cell(methods,total, i)


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

	cell="1_cell.cell"
	out='intopt'
	"""
	Integration parameters optimization
	"""
	#intopt(inp,out, curves,geom,cell)

	"""
	Indexing methods optimization
	"""

	out="tutorial"
	integ=['rings','nocen','4,5,7']
	r=1

	#labels=indexingopt(inp,out,curves,geom,cell,integ,r, 40)

	xlabel=['0','1','2','3','4','5','6','7','8','9','10','11']
	#crystplots.plot_idx(xlabel,1)
	

	###############End indexing optimization##################
	

if __name__ == '__main__':
    main()
