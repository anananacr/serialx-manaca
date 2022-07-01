import subprocess as sub
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt
import csv
from datetime import datetime
from itertools import product
import crystplots
import seaborn as sns
import pandas as pd

'''
This class saves all information obtained from each crystal founded by CrstFEL in a run of indexamajig.
'''
class crystal:
	def __init__(self, crystal_param):
		self.id=crystal_param[0]
		self.a=crystal_param[1]*10
		self.b=crystal_param[2]*10
		self.c=crystal_param[3]*10
		self.alf=crystal_param[4]
		self.bet=crystal_param[5]
		self.gam=crystal_param[6]
		self.latt=crystal_param[7]
		self.cent=crystal_param[8]
		self.detshift=crystal_param[9]
		self.idx=crystal_param[10]

'''
This function prints date in a file, which is suitable for processing time tracking
'''
def print_date(file):
	now = datetime.now()
	f=open(file,"a+") 
	f.write("\nClock "+str(now.year)+"."+str(now.month)+"."+str(now.day)+" "+str(now.hour)+":"+str(now.minute)+":"+str(now.second)+"\n")
	f.close()

'''
This function combines two lists one by one and return a list with combined components in a string separated by an space.
'''
def combine_param(curves, param):
	combined=[]
	for i in curves:
		for j in param:
			string=""
			string=str(i)+" "+str(j)
			combined.append(string)
	return combined
'''
This function runs indexamajig with my_string, save its stdout and stderr, which will give indexing rate and processing time.
Also, it searchs in the stream file all main parameters of crystals founded and write evrything in output.tab
'''
def grepindexamajig(my_string, out, count):
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f=open("output.tab","a+")
	print("\nIndexamajig finished command "+str(count)+"\n")
	f=open("output.tab","a+")
	f.write(my_string+"\n")
	f.write(resumoerr)
	f.write("\nEnd command "+str(count)+"\n")
	f.write(resumoout)
	print_date("output.tab")
	f.write("\nImages chunk\n")
	resumo=sub.run("grep -A11 'Begin chunk' "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.write("\nPeak-intensity\n")
	resumo=sub.run("./peak-intensity "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.write("\nCell parameters\n")
	resumo=sub.run("grep 'Cell parameters' "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.write("\nlattice_type\n")
	resumo=sub.run("grep 'lattice_type' "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.write("\ncentering\n")
	resumo=sub.run("grep 'centering' "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.write("\ndet_shift\n")
	resumo=sub.run("grep 'det_shift' "+out , shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	resumoerr=resumo.stderr.decode('utf-8')
	f.write(resumoerr)
	f.write(resumoout)
	f.close()

'''
This function searchs output.tab from grepindexamjig function and return indexing rate, number of crystal founded, processing time, MAP, MPP, TAP
'''
def filesearch(inp):
	my_string=" grep -B2 'End' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	index=[]
	crystals=[]
	time=[]
	for i in range(len(resumoout)):
		if resumoout[i]=='t'and resumoout[i+1]=='s' and resumoout[i+2]==',':
			number=""
			k=i+4
			while resumoout[k]!="%":
				number=number+resumoout[k]
				k+=1
			index.append(number)

	for i in range(len(resumoout)):
		if resumoout[i]=='l'and resumoout[i+1]=='l':
			string=''
			j=i+5
			while resumoout[j]!=" ":
				string=string+resumoout[j]
				j+=1
			crystals.append(string)
	my_string="grep -B12 'Clock' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	for i in range(len(resumoout)):
		if resumoout[i]==')'and resumoout[i+1]==':':
			string=''			
			j=i+11
			while resumoout[j]!=" ":
				string=string+resumoout[j]
				j+=1
			time.append(string)
	my_string="grep -A4 'Peak-intensity' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoout=resumo.stdout.decode('utf-8')
	mpp=[]
	map=[]
	tap=[]
	for i in range(len(resumoout)):
		if resumoout[i]=='a'and resumoout[i+1]=='n':
			string=''
			j=i+3
			while resumoout[j]!=" ":
				string=string+resumoout[j]
				j+=1
			if resumoout[j+1]=="p":mpp.append(string)
			if resumoout[j+5]=="p":map.append(string)
			if resumoout[j+5]=="t": tap.append(string)
	return(index, crystals, time,mpp,map,tap)

'''
This function searchs output.tab generated in grepindexamajig function, and creates a list of crystals founded with their main parameters associated, for every method tested in the peak searching optimization
'''

def filesearch_crystal(inp, cell_flag):
	my_string="grep 'Image serial number' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoid=resumo.stdout.decode('utf-8')
	my_string="grep 'Cell parameters' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumocell=resumo.stdout.decode('utf-8')
	my_string="grep 'lattice_type' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumolatt=resumo.stdout.decode('utf-8')
	my_string="grep 'centering' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumocent=resumo.stdout.decode('utf-8')
	my_string="grep 'det_shift' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoshift=resumo.stdout.decode('utf-8')
	my_string="grep 'indexed_by' "+inp
	resumo=sub.run(my_string, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)
	resumoidx=resumo.stdout.decode('utf-8')
	id=[]
	latt=[]
	cent=[]
	cell=[]
	shift=[]
	idx=[]
	for i in range(len(resumoid)):
		if resumoid[i]==':' and resumoid[i+1]==' ':
			number=""
			k=i+2
			while resumoid[k]!="\n":
				number=number+resumoid[k]
				k+=1
			id.append(int(number))
	for i in range(len(resumolatt)):
		if resumolatt[i]=='=' and resumolatt[i+1]==' ':
			number=""
			k=i+2
			while resumolatt[k]!="\n":
				number=number+resumolatt[k]
				k+=1
			latt.append(number)
	for i in range(len(resumocent)):
		if resumocent[i]=='=' and resumocent[i+1]==' ':
			number=""
			k=i+2
			while resumocent[k]!="\n":
				number=number+resumocent[k]
				k+=1
			cent.append(number)
	for i in range(len(resumoidx)):
                if resumoidx[i]=='y' and resumoidx[i+2]=='=':
                        number=""
                        k=i+4
                        while resumoidx[k]!="\n":
                                number=number+resumoidx[k]
                                k+=1
                        idx.append(number)
	for i in range(2,len(resumoshift)):
		if resumoshift[i]=='=' and resumoshift[i-2]=='x':
			number=""
			numbery=""
			k=i+2
			while resumoshift[k]!=" ":
				number=number+resumoshift[k]
				k+=1
			k+=5
			while resumoshift[k]!=" ":
				numbery=numbery+resumoshift[k]
				k+=1
			shift.append([float(number),float(numbery)])
			number=""
			numbery=""

	for i in range(len(resumocell)):
		if resumocell[i]=='s' and resumocell[i+1]==' ':
			number=""
			k=i+2
			while resumocell[k]!="\n":
				number=number+resumocell[k]
				k+=1
			cell.append(number)
	a=[]
	b=[]
	c=[]
	alf=[]
	bet=[]
	gam=[]
	for i in cell:
		number=""
		k=0
		a_flag=1
		b_flag=0
		c_flag=0
		alf_flag=0
		bet_flag=0
		gam_flag=0
		while i[k]!=" ":
			number=number+i[k]
			k+=1
			if i[k+1]==" " and a_flag==1:
				k+=2
				a.append(float(number))
				a_flag=0
				b_flag=1
				number=""
			if i[k+1]==" " and b_flag==1:
				k+=2
				b.append(float(number))
				b_flag=0
				c_flag=1
				number=""
			if i[k+1]==" " and c_flag==1:
				k+=6
				c.append(float(number))
				number=""
				c_flag=0
				alf_flag=1
			if i[k+1]==" " and alf_flag==1:
				k+=2
				alf.append(float(number))
				number=""
				alf_flag=0
				bet_flag=1
			if i[k+1]==" " and bet_flag==1:
				k+=2
				bet.append(float(number))
				number=""
				bet_flag=0
				gam_flag=1
			if i[k+1]==" " and gam_flag==1:
				gam.append(float(number))
				number=""
				gam_flag=0
	results=filesearch(inp)
	total=results[1]
	list_crystal=[]
	methods=[]
	k=0
	j=0
	off=0
	if cell_flag==1:
		latt.pop(0)
		cent.pop(0)
	print(len(idx), len(id), len(a),len(b), len(c), len(alf), len(cent), len(latt), len(shift))
	print(int(total[k]))
	for i in range(len(id)):
		if idx[i]=='none':
			i+=1
		else:
			x=crystal([id[i],a[j],b[j],c[j],alf[j],bet[j],gam[j],latt[j],cent[j],shift[j],idx[i]])
			list_crystal.append(x)
			if j==(int(total[k])-1+off):
				methods.append(list_crystal)
				list_crystal=[]
				off=off+int(total[k])
				k+=1
			j+=1
	return(methods, total)

'''
This function writes filesearch function output in a table and save in out.
'''

def fileformat(inp,out, curves, param):
	test=[]
	tmp=[]
	print(inp)
	for i in range(len(curves)):
		if len(curves[0])==2:
			test.append(curves[i][1])
		else:
			test.append(i+1)
	if len(curves[0])==4 and curves[0][0]=="zaef":
		for i in range(len(param)):
			tmp.append(i+1)
		param=tmp
	if len(curves[0])==6 and curves[0][0]=="peakfinder8":
		for i in range(len(param)):
			tmp.append(i+1)
		param=tmp
	combined=combine_param(test, param)
	print(combined)
	results=filesearch(inp)
	print(results)
	index=results[0]
	crystals=results[1]
	time=results[2]
	mpp=results[3]
	map=results[4]
	tap=results[5]
	f=open(out,"w+")
	f=open(out,"a+")
	if time==[]:
		f.write("thr sqrd index crystals mpp map tap\n")
		for i in range(len(combined)):
			f.write(combined[i]+" "+index[i]+" "+crystals[i]+" "+mpp[i]+" "+map[i]+" "+tap[i]+"\n")
	else:
		f.write("thr sqrd index crystals time mpp map tap\n")
		for i in range(len(combined)):
			f.write(combined[i]+" "+index[i]+" "+crystals[i]+" "+time[i]+" "+mpp[i]+" "+map[i]+" "+tap[i]+"\n")
	f.close()

'''
This function combines command line options of indexamajig with curves you already have tested and parameters you want to test in next peak seach optimization. 
It builds my_string that will be used to call indexamajig in the grepindexamaijg function.
'''

def optloop(inp,out, curves,param, geom):
	count=1
	f=open("output.tab","w+")
	#["zaef",100,5000,2]
	step=len(curves[0])
	#print(step)
	if curves[0][0]=="zaef":
		shelf=[" --peaks=", " --threshold="," --min-squared-gradient="," --min-snr=", " --peak-radius=", " --median-filter="]
		if step==6 and param[0]==" ":
			shelf=[" --peaks=", " --threshold="," --min-squared-gradient="," --min-snr=", " --peak-radius="," --filter-noise", " "]
			#print("only filter")
		if step==6 and param[0]!=" ":
			shelf=[" --peaks=", " --threshold="," --min-squared-gradient="," --min-snr=", " --peak-radius="," --filter-noise", " --median-filter="]
			#print("MF+NOI")
	if curves[0][0]=="peakfinder8":
		shelf=[" --peaks=", " --threshold="," --min-snr="," --min-pix-count="," --max-pix-count="," --local-bg-radius="," --min-res="," --max-res="]
	for i in curves:
		k=0
		my_string="indexamajig -i "+inp+" -g "+geom+" --indexing=mosflm -o "+out+" -j 28 --profile --no-refls-in-stream"
		while k<step:
			my_string=my_string+shelf[k]+str(i[k])
			k+=1
			save_string=my_string
		for j in param:
			if step==6 and curves[0][0]=='peakfinder8':
				my_string=my_string+shelf[step]+str(j[0])+shelf[step+1]+str(j[1])
			else:
				my_string=my_string+shelf[step]+str(j)
			grepindexamajig(my_string, out,count)
			count+=1
			my_string=save_string

	#sub.call("mkdir peakopt_"+str(curves[0][0])+"_"+str(step)+"/; mv output.tab peakopt_"+str(curves[0][0])+"_"+str(step)+"/;", shell=True)


'''
This function saves stream sfiles for all final curves you have tested enabling their perfomance comparison.
'''

def finalpeakopt(inp,curves,geom,cell,proc):
	count=1
	f=open("output.tab","w+")
	i=0

	while i<(len(curves)):
		#print(i)
		out="calibration_"+str(i)+".stream"
		if curves[i][0]=='zaef':
			my_string="indexamajig -i "+inp+" -g "+geom+" --indexing=mosflm-latt-nocell --min-peaks=4 --peaks=zaef -o "+out+" -j "+str(proc)+" --threshold="+str(curves[i][1])+" --min-squared-gradient="+str(curves[i][2])+" --min-snr="+str(curves[i][3])+" --peak-radius="+str(curves[i][4])+" --profile"
			if curves[i][5]==0: my_string=my_string+" --filter-noise"
			elif curves[i][5]==' ': my_string=my_string
			else: my_string=my_string+" --median-filter="+str(curves[i][5])
		if curves[i][0]=='peakfinder8':
			my_string="indexamajig -i "+inp+" -g "+geom+" --indexing=mosflm-latt-nocell --min-peaks=4 --peaks=peakfinder8 -o "+out+" -j "+str(proc)+" --threshold="+str(curves[i][1])+" --min-snr="+str(curves[i][2])+" --min-pix-count="+str(curves[i][3])+" --max-pix-count="+str(curves[i][4])+" --local-bg-radius="+str(curves[i][5])+" --profile"
		if cell!=0:
			my_string=my_string+" -p "+cell
		grepindexamajig(my_string,out, count)
		count+=1
		i+=1
'''
This function is an attempt to communicate with users asking which peak search parameters in zaef method they want to test. In practice is better to pass parameters directly in the main function.
'''

def zaef_peakopt(inp,out,geom):
	#zaef method
	method='zaef'
	list_param=set_param(method)

	print("Peak search optmization: "+method+" method")
	count=int(input("From which step do you wanna begin? 1, 2, 3 or 4"))-1
	erro=input("Peak search optmization from step "+str(count+1)+". Do you wanna retype? y/n ")
	while erro!='n':
		count=int(input("From which step do you wanna begin? 1, 2, 3 or 4"))-1
		erro=input("Peak search optmization from step "+str(count+1)+". Do you wanna retype? y/n ")
	#count=0 start from the beggining
	ele=[]
	list_curves=[]
	#print(list_param[0:-2])
	for i in list_param[count:-2]:
		curves = []
		n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
		erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
		while erro!='n':
			 n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
			 erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")

		for j in range(0, n):
			print("Enter your peaksearch parameters option number "+str(j+1)+":")
			for k in range(0, count+2):
				if k==0: ele.append(method)
				else: ele.append(input())
			curves.append(ele)
			ele=[]

		print("Check your options:")
		print(curves)
		erro=input("Do you wanna retype your parameters for this step? y/n ")
		while erro!='n':
			curves = []
			n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
			erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
			while erro!='n':
				n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
				erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
			for j in range(0, n):
				print("Enter your peaksearch parameters option number "+str(j+1)+":")
				for k in range(0, count+2):
					if k==0: ele.append(method)
					else: ele.append(input())
				curves.append(ele)
				ele=[]
			print("Check your options:")
			print(curves)
			erro=input("Do you wanna retype your parameters for this step? y/n ")

		list_curves.append(curves)

		optloop(inp,out,curves,i,geom)
		fileformat("output.tab", "peakopt.tab", curves, i)
		
		if count!=3:
			crystplots.plot("zaef", "mpp",curves,i)
			crystplots.plot_thr("zaef","mpp",curves,i)
			sub.call("mv peakopt_zaef_"+str(count+1)+"/peakopt.tab .", shell=True)
			crystplots.plot("zaef", "map",curves,i)
			crystplots.plot_thr("zaef","map",curves,i)
			sub.call("mv peakopt_zaef_"+str(count+1)+"/peakopt.tab .", shell=True)
			crystplots.plot("zaef", "index",curves,i)
			crystplots.plot_thr("zaef","index",curves,i)

		count+=1
	sub.call("mkdir peakopt_zaef_4; mv *t.tab peakopt_zaef_4", shell=True)
	save=list_curves[-1]
	tmp2=[]
	curves_tmp=[]
	for i in curves:
		tmp2=i
		tmp2.append(" ")
		curves_tmp.append(tmp2)
		tmp2=[]
	count=0
	for i in save:
		curves[count]=i[0:-1]
		count+=1

	print("Start MF+NOI tests")
	list_curves[-1]=curves
	list_curves.append(curves_tmp)
	
	#curves=[['zaef', '80', '5000', '4', '2,3,4', ' '], ['zaef', '80', '5000', '4', '3,5,7', ' '], ['zaef', '100', '5000', '4', '3,4,5', ' '], ['zaef', '100', '5000', '4', '2,3,4', ' '], ['zaef', '200', '25000', '5', '2,3,4',' ']]
	optloop(inp,out,curves_tmp,list_param[-2],geom)
	fileformat("output.tab", "peakopt.tab", curves_tmp, list_param[-2])
	sub.call("mkdir peakopt_zaef_5; mv *t.tab peakopt_zaef_5", shell=True)
	list_curves.append(curves_tmp)
	print("Start NOI tests")
	#print(len(curves_tmp),curves_tmp,list_param[-1], list_param[-1][0])
	optloop(inp,out,curves_tmp,list_param[-1],geom)
	param=[0]
	fileformat("output.tab", "peakopt.tab", curves_tmp, param)
	sub.call("mkdir peakopt_zaef_6; mv *t.tab peakopt_zaef_6", shell=True)
	
	curves=list_curves[-3]
	param=list_param[-2]
	crystplots.plot_median("mpp",curves,param)
	crystplots.plot_median_thr("mpp",curves,param)
	crystplots.plot_median("map",curves,param)
	crystplots.plot_median_thr("map",curves,param)
	
	return(list_curves)


'''
This function is an attempt to communicate with users asking which peak search parameters in peakfinder8 method they want to test. In practice is better to pass parameters directly in the main function.
'''
def peakfinder8_peakopt(inp,out,geom):
	#peakfinder8 method
	method='peakfinder8'
	list_param=set_param(method)
	print("Peak search optmization: "+method+" method")
	count=int(input("From which step do you wanna begin? 1, 2, 3, 4 or 5: "))-1
	erro=input("Peak search optmization from step "+str(count+1)+". Do you wanna retype? y/n ")
	while erro!='n':
		count=int(input("From which step do you wanna begin? 1, 2, 3 or 4"))-1
		erro=input("Peak search optmization from step "+str(count+1)+". Do you wanna retype? y/n ")
	#count=0 to start from the beggining
	ele=[]
	list_curves=[]
	for i in list_param[count:]:
		curves = []
		n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
		erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
		while erro!='n':
			 n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
			 erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")

		for j in range(0, n):
			print("Enter your peaksearch parameters option number "+str(j+1)+":")
			for k in range(0, count+2):
				if k==0: ele.append(method)
				else: ele.append(input())
			curves.append(ele)
			ele=[]

		print("Check your options:")
		print(curves)
		erro=input("Do you wanna retype your parameters for this step? y/n")
		while erro!='n':
			curves = []
			n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
			erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
			while erro!='n':
				n = int(input("Enter the number of tests for peaksearch step number "+str(count+1)+"- (You should have "+str(count+1)+" parameters for each option): "))
				erro=input("You have "+str(n)+" options for this step. Do you wanna retype? y/n ")
			for j in range(0, n):
				print("Enter your peaksearch parameters option number "+str(j+1)+":")
				for k in range(0, count+2):
					if k==0: ele.append(method)
					else: ele.append(input())
				curves.append(ele)
				ele=[]
			print("Check your options:")
			print(curves)
			erro=input("Do you wanna retype your parameters for this step? y/n ")

		list_curves.append(curves)
		optloop(inp,out,curves,i,geom)
		fileformat("output.tab", "peakopt.tab", curves, i)
		crystplots.plot("peakfinder8","mpp",curves,i)
		crystplots.plot_thr("peakfinder8","mpp", curves,i)
		sub.call("mv peakopt_peakfinder8_"+str(count+1)+"/peakopt.tab .")
		crystplots.plot("peakfinder8","map",curves,i)
		crystplots.plot_thr("peakfinder8","map", curves,i)
		sub.call("mv peakopt_peakfinder8_"+str(count+1)+"/peakopt.tab .")
		crystplots.plot("peakfinder8","index",curves,i)
		crystplots.plot_thr("peakfinder8","index", curves,i)
		count+=1

	return(list_curves)

'''
This function reads list of parameters changed from default by users. In practice is better to pass parameters directly in the main function.
'''
def read_list_param(method):
	if method=="zaef": method=0
	else: method=1
	list_param=[]
	f=open("peakopt_param.tab","r")
	string=""
	param=[]
	shell=[]
	line=f.readline()
	line=f.readline()
	line=f.readline(12)
	while line=="Peakopt_step":
		line=f.readline()
		line=f.readline()
		#print(line)
		count=1
		if len(list_param)!=2 and method==0 or len(list_param)!=4 and method==1:
			while line[count]!="," and line[count]!="]":
				string=string[0:]+line[count]
				if line[count+1]==",":
					if string[0]!=" ":
						param.append(string)
					else: param.append(string[1:])
					string=""
					count+=2
				else:
					count+=1
			param.append(string[1:])
			#print(string)
			string=""
			list_param.append(param)
			param=[]
			line=f.readline(12)
		if len(list_param)==2 and method==0 and count==1:
			while (line[count]+line[count+1])!=",'" and count<len(line[:-2]):
				string=string[0:]+line[count]
				if line[count]=="," and line[count+1]==" ":
					param.append(string[0:-1])
					string=""
					count+=2
				else:
					count+=1
			param.append(string)
			#print(string)
			string=""
			list_param.append(param)
			param=[]
			line=f.readline(12)
		if len(list_param)==4 and method==1 and count==1:
			count=2
			while (line[count]+line[count+1]+line[count+2])!="], " and count<len(line[:-3]):
				while line[count]!="," and line[count]!="]" and count<len(line[:-3]):
					string=string[0:]+line[count]
					if line[count+1]==",":
						shell.append(string)
						string=""
						count+=2
					if line[count+1]=="]" and  count+5<len(line[:-1]):
						shell.append(string)
						param.append(shell)
						string=""
						shell=[]
						count+=5
					else:
						count+=1
			#print(string)
			shell.append(string)
			param.append(shell)
			string=""
			shell=""
			list_param.append(param)
			param=[]
			line=f.readline(12)
	print(list_param)
	return(list_param)

'''
This function sets default test parameters for zaef and peakfinder8 peak search methods.
'''
def set_param(method):
	#zaef method
	list_param=[]
	if method=="zaef":
		param=[500,5000,25000,50000,100000,150000]
		list_param.append(param)
		param=[2,3,4,5,6,8,10,15,20]
		list_param.append(param)
		param=["2,3,4","2,4,5","3,4,5","3,5,7","4,5,7"]
		list_param.append(param)
		param=[2,4,6,8,10]
		list_param.append(param)
		param=[2,4,6,8,10]
		list_param.append(param)
		param=[" "]
		list_param.append(param)
	if method=="peakfinder8":
		param=[2,3,4,5,6,8,10,15,20]
		list_param.append(param)
		param=[2,3,4,5,6,7,8,9]
		list_param.append(param)
		param=[5,10,20,50,100,200,300,500]
		list_param.append(param)
		param=[2,3,4,5,6,7,8,9,10,15,20]
		list_param.append(param)
		param=[[0,1200],[0,450],[450,1200]]
		list_param.append(param)
	f=open("peakopt_param.tab","w+")
	f=open("peakopt_param.tab","a+")
	count=1
	print("Default parameters:\n")
	f.write("Peakopt_steps: Zaef 1-sqrd 2-snr 3-radii 4-median filter 5- medianfilter + noise filter 6-noise filter\nPeakfinder8 1-snr 2-minpixcount 3-maxpixcount 4-local-bg-radius 5-min max resolution\n")
	for i in list_param:
		f.write("Peakopt_step_"+str(count)+"\n"+str(i)+"\n")
		print("Peakopt_step_"+str(count)+"\n"+str(i)+"\n")
		count+=1
	f.close()
	if method=="zaef":
		print("Peakopt_steps: Zaef 1-sqrd 2-snr 3-radii 4-median filter 5- medianfilter+noise filter 6-noise filter\n")
	if method=="peakfinder8":
		print("Peakopt_steps: Peakfinder8 1-snr 2-minpixcount 3-maxpixcount 4-local-bg-radius 5-min max resolution\n")

	default=input("Keep default parameters? y/n")
	if default=="n":
		print("Change parameters selected, don't forget to save your corrections on peakopt_param.tab") 
		sub.call("pluma peakopt_param.tab", shell=True)
		list_param=read_list_param(method)
	return(list_param)

def main():
	######## Start peak search optimization ###########
	
	inp="files.lst"
	geom="pilatus2mpanel.geom"
	out="calibration.stream"

	#Users communication disbabled.
	"""
	method=input("Select your peaksearch method: 1- zaef 2-peakfinder8 ")
	print("Peak search optmization: "+method+" method")
	if method=="1":
		method='zaef'
		print("Peak search optmization: "+method+" method")
		zaef_peakopt(inp,out,geom)
	if method=="2":
		method='peakfinder8'
		print("Peak search optmization: "+method+" method")
		peakfinder8_peakopt(inp,out,geom)
	"""


	method=2
	if method=="1":
		method='zaef'
		print("Peak search optmization: "+method+" method")
		zaef_peakopt(inp,out,geom)
	if method=="2":
		method='peakfinder8'
		print("Peak search optmization: "+method+" method")
		peakfinder8_peakopt(inp,out,geom)

	#Final curves comparison

	#curves=[['zaef',70,5000,10,[4,5,7],1],['zaef',70,5000,8,[2,3,4],0],['zaef',90,5000,5,[4,5,7],0],['peakfinder8',40,10,2,200,3],['peakfinder8',70,10,1,200,2],['peakfinder8',70,10,1,200,3]]
	#curves=[['zaef',80,5000,4,'2,3,4',0],['zaef',80,5000,4,'2,3,4',6],['zaef',100,5000,4,'3,4,5',0],['peakfinder8',80,5,2,50,6], ['peakfinder8',100,5,2,50,6]]
	curves=[['zaef',100,5000,4,'4,5,7',' ']]

	#Histograms, shift maps and unit cell parameters according to the image ID plots for each curve. Here no unit cell is included.
	"""
	finalpeakopt(inp,curves,geom,0)
	fileformat("output.tab","peakopt.tab",curves, [0])
	methods,total=filesearch_crystal("output.tab")
	for i in range(len(curves)):
		crystplots.shift_map(methods,total,i)
		crystplots.plot_hist(methods,total,i)
		crystplots.plot_crystals_cell(methods,total,i)
	sub.call("mkdir finalopt; mv peakopt.tab finalopt; mv *.png finalopt", shell=True)
	"""
	
	#Histograms, shift maps and unit cell parameters according to the image ID plots for each curve. Here unit cell is included, substitute name.cell with your unit cell file.

	"""
	finalpeakopt(inp,curves,geom,'name.cell')
	fileformat("output.tab","peakopt.tab",curves, [0])
	methods,total=filesearch_crystal("output.tab")
	
	for i in range(len(curves)):
		crystplots.shift_map(methods,total,i)
		crystplots.plot_hist(methods,total,i)
		crystplots.plot_crystals_cell(methods,total,i)
	sub.call("mkdir finalopt_cell; mv peakopt.tab finalopt; mv *.png finalopt_cell", shell=True)

	"""

	####### End peak search optimization #########
