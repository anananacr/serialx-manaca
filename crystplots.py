#Importa pacotes relevantes
from mpl_toolkits.axes_grid1 import host_subplot
import math
import re
import scipy
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
from colorspacious import cspace_converter
from collections import OrderedDict
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy import stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import csv
import subprocess as sub
import seaborn as sns
import mpl_toolkits.axisartist as AA
#from PIL import Image, ImageDraw
import os
import matplotlib.image as mpimg

def remove_repetidos(lista):
	l = []
	for i in lista:
		if i not in l:
			l.append(i)
	l.sort()
	return l

def normal(mean, std, color="black"):
    x = np.linspace(mean-4*std, mean+4*std, 200)
    p = stats.norm.pdf(x, mean, std)
    z = plt.plot(x, p, color, linewidth=1, ls='-')


def calc_vol(ra,rb,rc,ralf,rbet,rgam):
	#if abc in A
	fac_a=0.1
	#fac_a=1 in nm
	#if alf bet gam in deg
	fac=math.pi/180
	#fac=1 in rad
	ra=fac_a*ra
	rb=fac_a*rb
	rc=fac_a*rc
	ralf=fac*ralf
	rbet=fac*rbet
	rgam=fac*rgam
	vol=ra*rb*rc*math.sqrt(1-pow(math.cos(ralf),2)-pow(math.cos(rbet),2)-pow(math.cos(rgam),2) +2*math.cos(ralf)*math.cos(rbet)*math.cos(rgam))
	return(vol)

def plot(method,rate,labels, labels_par):
	step=len(labels[0])
    #sns.set_style("darkgrid", {"axes.facecolor": ".9"})
	#color=sns.color_palette("Paired")
	#sns.set_palette(color)
	sns.set_context("paper")
    #plot index rate vs threshold for each par (sqrd or snr)
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)

    #upa o dado na mesma pasta e troca o nome do arquivo
	data_dir = "./"
	data = os.path.join(data_dir, 'peakopt.tab')
    #carrega dados
	df= pd.read_csv(data, delimiter=' ', skiprows=0,usecols=(0,1,2,3,4,5,6,7))
	thr=df['thr'].tolist()
	param=df['sqrd'].tolist()
    #rate=df['index'].tolist()
	score=df[rate].tolist()
	#ver quais curvas vai ter que plotar quantos valores de sqrd/snr foram testados
	curves=remove_repetidos(param)
    #print(curves)
	fig = plt.figure(figsize=(8, 7), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	if rate=="mpp":
		ax.set_ylabel('Mean peaks per pattern', fontsize=10)
	if rate=="map":
		ax.set_ylabel('Mean ADU intensity per peak', fontsize=10)
	if rate=="tap":
		ax.set_ylabel('Total ADU intensity per pattern', fontsize=10)
	if rate=='crystals':
		ax.set_ylabel('# Crystals indexed', fontsize=10)
	if rate=="index":
		ax.set_ylabel('Indexing rate (overall) %', fontsize=10)
	if rate=="time":
		ax.set_ylabel('Processing time', fontsize=10)
	if method=='zaef':
		sup_titles=['Sqrd opt', 'Snr opt', 'Radii opt']
		x_label=['Thr', 'Best thr/sqrd option', "Best thr/sqrd/snr option", "Best thr/sqrd/snr/radii"]
	if method=='peakfinder8':
		sup_titles=['Snr opt', 'Min pixel count opt', 'Max pixel count opt', 'Local bg radius opt', 'Resol opt']
		x_label=['Thr', 'Best thr/snr option', "Best thr/snr/minpix option", "Best thr/snr/minpix/maxpix", "Best thr/snr/minpix/maxpix/lbg"]
	ax.set_xlabel(x_label[step-2], fontsize=10)
	fig.suptitle(sup_titles[step-2], fontsize=16)
	lines=[]
	titles=[]
	colors=['#1f1fff', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#918e8e', '#e41a1c', '#dede00', '#2b2929', '#377eb8','#c9c9c9'  ]
	for j in range(len(curves)):
		pointsy=[]
		pointsx=[]
		for i in range(len(param)):
			if param[i]==curves[j]:
				pointsy.append(score[i])
				pointsx.append(thr[i])
		plt.plot(pointsx, pointsy, '--o',color=colors[j],label=str(labels_par[j]))
	plt.legend()
    
    
	plt.savefig("plot_"+rate+".png")
    #sub.call("mkdir peakopt_"+str(method)+"_"+str(step-1)+"/; mv *.tab peakopt_"+str(method)+"_"+str(step-1)+"/; mv *.png peakopt_"+str(method)+"_"+str(step-1)+"/;", shell=True)

def plot_thr(method,rate,labels, labels_param):
	step=len(labels[0])
	art=0
    #plot index rate vs par (sqrd or snr) for each threshold
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)
	sns.set_context("paper")
    #upa o dado na mesma pasta e troca o nome do arquivo
	data_dir = "./"
	data = os.path.join(data_dir, 'peakopt.tab')
    #carrega dados
	df= pd.read_csv(data, delimiter=' ', skiprows=0,usecols=(0,1,2,3,4,5,6,7))
	thr=df['thr'].tolist()
	param=df['sqrd'].tolist() 
    #rate=df['index'].tolist()
	score=df[rate].tolist()
	#ver quais curvas vai ter que plotar quantos valores de sqrd/snr foram testados
	curves=remove_repetidos(thr)
	fig = plt.figure(figsize=(8, 7), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	if rate=="mpp":
		ax.set_ylabel('Mean peaks per pattern', fontsize=10)
	if rate=="map":
		ax.set_ylabel('Mean ADU intensity per peak', fontsize=10)
	if rate=="tap":
		ax.set_ylabel('Mean ADU intensity total per pattern', fontsize=10)
	if rate=='crystals':
		ax.set_ylabel('# Crystals indexed', fontsize=10)
	if rate=="index":
		ax.set_ylabel('Indexing rate (overall) %', fontsize=10)
	if rate=="time":
		ax.set_ylabel('Processing time', fontsize=10)
	if method=='zaef':
		sup_titles=['Sqrd opt', 'Snr opt', 'Radii opt']
		x_label=['Sqrd', 'Snr', "Radii option"]
	if method=='peakfinder8':
		sup_titles=['Snr opt', 'Min pixel count opt', 'Max pixel count opt', 'Local bg radius opt', 'Resol opt']
		x_label=['Snr', 'Min pix count', "Max pix count", "Local bg radius", "Resol options"]
	ax.set_xlabel(x_label[step-2], fontsize=10)
	fig.suptitle(sup_titles[step-2], fontsize=16)
	lines=[]
	titles=[]
	colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']
    
	for j in range(len(curves)):
		pointsy=[]
		pointsx=[]
		for i in range(len(param)):
			if thr[i]==curves[j]:
				pointsy.append(score[i])
				pointsx.append(param[i])
		plt.plot(pointsx, pointsy, '--o', color=colors[j], label=str(labels[j]))
	plt.legend()
	plt.savefig("plot_thr_"+rate+".png")
	sub.call("mkdir peakopt_"+str(method)+"_"+str(step-1)+"/; mv *t.tab peakopt_"+str(method)+"_"+str(step-1)+"/; mv *.png peakopt_"+str(method)+"_"+str(step-1)+"/;", shell=True)

def plot_median(rate,labels,labels_par):
	#plot index rate vs threshold for each par (sqrd or snr)
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)
	sns.set_context("paper")
	#upa o dado na mesma pasta e troca o nome do arquivo
	data_dir = "./"
	fig = plt.figure(figsize=(12, 10), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	folder=['peakopt_zaef_4',"peakopt_zaef_5", "peakopt_zaef_6"]
	title=["MF", "MF+NOI", "NOI"]
	marker=["--o", "--^", "--s"]
	count=0
	for k in folder:
		data = os.path.join(data_dir+k, 'peakopt.tab')
	#carrega dados
		df= pd.read_csv(data, delimiter=' ', skiprows=0,usecols=(0,1,2,3,4,5,6,7))
		thr=df['thr'].to_list()
		param=df['sqrd'].to_list()
		score=df[rate].to_list()
	#ver quais curvas vai ter que plotar quantos valores de sqrd/snr foram testados
		curves=remove_repetidos(param)
	#print(curves)
		if rate=="mpp":
			ax.set_ylabel('Mean peaks per pattern', fontsize=10)
		if rate=="map":
			ax.set_ylabel('Mean ADU intensity per peak', fontsize=10)
		if rate=="tap":
			ax.set_ylabel('Mean ADU intensity total per pattern', fontsize=10)
		if rate=='crystals':
			ax.set_ylabel('# Crystals indexed', fontsize=10)
		if rate=="index":
			ax.set_ylabel('Indexing rate (overall) %', fontsize=10)
		if rate=="time":
			ax.set_ylabel('Processing time', fontsize=10)
		ax.set_xlabel('Best thr/sqrd/snr/r', fontsize=10)
		fig.suptitle("Median filter and noise filter optmization", fontsize=16)
		colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']
		for j in range(len(curves)):
			pointsy=[]
			pointsx=[]
			for i in range(len(param)):
				if param[i]==curves[j]:
					pointsy.append(score[i])
					pointsx.append(thr[i])
			plt.plot(pointsx, pointsy, marker[count], color=colors[j], markersize=10, label=str(labels_par[j])+" "+title[count])
		count+=1
	plt.legend()


	plt.show()

def plot_median_thr(rate,labels, labels_par):
	#plot index rate vs par (sqrd or snr) for each threshold
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)
	sns.set_context("paper")
	#upa o dado na mesma pasta e troca o nome do arquivo
	data_dir = "./"
	folder=['peakopt_zaef_4',"peakopt_zaef_5", "peakopt_zaef_6"]
	title=["MF", "MF+NOI", "NOI"]
	marker=["--o", "--^", "--s"]
	count=0
	fig = plt.figure(figsize=(12, 10), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	for k in folder:
		data = os.path.join(data_dir+k, 'peakopt.tab')
		#carrega dados
		df= pd.read_csv(data, delimiter=' ', skiprows=0,usecols=(0,1,2,3,4,5,6,7))
		thr=df['thr'].to_list()
		param=df['sqrd'].to_list()
		score=df[rate].to_list()
    #ver quais curvas vai ter que plotar quantos valores de sqrd/snr foram testados
		curves=remove_repetidos(thr)
		if rate=="mpp":
			ax.set_ylabel('Mean peaks per pattern', fontsize=10)
		if rate=="map":
			ax.set_ylabel('Mean ADU intensity per peak', fontsize=10)
		if rate=="tap":
			ax.set_ylabel('Total ADU intensity per pattern', fontsize=10)
		if rate=='crystals':
			ax.set_ylabel('# Crystals indexed', fontsize=10)
		if rate=="index":
			ax.set_ylabel('Indexing rate (overall) %', fontsize=10)
		if rate=="time":
			ax.set_ylabel('Processing time', fontsize=10)
		ax.set_xlabel("MF and noise", fontsize=10)
		fig.suptitle("Median and noise filter optmization", fontsize=16)
		colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']

		for j in range(len(curves)):
			pointsy=[]
			pointsx=[]
			for i in range(len(param)):
				if thr[i]==curves[j]:
					pointsy.append(score[i])
					pointsx.append(param[i])
			plt.plot(pointsx, pointsy, marker[count], color=colors[j], markersize=10, label=str(labels[j])+" "+title[count])
		count+=1
	plt.legend()


	plt.show()

def plot_crystals_cell(methods,total, opt):
	param=['a','b','c','alf','bet','gam']
	colors=['r','b','g','darkorange','darkviolet', 'grey']
	fig = plt.figure(figsize=(9, 7), tight_layout=True)
	k=0
	for j in param:
		pointsx=[]
		pointsy=[]
		for i in range(int(total[opt])-1):
			pointsx.append(methods[opt][i].id)
			if j=='a':
				pointsy.append(methods[opt][i].a)
				unit='A'
				k=0
			if j=='b':
				pointsy.append(methods[opt][i].b)
				unit='A'
				k=1
			if j=='c':
				pointsy.append(methods[opt][i].c)
				unit='A'
				k=2
				plt.title("Cell parameters serial image")
			if j=='alf':
				pointsy.append(methods[opt][i].alf)
				unit='deg'
				k=3
			if j=='bet':
				pointsy.append(methods[opt][i].bet)
				unit='deg'
				k=4
			if j=='gam':
				pointsy.append(methods[opt][i].gam)
				unit='deg'
				k=5
		ax = fig.add_subplot(2,3,k+1)
		ax.set_xlabel('Crystal ID', fontsize=10)
		ax.set_ylabel(j+" ("+unit+")", fontsize=10)
		plt.scatter(pointsx,pointsy,color=colors[k], marker='.',linewidths=1)
	plt.savefig('cell_id_'+str(opt)+'.png')
	plt.close()

def plot_crystals_vol(methods,total,labels):
	colors=['r','b','g','darkorange','darkviolet', 'grey']
	fig = plt.figure(figsize=(9, 7), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	ax.set_xlabel('Crystal ID', fontsize=10)
	ax.set_ylabel("Unit cell volume (nm^3)", fontsize=10)
	plt.title("Unit cell volume serial image")
	k=0
	
	for j in methods:
		pointsx=[]
		pointsy=[]
		for i in range(int(total[k])-1):
			pointsx.append(j[i].id)
			ra=(j[i].a)
			rb=(j[i].b)
			rc=(j[i].c)
			ralf=j[i].alf
			rbet=j[i].bet
			rgam=j[i].gam
			vol=calc_vol(ra,rb,rc,ralf,rbet,rgam)
			vol=round(vol,5)
			pointsy.append(vol)
		#print(pointsx[0:10],pointsy[0:10])
		plt.scatter(pointsx,pointsy,color=colors[k],marker=".", label=labels[k])
		k+=1
	plt.legend()
	#plt.show()
	plt.savefig('vol_id_ringsnocen.png')
	#plt.close()

def shift_map(methods,total, opt):
	pointsx=[]
	pointsy=[]
	occur=[]
	for i in range(int(total[opt])-1):
		pointsx.append(methods[opt][i].detshift[0])
		pointsy.append(methods[opt][i].detshift[1])
	#det-shift code
	nbins = 200
	H, xedges, yedges = np.histogram2d(pointsx,pointsy,bins=nbins)
	H = np.rot90(H)
	H = np.flipud(H)
	Hmasked = np.ma.masked_where(H==0,H)

	fig2 = plt.figure()
	cmaps=OrderedDict()
	plt.pcolormesh(xedges,yedges,Hmasked, cmap='plasma',linewidth=2)
	plt.title('Detector shifts according to prediction refinement')
	plt.xlabel('x shift / mm')
	plt.ylabel('y shift / mm')
	plt.plot(0, 0, 'bH', color='c')
	fig = plt.gcf()
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Counts')
	#plt.show()
	plt.savefig("detshift_"+str(opt)+".png")

def plot_hist(methods,total, opt):
	#print(methods, total)
	pointsx=[]
	k=0
	nbins=200
	param=['a','b','c','alf','bet','gam']
	title_param=['a','b','c','α', 'β', 'γ']
	colors=['r','b','g','darkorange','darkviolet', 'grey']
	#ame peak
	#peak=[[37,38.5],[45,47],[78,80],[89.5,91],[89,91],[99,103]]
	#liso peak triclinic
	peak=[[77,80],[76,80],[36,40],[88.5,91.5],[88.5,91.5],[88.5,91.5]]
	#liso peak tetragonal
	#peak=[[60,80],[60,80],[20,40],[88,92],[88,92],[88,92]]

	#axes=[[0,150],[0,120],[0,150],[60,120],[60,120],[60,120]]
	axes=[[75,83],[75,83],[34,38],[88,92],[88,92],[88,92]]
	#axes=[[35,45],[75,85],[40,50],[85,95],[85,95],[75,115]]
	#ylim=[0.6, 0.6, 0.8, 0.5, 0.5, 0.5]
	fig = plt.figure(figsize=(9, 7), tight_layout=True)
	count=0
	for j in param:
		pointsx=[]
		pointsp=[]
		indexed=0
		selected=0
		for i in range(int(total[opt])-1):
			indexed+=1
			if j=='a':
				k=0
				unit='Å'
				pointsx.append(methods[opt][i].a)
				if peak[k][0]<(methods[opt][i].a)<peak[k][1]:
					pointsp.append(methods[opt][i].a)
					selected+=1
			if j=='b':
				k=1
				pointsx.append(methods[opt][i].b)
				unit='Å'
				if peak[k][0]<(methods[opt][i].b)<peak[k][1]:
					pointsp.append(methods[opt][i].b)
					selected+=1
			if j=='c':
				k=2
				pointsx.append(methods[opt][i].c)
				if peak[k][0]<(methods[opt][i].c)<peak[k][1]:
					pointsp.append(methods[opt][i].c)
					selected+=1
				plt.title("Unit cell parameters")
				unit='Å'
			if j=='alf':
				k=3
				pointsx.append(methods[opt][i].alf)
				unit='degrees'
				if peak[k][0]<(methods[opt][i].alf)<peak[k][1]:
					pointsp.append(methods[opt][i].alf)
					selected+=1
			if j=='bet':
				k=4
				pointsx.append(methods[opt][i].bet)
				unit='degrees'
				if peak[k][0]<(methods[opt][i].bet)<peak[k][1]:
					pointsp.append(methods[opt][i].bet)
					selected+=1
			if j=='gam':
				k=5
				pointsx.append(methods[opt][i].gam)
				unit='degrees'
				if peak[k][0]<(methods[opt][i].gam)<peak[k][1]:
					pointsp.append(methods[opt][i].gam)
					selected+=1
		ax = fig.add_subplot(2,3,k+1)
		'''
		n,bins,_=plt.hist(pointsx,nbins,color=colors[k], alpha=0.85, density=True, histtype='stepfilled')
		mu, sigma = scipy.stats.norm.fit(pointsx)
		best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
		#fact=max(n)/max(best_fit_line)
		fact=1
		plt.plot(bins, fact*best_fit_line, '--k')
		'''
		mu = np.mean(pointsp)
		sigma = np.std(pointsp)
		sns.histplot(pointsx,stat='density', color=colors[k], kde=True)
		normal(mu, sigma)
		plt.xlabel(title_param[count]+" ("+unit+")")
		plt.gca()
		#plt.ylabel('Counts')
		plt.text(0.3,0.3,r'$'+str(np.round(mu,2))+'\pm'+str(np.round(sigma,2))+'$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=14)
		ax.set_xlim(left=axes[k][0],right=int(axes[k][1]))
		#ax.set_ylim(0,float(ylim[k]))
		print('Indexed: '+str(indexed)+' Selected for fitting: '+str(selected))
		count+=1
	plt.savefig("hist_all_"+str(opt)+".png")
	plt.close()

def plot_idx(x_labels, step):

	majorLocator = MultipleLocator(2)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(1)

	rx = re.compile(r'peakopt\.(tab)')
	r = []
	for path, dnames, fnames in os.walk('./idx_'+str(step)):
    		r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])

	count=0
	rate=["index","mpp", "map"]
	ylab=['Indexing rate (% of hits)','MPP', 'MAP']
	fig= plt.figure(figsize=(10, 6), tight_layout=True)
	host = host_subplot(111,axes_class=AA.Axes)
	plt.subplots_adjust(right=0.75)
	par = host.twinx()
	par2 = host.twinx()
	new_fixed_axis = par2.get_grid_helper().new_fixed_axis
	offset = 60
	par2.axis["right"] = new_fixed_axis(loc="right",axes=par2,offset=(offset, 0))
	par2.axis["right"].toggle(all=True)
	new_fixed_axis = par.get_grid_helper().new_fixed_axis
	par.axis["right"].toggle(all=True)
	labels=x_labels

	host.set_xlabel("method ID", size='x-large', labelpad=25)
	host.set_ylabel(ylab[0], size='x-large')
	par.set_ylabel(ylab[1], size='x-large')
	par2.set_ylabel(ylab[2], size='x-large')

	fig.subplots_adjust(top=0.8)
	colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']
	for l in range(len(rate)):
		curve=[]
		point=[]
		for k in r:
			data = os.path.join(k)
			#carrega dados
			df= pd.read_csv(data, delimiter=' ', skiprows=0,usecols=(0,1,2,3,4,5,6))
			thr=df['thr'].to_list()
			score=df[rate[l]].to_list()
    			#ver quais curvas vai ter que plotar quantos valores de sqrd/snr foram testados
			#ax.set_ylabel(ylab[l], fontsize=10)
			fig.suptitle("Indexing methods", size='x-large')
			for j in range(len(score)):
				point.append(thr[j])
				point.append(score[j])
				curve.append(point)
				point=[]
		curve= sorted(curve, key=lambda y: y[0])
		a = [x for x,y in curve]
		b = [y for x,y in curve]
		if l==0:
			p1, =host.plot(a,b, "--o", color=colors[l], markersize=5, label=rate[l])
			host.axis["left"].label.set_color(p1.get_color())
			host.axis["left"].label.set_fontsize(16)
			host.axis["bottom"].label.set_fontsize(16)
		if l==1:
			p2, =par.plot(a,b, "--o", color=colors[l], markersize=5, label=rate[l])
			par.axis["right"].label.set_color(p2.get_color())
			par.axis["right"].label.set_fontsize(16)
		if l==2:
			p3, =par2.plot(a,b, "--o", color=colors[l], markersize=5, label=rate[l])
			par2.axis["right"].label.set_color(p3.get_color())
			par2.axis["right"].label.set_fontsize(16)
	ax=plt.gca()
	#plt.legend(fontsize=12)
	plt.xticks(a,labels, fontsize=16)
	host.set_ylim(0,100)

	plt.margins(0.7)
	plt.setp(host.axis["bottom"].major_ticklabels, rotation=45, rotation_mode='anchor', size='large')
	plt.setp(par.axis["right"].major_ticklabels, size='large')
	plt.setp(par2.axis["right"].major_ticklabels, size='large')
	plt.setp(host.axis["left"].major_ticklabels, size='large')
	#plt.show()
	plt.savefig("indx_"+str(step)+".png")

def plot_check(labels,method,run,dir):
	#plot index rate vs threshold for each par (sqrd or snr)
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)
	sns.set_context("paper")
	#upa o dado na mesma pasta e troca o nome do arquivo
	data_dir = dir
	fig = plt.figure(figsize=(10, 7), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	y_label=["Completeness (%)", 'SNR']
	#rate=['refs', 'Meas']
	rate=['refs', 'Meas']
	name=["compl","snr"]
	'''
	if method=='part':
		rx = re.compile(r'_'+run+'\.(dat)')
	if method=='proc':
		rx = re.compile(r'0\.(dat)')

	r = []
	for path, dnames, fnames in os.walk(data_dir):
		r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])
	#print(r)
	#r=sorted(r)
	#print(r)
	'''
	#colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']
	param=0
	for i in rate:
		fig = plt.figure(figsize=(10, 7), tight_layout=True)
		ax = fig.add_subplot(1,1,1)
		ax2 = ax.twiny()
		count=0
		for j in dir:
			print(count)
			data_dir=j
			if method=='part':
				rx = re.compile(r'_'+run+'\.(dat)')
			if method=='proc':
				rx = re.compile(r'0\.(dat)')
			r = []
			for path, dnames, fnames in os.walk(data_dir):
				r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])

			tests=[]
			number=""
			f=open(data_dir+'resume'+name[param]+'.txt','w+')
			for k in r:
				if method=='part':
					number=k[-8:-6]
				if method=='proc':
					number=k[-6:-4]

				if number[0]=="_":
					tests.append(number[1:])
				else:
					tests.append(number)

				data_path = os.path.join(k)
				data= pd.read_csv(data_path,delimiter=' ', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), skipinitialspace=True)
				print(data)
				df=pd.DataFrame(data, columns=['Center','Std',i])
				res=df['Std'].to_list()
				d=df['Center'].to_list()
				score=df[i].to_list()
				new=[]
				for l in d:
					new.append(l/10)
				d=new
				print(score)

				f=open(data_dir+'resume'+name[param]+'.txt','a+')
				ax.plot(d,score,marker='o', linestyle='-', label=labels[count][0])
				f.write(str(tests[-1])+" "+str(labels[int(tests[-1])])+'\n')
				ax.set_ylabel(y_label[param], fontsize=16)
				ax.set_xlabel('1/d (1/Å)', fontsize=16)
				#ax.set_xlim(0.15,0.75)
				ax.yaxis.label.set_size(16)
				ax2.set_xlim(ax.get_xlim())
				ax2.set_xticks(d)
				ax2.set_xticklabels(res)
				ax.tick_params(axis='both', which='major', labelsize=16)
				ax.tick_params(axis='both', which='minor', labelsize=12)
				ax2.set_xlabel("Resolution (Å)",fontsize=14)
				ax2.tick_params(axis='both', which='major', labelsize=12)
				ax2.tick_params(axis='both', which='minor', labelsize=10)
				"""

				if (score[0])>=90 and param==0:
					plt.plot(res,score,marker='o', linestyle='-', label=count)
					f.write(str(count)+" "+str(labels[count]))
					ax.set_ylabel(y_label[param], fontsize=10)
					ax.set_xlabel('Resolution (A)', fontsize=10)
				"""


				f.close()
			count+=1
		#print(tests)
		ax.legend(fontsize=16)
		#plt.ylim(0,100)
		#plt.xlim(0,2)
		plt.savefig(name[param]+"_liso_"+str(run)+".png")
		plt.close()
		#plt.show()
		param+=1
	sub.call("mv *.png "+data_dir+";mv *resume*.txt "+data_dir, shell=True)

def plot_compare(labels,method,run, dir):
	#plot index rate vs threshold for each par (sqrd or snr)
	majorLocator = MultipleLocator(20)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(5)
	sns.set_context("paper")
	#upa o dado na mesma pasta e troca o nome do arquivo
	#data_dir = dir
	fig = plt.figure(figsize=(10,7), tight_layout=True)
	ax = fig.add_subplot(1,1,1)
	merit=['Rsplit','CC','CCstar']
	#colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00', 'darkcyan','maroon','indigo','magenta', 'seagreen', 'saddlebrown','yellow','springgreen','coral']
	for i in merit:
		fig = plt.figure(figsize=(10,7), tight_layout=True)
		ax = fig.add_subplot(1,1,1)
		ax2 = ax.twiny()
		count=0
		for j in dir:
			data_dir=j
			if method=='proc':
				rx = re.compile(r'_'+i+'\.(dat)')
			if method=='part':
				rx = re.compile(r'_'+run+'_'+i+'\.(dat)')
			r = []
			for path, dnames, fnames in os.walk(data_dir):
				r.extend([os.path.join(path, x) for x in fnames if rx.search(x)])
			print(r)
			tests=[]
			number=""
			#fig = plt.figure(figsize=(9, 7), tight_layout=True)
			#ax = fig.add_subplot(1,1,1)
			#ax2 = ax.twiny()
			f=open(data_dir+'resume'+i+'.txt','w+')
			for k in r:
				print(k)
				if method=='part':
					number=k[-8-len(i)-1:-6-len(i)-1]

				if method=='proc':
					number=k[-6-len(i)-1:-4-len(i)-1]

				if number[0]=="_":
					tests.append(number[1:])
				else:
					tests.append(number)

				#print(tests)
				data_path = os.path.join(k)
				data= pd.read_csv(data_path, delimiter=' ', usecols=(0,1,2,3), skipinitialspace=True)
				print(data)
				df=pd.DataFrame(data, columns=['1/d','centre', 'nref'])
				df=df.fillna(0)
				res=df['nref'].to_list()
				d=df['1/d'].to_list()
				score=df['centre'].to_list()
				new=[]
				for l in d:
					new.append(l/10)
				d=new
				print(d,res,score)
				f=open(data_dir+'resume'+i+'.txt','a+')
				#ax.plot(d,score, color=colors[count],marker='o', linestyle='-', label=labels[int(tests[-1])])
				ax.plot(d,score,marker='o', linestyle='-', label=labels[count][0])
				#f.write(str(tests[-1])+" "+str(labels[int(tests[-1])])+'\n1/d\n'+str(d)+'\nScore\n'+str(score))
				ax.set_ylabel(i, fontsize=16)
				ax.set_xlabel('1/d (1/Å)', fontsize=16)
				ax.yaxis.label.set_size(16)
				#ax.set_xlim(0.15,0.75)
				ax2.set_xlim(ax.get_xlim())
				ax2.set_xticks(d)
				ax2.set_xticklabels(res)
				ax.tick_params(axis='both', which='major', labelsize=16)
				ax.tick_params(axis='both', which='minor', labelsize=12)
				ax2.set_xlabel("Resolution (Å)",fontsize=16)
				ax2.tick_params(axis='both', which='major', labelsize=12)
				ax2.tick_params(axis='both', which='minor', labelsize=10)
				"""
				if (score[0])>=90 and param==0:
					plt.plot(res,score,marker='o', linestyle='-', label=count)
					f.write(str(count)+" "+str(labels[count]))
					ax.set_ylabel(y_label[param], fontsize=10)
					ax.set_xlabel('Resolution (A)', fontsize=10)
				"""
				f.close()
			count+=1
		#print(tests)
		ax.legend()
		if i=='CCstar':
			ax.set_ylabel('CC*', fontsize=16)
		if i=='Rsplit':
			plt.ylim(0,1000)
			#plt.xlim(0.3,0.75)
		
		plt.savefig(i+"_liso_compare_"+str(run)+".png")
		plt.close()
		#plt.show()
	sub.call("mv *_compare*.png "+data_dir+";mv *resume*.tab "+data_dir, shell=True)

def compare_hist(methods,total):
        #print(methods, total)
        pointsx=[]
        pointsy=[]
        k=0
        nbins=200
        param=['a','b','c','alf','bet','gam']
        colors_2=['r','b','g','orange','violet', 'silver']
        colors=['firebrick','navy','darkgreen','darkorange','darkviolet','dimgray']
        #ame peak
        #peak=[[35,42],[78,80],[45,55],[88,92],[95,108],[88,92]]
        #liso peak triclinic
        #peak=[[20,60],[60,80],[60,80],[88,92],[88,92],[88,92]]
        #liso peak tri
        peaka=[[37,38.5],[45,47],[78,80],[89.5,91],[89,91],[99,103]]
        peakb=[[37,38.5],[45,47],[78,80],[89.5,91],[89,91],[99,103]]

        axes=[[35,40],[40,50],[75,80],[85,95],[85,95],[95,110]]
        #axes=[[70,90],[70,90],[30,50],[85,95],[85,95],[85,95]]
        #ylim=[0.6, 0.6, 0.8, 0.5, 0.5, 0.5]
        fig = plt.figure(figsize=(9, 7), tight_layout=True)
        for j in param:
                pointsx=[]
                pointsy=[]
                pointsa=[]
                pointsb=[]
                for i in range(int(total[0])-1):
                        if j=='a':
                                k=0
                                unit='Å'
                                pointsx.append(methods[0][i].a)
                                if peaka[k][0]<(methods[0][i].a)<peaka[k][1]:
                                        pointsa.append(methods[0][i].a)
                        if j=='b':
                                k=1
                                pointsx.append(methods[0][i].b)
                                unit='Å'
                                if peaka[k][0]<(methods[0][i].b)<peaka[k][1]:
                                        pointsa.append(methods[0][i].b)
                        if j=='c':
                                k=2
                                pointsx.append(methods[0][i].c)
                                if peaka[k][0]<(methods[0][i].c)<peaka[k][1]:
                                        pointsa.append(methods[0][i].c)
                                plt.title("Unit cell parameters")
                                unit='Å'
                        if j=='alf':
                                k=3
                                pointsx.append(methods[0][i].alf)
                                unit='degrees'
                                if peaka[k][0]<(methods[0][i].alf)<peaka[k][1]:
                                        pointsa.append(methods[0][i].alf)
                        if j=='bet':
                                k=4
                                pointsx.append(methods[0][i].bet)
                                unit='degrees'
                                if peaka[k][0]<(methods[0][i].bet)<peaka[k][1]:
                                        pointsa.append(methods[0][i].bet)
                        if j=='gam':
                                k=5
                                pointsx.append(methods[0][i].gam)
                                unit='degrees'
                                if peaka[k][0]<(methods[0][i].gam)<peaka[k][1]:
                                        pointsa.append(methods[0][i].gam)
                for i in range(int(total[1])-1):
                        if j=='a':
                                k=0
                                unit='Å'
                                pointsy.append(methods[1][i].a)
                                if peakb[k][0]<(methods[1][i].a)<peakb[k][1]:
                                        pointsb.append(methods[1][i].a)
                        if j=='b':
                                k=1
                                pointsy.append(methods[1][i].b)
                                unit='Å'
                                if peakb[k][0]<(methods[1][i].b)<peakb[k][1]:
                                        pointsb.append(methods[1][i].b)
                        if j=='c':
                                k=2
                                pointsy.append(methods[1][i].c)
                                if peakb[k][0]<(methods[1][i].c)<peakb[k][1]:
                                        pointsb.append(methods[1][i].c)
                                plt.title("Unit cell parameters")
                                unit='Å'
                        if j=='alf':
                                k=3
                                pointsy.append(methods[1][i].alf)
                                unit='degrees'
                                if peakb[k][0]<(methods[1][i].alf)<peakb[k][1]:
                                        pointsb.append(methods[1][i].alf)
                        if j=='bet':
                                k=4
                                pointsy.append(methods[1][i].bet)
                                unit='degrees'
                                if peakb[k][0]<(methods[1][i].bet)<peakb[k][1]:
                                        pointsb.append(methods[1][i].bet)
                        if j=='gam':
                                k=5
                                pointsy.append(methods[1][i].gam)
                                unit='degrees'
                                if peakb[k][0]<(methods[1][i].gam)<peakb[k][1]:
                                        pointsb.append(methods[1][i].gam)
                ax = fig.add_subplot(2,3,k+1)
                '''
                n,bins,_=plt.hist(pointsx,nbins,color=colors[k], alpha=0.85, density=True, histtype='stepfilled')
                mu, sigma = scipy.stats.norm.fit(pointsx)
                best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
                #fact=max(n)/max(best_fit_line)
                fact=1
                plt.plot(bins, fact*best_fit_line, '--k')
                '''
                mu = np.mean(pointsa)
                sigma = np.std(pointsa)
                sns.histplot(pointsx,stat='density', color=colors[k], kde=True)
                sns.histplot(pointsy,stat='density', color=colors_2[k], kde=False)
                #normal(mu, sigma)
                plt.text(0.8,0.8,r'$'+str(np.round(mu,2))+'\pm'+str(np.round(sigma,2))+'$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=14)
                mu = np.mean(pointsb)
                sigma = np.std(pointsb)
                #normal(mu, sigma)
                plt.xlabel(j+" ("+unit+")")
                plt.gca()
                plt.text(0.3,0.3,r'$'+str(np.round(mu,2))+'\pm'+str(np.round(sigma,2))+'$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=14)
                ax.set_xlim(left=axes[k][0],right=int(axes[k][1]))
                plt.savefig("hist_comp.png")
        plt.close()
