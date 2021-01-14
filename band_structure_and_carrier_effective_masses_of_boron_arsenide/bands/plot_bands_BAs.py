#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

def parse_kpath(string):
	pts = string.split(',')
	labels = []
	xticks = []
	for pt in pts:
		blank,l,x = pt.split('"')
		if l.strip() == 'G':
			labels.append(r'$\Gamma$')
		elif l.strip() == 'Kla':
			labels.append(r'$K\leftarrow$')
		elif l.strip() == 'raL':
			labels.append(r'${\rightarrow}L$')
		elif l.strip() == 'raX':
			labels.append(r'${\rightarrow}X$')
		elif l.strip() == 'Gla':
			labels.append(r'$G\leftarrow$')
		else:
			labels.append(l.strip())
		xticks.append(float(x.strip()))
	return (xticks,labels)

def add_label(label,ax,loc):
	ax.annotate(label,xy=loc,xycoords='axes fraction',backgroundcolor='w',weight='bold')

def add_annotations(labels,ax,locs):
	for la,lo in zip(labels,locs):
		ax.annotate(la,xy=lo,xycoords='axes fraction',size=7)

def plot_bands(ax,paths,label_string,x_range,y_range,y_shift,formats,legend_data,outname):
	font = {'family':'sans-serif','size': 10}
	plt.rc('font', **font)

	bands = []
	for bandfile in paths:
		with open(bandfile,'r') as fin:
			xcoord = []
			eng = []
			appended = False
			for row in fin:
				try:
					x,e = row.split()
					xcoord.append(float(x))
					eng.append(float(e)-y_shift)
					appended = False
				except:
					bands.append([xcoord,eng])
					appended = True
					xcoord = []
					eng = []
					continue
			if not appended:
				bands.append([xcoord,eng])

	if formats == 'default':
		c='k'
		lw=0.95
		ls='-'
	
	for band,fmt in zip(bands,formats):
		if formats != 'default':
			c=fmt['c']
			lw=fmt['lw']
			ls=fmt['ls']
		ax.plot(band[0],band[1],c=c,lw=lw,ls=ls)

	# style
	xticks,labels = parse_kpath(label_string)
	for label in labels:
		if label == 'G':
			label = r'$\Gamma$'
	ax.set_ylabel('Energy (eV)')
	ax.grid(axis='x',c='k')

	ax.set_xticks(xticks)
	ax.set_xticklabels(labels)
	if x_range == 'default':
		ax.set_xlim(0, xticks[-1])
	else:
		ax.set_xlim(x_range)
	if y_range != 'default':
		ax.set_ylim(y_range)

	if legend_data:
		ax.legend(**legend_data)


if __name__ == "__main__":

	path_gw = ['gw/band{i}_gw.csv'.format(i=i) for i in range(1,9)]
	path_gw_so = ['gw_so/band{i}_so_gw.csv'.format(i=i) for i in range(1,17)]
	path_lda = ['lda/band{i}_lda.csv'.format(i=i) for i in range(1,9)]
	path_lda_so = ['lda_so/band{i}_so_lda.csv'.format(i=i) for i in range(1,17)]
	gw_all = path_gw_so + path_gw
	lda_all =  path_lda_so + path_lda
	
	formats_lda = [{'c':'xkcd:blue','lw':0.95,'ls':'-'}]*16+[{'c':'xkcd:royal blue','lw':0.95,'ls':'--'}]*8
	formats_gw = [{'c':'xkcd:red','lw':0.95,'ls':'-'}]*16+[{'c':'xkcd:deep red','lw':0.95,'ls':'--'}]*8
	
	gw_so_vbm = 7.34345955276091
	lda_so_vbm = 8.39021286737111
	zoom_range = (3.58258,4.08258)
	
	labels = '" G "  0.00000," X "  1.31513," W "  1.97270," K "  2.43767," G "  3.83258," L "  4.97152," U "  5.77687," W "  6.24184," L "  7.17178,"K/U"  7.97714," X "  8.44211'
	zoom_labels = '" Kla " {min}," G "  3.83258," raL "  {max}'.format(min=zoom_range[0],max=zoom_range[1])

	gw_legend_elements = [Line2D([0],[0],c='xkcd:deep red',ls='--',label='GW'), Line2D([0],[0],c='xkcd:red',ls='-',label='GW+SO')]
	gw_legend = {'loc':4,'handles':gw_legend_elements,'prop':{'size':8},'frameon':True,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.7}
	gw_zoom_legend = {'loc':8,'handles':gw_legend_elements,'prop':{'size':8},'frameon':True,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.7}

	lda_legend_elements = [Line2D([0],[0],c='xkcd:royal blue',ls='--',label='LDA'), Line2D([0],[0],c='xkcd:blue',ls='-',label='LDA+SO')]
	lda_legend = {'loc':4,'handles':lda_legend_elements,'prop':{'size':8},'frameon':True,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.7}
	lda_zoom_legend = {'loc':8,'handles':lda_legend_elements,'prop':{'size':8},'frameon':True,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.7}
	
	fig = plt.figure(figsize=(6.5,4))
	gs = gridspec.GridSpec(2,2,width_ratios = [2,1])
	gw_full_ax = fig.add_subplot(gs[1,0])
	lda_full_ax = fig.add_subplot(gs[0,0],sharex=gw_full_ax)
	plt.setp(lda_full_ax.get_xticklabels(), visible=False)
	gw_zoom_ax = fig.add_subplot(gs[1,1])
	lda_zoom_ax = fig.add_subplot(gs[0,1],sharex=gw_zoom_ax)
	plt.setp(lda_zoom_ax.get_xticklabels(), visible=False)

	minorLocator_full = MultipleLocator(1)
	minorLocator_zoom = MultipleLocator(0.1)

	plot_bands(lda_full_ax,lda_all,labels,'default',(-16,5.5),lda_so_vbm,formats_lda,None,'BAs_full_lda')
	lda_full_ax.yaxis.set_minor_locator(minorLocator_full)
	lda_full_ax.tick_params(axis='x',bottom=False)
	add_label('(a)',lda_full_ax,(-0.195,0.94))

	plot_bands(gw_full_ax,gw_all,labels,'default',(-16,6.5),gw_so_vbm,formats_gw,None,'BAs_full_gw')
	gw_full_ax.yaxis.set_minor_locator(minorLocator_full)
	add_label('(c)',gw_full_ax,(-0.195,0.94))
	ann_labels = [r'$\Delta$',r'$X_{6v,1}$',r'$X_{7v,1}$',r'$X_{6v,2}$',r'$X_{7v,2}$',
		r'$X_{6c}$',r'$X_{7c}$',r'$\Gamma_{6v}$',r'$\Gamma_{7v}$',r'$\Gamma_{8v}$',
		r'$\Gamma_{7c}$',r'$\Gamma_{6c}$',r'$L_{6v,1}$',r'$L_{6v,2}$',
		r'$L_{6v,3}$',r'$L_{4v},L_{5v}$',r'$L_{6c,1}$',r'$L_{6c,2}$',r'$L_{4c},L_{5c}$']
	xx = 1.003#0.976
	gx = 0.458 #0.436
	lx = 0.802
	ann_locs = [(0.10,0.73),(xx,0.22),(xx,0.34),(xx,0.46),(xx,0.55),(xx,0.75),(xx,0.84),
	(gx,0.06),(gx,0.613),(gx,0.74),(gx,0.83),(gx,1.02),
	(lx,0.092),(lx,0.27),(lx,0.56),(lx,0.645),(lx,0.82),(lx,0.926),(lx,1.02)]
	add_annotations(ann_labels,gw_full_ax,ann_locs)
	gw_full_ax.annotate(r'$\Gamma_{8c}$',xy=(0.453,0.906),xytext=(gx+0.055,1.02),xycoords='axes fraction',textcoords='axes fraction',
		size=7,arrowprops=dict(facecolor='black',shrink=0.000,width=0.5,headwidth=2.9,headlength=4.5))
	plot_bands(lda_zoom_ax,lda_all,zoom_labels,zoom_range,(-1.5,0.5),lda_so_vbm,formats_lda,lda_zoom_legend,'BAs_zoom_lda') 
	lda_zoom_ax.yaxis.set_minor_locator(minorLocator_zoom)
	lda_zoom_ax.tick_params(axis='x',bottom=False)
	add_label('(b)',lda_zoom_ax,(-0.43,0.94))
	
	plot_bands(gw_zoom_ax,gw_all,zoom_labels,zoom_range,(-1.5,0.5),gw_so_vbm,formats_gw,gw_zoom_legend,'BAs_zoom_gw') 
	gw_zoom_ax.yaxis.set_minor_locator(minorLocator_zoom)
	ann_labels = [r'$\Gamma_{7v}$',r'$\Gamma_{8v}$']
	ann_locs = [(0.51,0.57),(0.51,0.78)]
	add_annotations(ann_labels,gw_zoom_ax,ann_locs)
	add_label('(d)',gw_zoom_ax,(-0.43,0.94))
	
	plt.tight_layout()
	#gs.tight_layout(fig,w_pad=2)
	
	plt.savefig('BAs_full'+'.png', format='png',dpi=600)

