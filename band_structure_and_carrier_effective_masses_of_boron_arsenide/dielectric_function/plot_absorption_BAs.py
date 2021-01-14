#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

def plot_data(ax,xdata,ydata,x_range,y_range,x_label,y_label,formats,legend_data):
	font = {'family':'sans-serif','size': 10}
	plt.rc('font', **font)

	if formats == 'default':
		c='k'
		lw=0.95
		ls='-'
	else:
		c=formats['c']
		lw=formats['lw']
		ls=formats['ls']
	
	ax.plot(xdata,ydata,c=c,lw=lw,ls=ls)

	# style
	ax.set_xlabel(x_label)
	ax.set_ylabel(y_label)

	if x_range != 'default':
		ax.set_xlim(x_range)
	if y_range != 'default':
		ax.set_ylim(y_range)

	if legend_data:
		ax.legend(**legend_data)


if __name__ == "__main__":

	data = pd.read_csv('abs_18.dat',delim_whitespace=True,header=None,names=['eng','eps2_noeh','eps2_eh'])
	formats_gw_noeh = {'c':'k','lw':1.05,'ls':'--'}
	formats_gw_eh = {'c':'r','lw':1.05,'ls':'-'}
	
	legend_elements = [Line2D([0],[0],c='k',ls='--',label='Without excitons'), Line2D([0],[0],c='r',ls='-',label='With excitons')]
	legend = {'loc':0,'handles':legend_elements,'prop':{'size':7},'frameon':True,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':1.8}

	fig,ax = plt.subplots(figsize=(3.36,3.36))
	plot_data(ax,data['eng'],data['eps2_noeh'],'default','default','Energy (eV)',r'$\epsilon_2$',formats_gw_noeh,None)
	plot_data(ax,data['eng'],data['eps2_eh'],(2,10),(0,40),'Energy (eV)',r'$\epsilon_2$',formats_gw_eh,legend)

	minorLocator_x = MultipleLocator(0.2)
	minorLocator_y = MultipleLocator(1)
	ax.xaxis.set_minor_locator(minorLocator_x)
	ax.yaxis.set_minor_locator(minorLocator_y)

	plt.tight_layout()
	plt.savefig('BAs_eps2.png', format='png',dpi=600)
	
