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

font = {'family':'Helvetica','size': 10}
plt.rc('font', **font)

def add_label(label,color,size,ax,loc):
   ax.annotate(label,xy=loc,size=size,color=color,xycoords='axes fraction')

def fx(xdata,an):
   return [an[0]*x**an[1] for x in xdata]

def plot_data(ax,xdata,ydata,x_range,y_range,x_label,y_label,formats,legend_data):
   if formats == 'default':
      c='k'
      lw=0.95
      ls='-'
      marker = 'None'
   else:
      c=formats['c']
      lw=formats['lw']
      ls=formats['ls']
      marker=formats['marker']

   if 'alpha' in formats.keys():
      ax.plot(xdata,ydata,c=c,lw=lw,ls=ls,alpha=formats['alpha'],marker=marker,markersize=3)
   else:
      ax.plot(xdata,ydata,c=c,lw=lw,ls=ls,marker=marker,markersize=3)

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

   date = 'plot_date'

   unstrain_data = pd.read_csv('converged_unstrained/temp_72_data.dat',delim_whitespace=True,header=None,names=['t','emx','emy','emz','ema','hmx','hmy','hmz','hma','emf','hmf','esprd','hsprd'])
   strain_data = pd.read_csv('converged_strained/temp_72_data_strain.dat',delim_whitespace=True,header=None,names=['t','emx','emy','emz','ema','hmx','hmy','hmz','hma','emf','hmf','esprd','hsprd'])
   liu_elec = pd.read_csv('liu_temp_elecs.csv',delim_whitespace=True,header=None,names=['t','ema'])
   liu_hole = pd.read_csv('liu_temp_holes.csv',delim_whitespace=True,header=None,names=['t','hma'])
   
   formats_kc = {'c':'k','lw':1.05,'ls':'None','marker':'o'}
   formats_ks = {'c':'k','lw':1.05,'ls':'None','marker':'s'}
   formats_kt = {'c':'k','lw':1.05,'ls':'None','marker':'^'}
   formats_rc = {'c':'r','lw':1.05,'ls':'None','marker':'o'}
   formats_rs = {'c':'r','lw':1.05,'ls':'None','marker':'s'}
   formats_rt = {'c':'r','lw':1.05,'ls':'None','marker':'^'}
   
   formats_liu = {'c':'xkcd:bright blue','lw':1,'ls':'-','marker':None}
   formats_trans = {'c':'k','lw':0.7,'ls':'-','marker':None,'alpha':0.6}

   formats_kl = {'c':'k','lw':1.05,'ls':'-','marker':'o'}
   formats_kd = {'c':'k','lw':1.05,'ls':'--','marker':'s'}
   formats_kdd = {'c':'k','lw':1.05,'ls':':','marker':'^'}
   formats_rl = {'c':'r','lw':1.05,'ls':'-','marker':'o'}
   formats_rd = {'c':'r','lw':1.05,'ls':'--','marker':'s'}
   formats_rdd = {'c':'r','lw':1.05,'ls':':','marker':'^'}

   fig = plt.figure(figsize=(6.6,6))
   gs = gridspec.GridSpec(2,2)
   e_ax = fig.add_subplot(gs[0,0])
   h_ax = fig.add_subplot(gs[0,1])
   mlx = MultipleLocator(50)
   mly = MultipleLocator(500)
   
   e_legend_elements = [Line2D([0],[0],c='k',ls='-',marker='o',markersize=4,label='No strain\n'+r'($\mathrm{m_e^*=(2m_t+m_l)/3}$'+'\n       '+r'$\mathrm{=0.52m_0}$)'),Line2D([0],[0],c='k',ls='--',marker='s',markersize=4,label='Strain, in-plane\n'+r'($\mathrm{m_e^*=m_t=0.24m_0}$)'),
   Line2D([0],[0],c='k',ls=':',marker='^',markersize=4,label='Strain, out-of-plane\n'+r'($\mathrm{m_e^*=m_l=1.09m_0}$)'),Line2D([0],[0],c='xkcd:bright blue',ls='-',label='No strain, Liu et al.')]
   h_legend_elements = [Line2D([0],[0],c='r',ls='-',marker='o',markersize=4,label='No strain\n'+'(w/o LS, 3 bands)'),Line2D([0],[0],c='r',ls='--',marker='s',markersize=4,label='Strain, in-plane\n'+'(w/o LS, 1 band)'),
   Line2D([0],[0],c='r',ls=':',marker='^',markersize=4,label='Strain, out-of-plane\n'+'(w/o LS, 1 band)'),Line2D([0],[0],c='xkcd:bright blue',ls='-',label='No strain, Liu et al.\n'+'(w/ LS, 2 bands)')]

   e_legend = {'loc':0,'handles':e_legend_elements,'prop':{'size':7.5},'frameon':False,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.85,'labelspacing':0.75} 
   h_legend = {'loc':0,'handles':h_legend_elements,'prop':{'size':7.5},'frameon':False,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2.85,'labelspacing':0.75} 

   plot_data(e_ax,unstrain_data['t'],unstrain_data['ema'],'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_kl,None)
   plot_data(e_ax,strain_data['t'],(strain_data['emy']+strain_data['emz'])/2,'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_kd,None)
   plot_data(e_ax,strain_data['t'],strain_data['emx'],'default','default',r'T (K)',r'$\mathrm{\mu_e}$ (cm$^2$/Vs)',formats_kdd,e_legend)
   plot_data(e_ax,liu_elec['t'],liu_elec['ema'],'default','default',r'$T$'+' [K]',r'$\mu_e$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_liu,None)

   plot_data(h_ax,liu_hole['t'],liu_hole['hma'],'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_liu,None)
   plot_data(h_ax,unstrain_data['t'],unstrain_data['hma'],'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_rl,None)
   plot_data(h_ax,strain_data[strain_data['t']>151]['t'],(strain_data[strain_data['t']>151]['hmy']+strain_data[strain_data['t']>151]['hmz'])/2,'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_rd,None)
   plot_data(h_ax,strain_data[strain_data['t']>151]['t'],strain_data[strain_data['t']>151]['hmx'],'default','default',r'$T$'+' [K]',r'$\mu_h$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_rdd,h_legend)

   for ax in [e_ax, h_ax]:
      ax.xaxis.label.set_size(10)
      ax.yaxis.label.set_size(10)
      ax.tick_params(axis='both', which='major', labelsize=9)
      ax.xaxis.set_minor_locator(mlx)
      ax.yaxis.set_minor_locator(mly)

   add_label('(a)','k',10,e_ax,(-0.295,0.95))
   add_label('(b)','k',10,h_ax,(-0.295,0.95))

   e_ax = fig.add_subplot(gs[1,0])
   h_ax = fig.add_subplot(gs[1,1])

   xticks = [100,200,300,400,500]
   labels = [str(i) for i in xticks]
   mlx = MultipleLocator(50)

   for ax in [e_ax, h_ax]:
      ax.set_xscale('log')
      ax.set_yscale('log')
      ax.set_xticks(xticks)
      ax.set_xticklabels(labels)
      ax.xaxis.label.set_size(10)
      ax.yaxis.label.set_size(10)
      ax.xaxis.set_minor_locator(mlx)
      ax.tick_params(axis='both', which='major', labelsize=9)
      ax.tick_params(axis='x', which='minor', labelbottom=False)

   e_legend_elements = [Line2D([0],[0],c='k',ls='None',marker='o',markersize=3,label='No strain'),Line2D([0],[0],c='k',ls='None',marker='s',markersize=3,label='Strain, in-plane'),
   Line2D([0],[0],c='k',ls='None',marker='^',markersize=3,label='Strain, out-of-plane'),Line2D([0],[0],c='xkcd:bright blue',ls='-',label='No strain, Liu et al.')]
   e_legend = {'loc':3,'handles':e_legend_elements,'prop':{'size':7.5},'frameon':False,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2} 
   h_legend_elements = [Line2D([0],[0],c='r',ls='None',marker='o',markersize=3,label='No strain'),Line2D([0],[0],c='r',ls='None',marker='s',markersize=3,label='Strain, in-plane'),
   Line2D([0],[0],c='r',ls='None',marker='^',markersize=3,label='Strain, out-of-plane'),Line2D([0],[0],c='xkcd:bright blue',ls='-',label='No strain, Liu et al.')]
   h_legend = {'loc':3,'handles':h_legend_elements,'prop':{'size':7.5},'frameon':False,'shadow':False,'framealpha':1,'edgecolor':'k','handlelength':2} 

   plot_data(e_ax,unstrain_data['t'],unstrain_data['ema'],'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_kc,None)
   plot_data(e_ax,strain_data['t'],(strain_data['emy']+strain_data['emz'])/2,'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_ks,None)
   plot_data(e_ax,strain_data['t'],strain_data['emx'],'default','default',r'T (K)',r'$\mathrm{\mu_e}$ (cm$^2$/Vs)',formats_kt,e_legend)

   ys = e_ax.get_ylim()
   xs = e_ax.get_xlim()

   plot_data(e_ax,liu_elec['t'],liu_elec['ema'],'default','default',r'T (K)',r'$\mathrm{\mu_e}$ (cm$^2$/Vs)',formats_liu,e_legend)

   u_e_l = [4.19567e7,-1.77703]
   u_e_h = [2.93531e8,-2.14745]
   ip_e_l = [1.24024e7,-1.47533]
   ip_e_h = [8.23208e7,-1.82769]
   op_e_l = [3.19376e6,-1.52366]
   op_e_h = [1.08321e7,-1.75199]
   temp_range_low = [50,225]
   temp_range_high = [180,600]
   for params in [u_e_l,ip_e_l,op_e_l]:
      plot_data(e_ax,temp_range_low,fx(temp_range_low,params),'default','default',r'T (K)',r'$\mathrm{\mu_e}$ (cm$^2$/Vs)',formats_trans,None)
   for params in [u_e_h,ip_e_h,op_e_h]:
      plot_data(e_ax,temp_range_high,fx(temp_range_high,params),'default','default',r'$T$'+' [K]',r'$\mu_e$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_trans,None)

   e_ax.set_ylim(ys)
   e_ax.set_xlim(xs)
   
   add_label(r'${\sim}\mathrm{T}^{-1.52}$','k',8,e_ax,(0.06,0.44))
   add_label(r'${\sim}\mathrm{T}^{-1.75}$','k',8,e_ax,(0.62,0.08))
   add_label(r'${\sim}\mathrm{T}^{-1.77}$','k',8,e_ax,(0.05,0.72))
   add_label(r'${\sim}\mathrm{T}^{-2.15}$','k',8,e_ax,(0.43,0.45))
   add_label(r'${\sim}\mathrm{T}^{-1.48}$','k',8,e_ax,(0.2,0.9))
   add_label(r'${\sim}\mathrm{T}^{-1.82}$','k',8,e_ax,(0.74,0.57))
   add_label(r'${\sim}\mathrm{T}^{-2.34}$','xkcd:bright blue',8,e_ax,(0.54,0.34))

   plot_data(h_ax,unstrain_data['t'],unstrain_data['hma'],'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_rc,None)
   plot_data(h_ax,strain_data[strain_data['t']>151]['t'],(strain_data[strain_data['t']>151]['hmy']+strain_data[strain_data['t']>151]['hmz'])/2,'default','default',r'T (K)',r'Mobility (cm$^2$/Vs)',formats_rs,None)
   plot_data(h_ax,strain_data[strain_data['t']>151]['t'],strain_data[strain_data['t']>151]['hmx'],'default','default',r'T (K)',r'$\mathrm{\mu_h}$ (cm$^2$/Vs)',formats_rt,h_legend)

   ys = h_ax.get_ylim()
   xs = h_ax.get_xlim()

   plot_data(h_ax,liu_hole['t'],liu_hole['hma'],'default','default',r'T (K)',r'$\mathrm{\mu_e}$ (cm$^2$/Vs)',formats_liu,h_legend)


   u_h_l = [3.77459e6,-1.3112]
   u_h_h = [1.16774e9,-2.40649]
   ip_h_h = [1.11914e11,-3.03446]
   op_h_h = [1.23933e11,-3.074]

   plot_data(h_ax,temp_range_low,fx(temp_range_low,u_h_l),'default','default',r'T (K)',r'$\mathrm{\mu_h}$ (cm$^2$/Vs)',formats_trans,None)
   for params in [u_h_h,ip_h_h,op_h_h]:
      plot_data(h_ax,temp_range_high,fx(temp_range_high,params),'default','default',r'$T$'+' [K]',r'$\mu_h$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_trans,None)

   h_ax.set_ylim(ys)
   h_ax.set_xlim(xs)

   add_label(r'${\sim}\mathrm{T}^{-1.31}$','k',8,h_ax,(0.06,0.7))
   add_label(r'${\sim}\mathrm{T}^{-2.41}$','k',8,h_ax,(0.5,0.28))
   add_label(r'${\sim}\mathrm{T}^{-3.07}$','k',8,h_ax,(0.28,0.85))
   add_label(r'${\sim}\mathrm{T}^{-3.03}$','k',8,h_ax,(0.72,0.65))
   add_label(r'${\sim}\mathrm{T}^{-2.81}$','xkcd:bright blue',8,h_ax,(0.34,0.75))

   add_label('(c)','k',10,e_ax,(-0.295,0.95))
   add_label('(d)','k',10,h_ax,(-0.295,0.95))

   gs.tight_layout(fig)
   plt.savefig('BAs_Fig4_'+date+'.pdf', format='pdf',dpi=600)

   
