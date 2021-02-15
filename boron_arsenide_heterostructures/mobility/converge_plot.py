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


def plot_data(ax,xdata,ydata,x_range,y_range,x_label,y_label,formats,legend_data):
   if formats == 'default':
      c='k'
      lw=0.95
      ls='-'
   else:
      c=formats['c']
      lw=formats['lw']
      ls=formats['ls']
      marker=formats['marker']
   
   ax.plot(xdata,ydata,c=c,lw=lw,ls=ls,marker=marker,markersize=3)

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

   electron_hole_combo_top3_0y = True
   use_this_min = 72

   unstrain_data = pd.read_csv('unstrained_convergence/convergence.dat',delim_whitespace=True,header=None,names=['kpt','kpt3','irq','emx','emy','emz','ema','hmx','hmy','hmz','hma','ex','ey','ez','hx','hy','hz','esprd','hsprd'])
   strain_data = pd.read_csv('strained_convergence/convergence.dat',delim_whitespace=True,header=None,names=['kpt','kpt3','irq','emx','emy','emz','ema','hmx','hmy','hmz','hma','ex','ey','ez','hx','hy','hz','esprd','hsprd'])

   formats_kl = {'c':'k','lw':1.05,'ls':'-','marker':'o'}
   formats_kd = {'c':'k','lw':1.05,'ls':'--','marker':'o'}
   formats_kdd = {'c':'k','lw':1.05,'ls':':','marker':'o'}
   formats_rl = {'c':'r','lw':1.05,'ls':'-','marker':'o'}
   formats_rd = {'c':'r','lw':1.05,'ls':'--','marker':'o'}
   formats_rdd = {'c':'r','lw':1.05,'ls':':','marker':'o'}

   e_lit = 1400
   h_lit = 2110

   unstrain_data = unstrain_data.drop([0,1,2])
   strain_data = strain_data.drop([0,1,2])
   e_unstrained = 1341
   e_low = 482
   e_high = 2417
   h_unstrained = 1387
   h_low = 2956
   h_high = 3550
   
   fig = plt.figure(figsize=(6.6,3))
   gs = gridspec.GridSpec(1,2)
   e_ax = fig.add_subplot(gs[0,0])
   h_ax = fig.add_subplot(gs[0,1])

   plot_data(e_ax,1/unstrain_data['kpt3'],unstrain_data['ema'],'default','default',r'$\frac{1}{\mathrm{k(q)\ point}^3}$',r'Mobility (cm$^2$/Vs)',formats_kl,None)
   plot_data(e_ax,1/strain_data['kpt3'],(strain_data['emy']+strain_data['emz'])/2,'default','default',r'$\frac{1}{\mathrm{k(q)\ point}^3}$',r'Mobility (cm$^2$/Vs)',formats_kd,None)
   plot_data(e_ax,1/strain_data['kpt3'],strain_data['emx'],'default','default',r'$\frac{1}{\mathrm{k/q\ point}^3}$',r'$\mu_e$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_kdd,None)

   add_label('0% strain\n'+r'($\mathrm{m_e^*=(2\times m_t+m_l)/3=0.52m_0}$)','k',8,e_ax,(0.25,0.37))
   add_label('1% strain, in-plane\n'+r'($\mathrm{m_e^*=m_t=0.24m_0}$)','k',8,e_ax,(0.42,0.78))
   add_label('1% strain, out-of-plane\n'+r'($\mathrm{m_e^*=m_l=1.09m_0}$)','k',8,e_ax,(0.32,0.042))
   add_label('Liu et al.','b',8,e_ax,(0.17,0.58))

   e_ax.scatter([0],[e_unstrained],c='k',marker='+',s=15) 
   e_ax.scatter([0],[e_low],c='k',marker='+',s=15) 
   e_ax.scatter([0],[e_high],c='k',marker='+',s=15) 
   bottom,top = e_ax.get_ylim()
   e_ax.vlines(x=0,ymin=0,ymax=top,lw=0.8,linestyles='dotted',zorder=-1)
   e_ax.set_ylim((0,top))
   left,right = e_ax.get_xlim()
   e_ax.hlines(y=e_lit,xmin=0,xmax=right,lw=0.8,linestyles='dotted',zorder=-2,color='b')
   e_ax.scatter([0],[e_lit],c='b',marker='+',s=15)
   e_ax.annotate('',xy=(-0.2/(120.0**3),e_high),xytext=(-0.2/(120.0**3),e_unstrained),arrowprops=dict(arrowstyle='-|>',facecolor='k'))
   e_ax.annotate('+{num:.1f}%'.format(num=(e_high-e_unstrained)/e_unstrained*100),xy=(-0.2/(120.0**3)*2.7,(e_high+e_unstrained)/2*1.07),size=8,color='k',xycoords='data',rotation=90)
   e_ax.set_xlim((left*3,right))

   plot_data(h_ax,1/unstrain_data['kpt3'],unstrain_data['hma'],'default','default',r'$\frac{1}{\mathrm{k(q)\ point}^3}$',r'Mobility (cm$^2$/Vs)',formats_rl,None)
   plot_data(h_ax,1/strain_data['kpt3'],(strain_data['hmy']+strain_data['hmz'])/2,'default','default',r'$\frac{1}{\mathrm{k(q)\ point}^3}$',r'Mobility (cm$^2$/Vs)',formats_rd,None)
   plot_data(h_ax,1/strain_data['kpt3'],strain_data['hmx'],'default','default',r'$\frac{1}{\mathrm{k/q\ point}^3}$',r'$\mu_h$'+' ['+r'$\mathrm{cm^{2}V^{-1}s^{-1}}$'+']',formats_rdd,None)

   add_label('0% strain\n(w/o LS, 3 bands)','r',8,h_ax,(0.45,0.22))
   add_label('1% strain, in-plane (w/o LS, 1 band)','r',8,h_ax,(0.22,0.855))
   add_label('1% strain, out-of-plane\n(w/o LS, 1 band)','r',8,h_ax,(0.35,0.66))
   add_label('Liu et al.\n(w/ LS, 2 bands)','b',8,h_ax,(0.17,0.45))

   h_ax.scatter([0],[h_unstrained],c='r',marker='+',s=15)
   h_ax.scatter([0],[h_low],c='r',marker='+',s=15) 
   h_ax.scatter([0],[h_high],c='r',marker='+',s=15) 
   bottom,top = h_ax.get_ylim()
   h_ax.vlines(x=0,ymin=0,ymax=top,lw=0.8,linestyles='dotted',zorder=-1)
   h_ax.set_ylim((0,top))
   left,right = h_ax.get_xlim()      
   h_ax.hlines(y=h_lit,xmin=0,xmax=right,lw=0.8,linestyles='dotted',zorder=-2,color='b')
   h_ax.scatter([0],[h_lit],c='b',marker='+',s=15)
   h_ax.annotate('',xy=(-0.2/(120.0**3),h_high),xytext=(-0.2/(120.0**3),h_lit),arrowprops=dict(arrowstyle='-|>',facecolor='k'))
   h_ax.annotate('+{num:.1f}%'.format(num=(h_high-h_lit)/h_lit*100),xy=(-0.2/(120.0**3)*2.7,(h_high+h_lit)/2*1.07),size=8,color='k',xycoords='data',rotation=90)
   h_ax.set_xlim((left*3,right))

   xticks = list(reversed(1/unstrain_data['kpt3'].values))
   labels = [r'$\frac{1}{'+str(i)+r'^3}$' for i in list(reversed(unstrain_data['kpt'].values))]
   xticks = [0] + xticks
   labels = [0] + labels

   for ax in [e_ax, h_ax]:
      ax.set_xticks(xticks)
      ax.set_xticklabels(labels)
      ax.xaxis.label.set_size(12)
      ax.yaxis.label.set_size(10)
      ax.tick_params(axis='both', which='major', labelsize=9)

   add_label('(a)','k',10,e_ax,(-0.295,0.95))
   add_label('(b)','k',10,h_ax,(-0.295,0.95))

   gs.tight_layout(fig)
   plt.savefig('BAs_FigS5_'+date+'.pdf', format='pdf',dpi=600)
