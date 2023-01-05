#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 11:12:11 2023

@author: Jonhatan Bernal Salinas
"""

from astropy.table import Table
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd



path_CIGALE_phot = '../Datos/Work/CIGALE_InputPhot/'
path_CIGALE_out = '../Datos/Work/CIGALE_Output/'
path_graphs = '../Datos/Final/Graphics/'
path_FINAL = '../Datos/Final/'

Ne_Sy_bib = Table.read(path_FINAL+'Ne_ratios_-_CIGALE_results_bib.tbl', format='ascii')

Ne_Sy1_bib = Ne_Sy_bib[Ne_Sy_bib['otype']=='Sy1']
Ne_Sy2_bib = Ne_Sy_bib[Ne_Sy_bib['otype']=='Sy2']
Ne_SyG_bib = Ne_Sy_bib[Ne_Sy_bib['otype']=='SyG']

# Graphics of Ne Line Ratios Histograms and Ne Line Ratios vs fracAGN but now by bibcode in each galaxy
### [NeV]/[NeII]

x_NeV_NeII_Sy1_bib = np.log10(Ne_Sy1_bib['[NeV]/[NeII]'])
x_NeV_NeII_Sy2_bib = np.log10(Ne_Sy2_bib['[NeV]/[NeII]'])
x_NeV_NeII_SyG_bib = np.log10(Ne_SyG_bib['[NeV]/[NeII]'])
bins_NeV_NeII_Sy1_bib = np.linspace(-2.5,1,50)
bins_NeV_NeII_Sy2_bib = np.linspace(-2.5,1,50)
bins_NeV_NeII_SyG_bib = np.linspace(-2.5,1,50)
plt.figure(figsize=(12,7))
plt.hist(x_NeV_NeII_Sy1_bib, bins=bins_NeV_NeII_Sy1_bib, histtype='step', label='Sy1')
plt.hist(x_NeV_NeII_Sy2_bib, bins=bins_NeV_NeII_Sy2_bib, histtype='step', label='Sy2')
plt.hist(x_NeV_NeII_SyG_bib, bins=bins_NeV_NeII_SyG_bib, histtype='step', label='SyG')
plt.xlabel('log([NeV]/[NeII])',fontsize=14)
plt.ylabel('Number of galaxies',fontsize=14)
plt.title('[NeV]/[NeII] Histogram (by bibcode)', fontsize=18)
plt.legend()
plt.savefig(path_graphs+'hist_NeV_NeII_Sy_bib.jpg')

plt.figure(figsize=(12,7))
plt.scatter(Ne_Sy1_bib['[NeV]/[NeII]'],Ne_Sy1_bib['bayes.agn.fracAGN'],200,color='blue',marker = '*',label='Sy1')
plt.errorbar(Ne_Sy1_bib['[NeV]/[NeII]'], Ne_Sy1_bib['bayes.agn.fracAGN'],\
             Ne_Sy1_bib['bayes.agn.fracAGN_err'], Ne_Sy1_bib['[NeV]/[NeII]_err'], fmt='k.')
plt.scatter(Ne_Sy2_bib['[NeV]/[NeII]'], Ne_Sy2_bib['bayes.agn.fracAGN'],120, color='red', marker = 'v', label='Sy2')
plt.errorbar(Ne_Sy2_bib['[NeV]/[NeII]'], Ne_Sy2_bib['bayes.agn.fracAGN'],\
             Ne_Sy2_bib['bayes.agn.fracAGN_err'], Ne_Sy2_bib['[NeV]/[NeII]_err'], fmt='k.')
plt.scatter(Ne_SyG_bib['[NeV]/[NeII]'], Ne_SyG_bib['bayes.agn.fracAGN'],200, color='black', marker = '.', label='SyG')
plt.errorbar(Ne_SyG_bib['[NeV]/[NeII]'], Ne_SyG_bib['bayes.agn.fracAGN'],\
             Ne_SyG_bib['bayes.agn.fracAGN_err'], Ne_SyG_bib['[NeV]/[NeII]_err'], fmt='k.')
plt.xlabel(r'[NeV]/[NeII]', fontsize=14)
plt.xscale('log')
plt.ylabel(r'fracAGN', fontsize=14)
plt.ylim(0.0,1.0)
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.title(r'fracAGN vs [NeV]/[NeII] (by bibcode)', fontsize=18)
plt.legend()
plt.savefig(path_graphs+'NeV_NeII_vs_fracAGN_Sy_bib.jpg')

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(16,6),sharey=True)
fig.suptitle(r'fracAGN vs [NeV]/[NeII] (by bibcode)', fontsize=18)
ax1.scatter(Ne_Sy1_bib['[NeV]/[NeII]'],Ne_Sy1_bib['bayes.agn.fracAGN'],200,color='blue',marker = '*')
ax1.errorbar(Ne_Sy1_bib['[NeV]/[NeII]'], Ne_Sy1_bib['bayes.agn.fracAGN'],\
             Ne_Sy1_bib['bayes.agn.fracAGN_err'], Ne_Sy1_bib['[NeV]/[NeII]_err'], fmt='k.')
ax1.set_xscale("log")
ax1.set_xlabel(r'[NeV]/[NeII]', fontsize=14)
ax1.set_ylabel(r'fracAGN', fontsize=14)
ax1.set_title('Seyfert 1', fontsize=14)
ax2.scatter(Ne_Sy2_bib['[NeV]/[NeII]'], Ne_Sy2_bib['bayes.agn.fracAGN'],120, color='red', marker = 'v')
ax2.errorbar(Ne_Sy2_bib['[NeV]/[NeII]'], Ne_Sy2_bib['bayes.agn.fracAGN'],\
             Ne_Sy2_bib['bayes.agn.fracAGN_err'], Ne_Sy2_bib['[NeV]/[NeII]_err'], fmt='k.')
ax2.set_xscale("log")
ax2.set_xlabel(r'[NeV]/[NeII]', fontsize=14)
ax2.set_title('Seyfert 2', fontsize=14)
ax1.grid(color='k', linestyle='--', linewidth=0.4)
ax2.grid(color='k', linestyle='--', linewidth=0.4)
plt.subplots_adjust(wspace=0) #Space between subplots
plt.savefig(path_graphs+'NeV_NeII_vs_fracAGN_Sy1_vs_Sy2_bib.jpg')

##Correlations:
    
#Seyfert 1
log_NeVNeII_Sy1_bib = np.log10(Ne_Sy1_bib['[NeV]/[NeII]'])
log_NeVNeIII_Sy1_bib = np.log10(Ne_Sy1_bib['[NeV]/[NeIII]'])
log_NeIIINeII_Sy1_bib = np.log10(Ne_Sy1_bib['[NeIII]/[NeII]'])
log_fracAGN_Sy1_bib = np.log10(Ne_Sy1_bib['bayes.agn.fracAGN'])
log_sfr_Sy1_bib = np.log10(Ne_Sy1_bib['bayes.sfh.sfr'])
log_sm_Sy1_bib = np.log10(Ne_Sy1_bib['bayes.stellar.m_star'])

Ne_Sy1_bib['log_[NeV]/[NeII]'] = log_NeVNeII_Sy1_bib
Ne_Sy1_bib['log_[NeV]/[NeIII]'] = log_NeVNeIII_Sy1_bib
Ne_Sy1_bib['log_[NeIII]/[NeII]'] = log_NeIIINeII_Sy1_bib
Ne_Sy1_bib['log_bayes.agn.fracAGN'] = log_fracAGN_Sy1_bib
Ne_Sy1_bib['log_bayes.sfh.sfr'] = log_sfr_Sy1_bib
Ne_Sy1_bib['log_bayes.stellar.m_star'] = log_sm_Sy1_bib

#Seyfert 2
log_NeVNeII_Sy2_bib = np.log10(Ne_Sy2_bib['[NeV]/[NeII]'])
log_NeVNeIII_Sy2_bib = np.log10(Ne_Sy2_bib['[NeV]/[NeIII]'])
log_NeIIINeII_Sy2_bib = np.log10(Ne_Sy2_bib['[NeIII]/[NeII]'])
log_fracAGN_Sy2_bib = np.log10(Ne_Sy2_bib['bayes.agn.fracAGN'])
log_sfr_Sy2_bib = np.log10(Ne_Sy2_bib['bayes.sfh.sfr'])
log_sm_Sy2_bib = np.log10(Ne_Sy2_bib['bayes.stellar.m_star'])

Ne_Sy2_bib['log_[NeV]/[NeII]'] = log_NeVNeII_Sy2_bib
Ne_Sy2_bib['log_[NeV]/[NeIII]'] = log_NeVNeIII_Sy2_bib
Ne_Sy2_bib['log_[NeIII]/[NeII]'] = log_NeIIINeII_Sy2_bib
Ne_Sy2_bib['log_bayes.agn.fracAGN'] = log_fracAGN_Sy2_bib
Ne_Sy2_bib['log_bayes.sfh.sfr'] = log_sfr_Sy2_bib
Ne_Sy2_bib['log_bayes.stellar.m_star'] = log_sm_Sy2_bib

df_Ne_Sy1_bib = Ne_Sy1_bib.to_pandas()
df_Ne_Sy2_bib = Ne_Sy2_bib.to_pandas()

NeVNeII_corr_Sy1_bib = df_Ne_Sy1_bib.corr()['log_[NeV]/[NeII]']
print(NeVNeII_corr_Sy1_bib,'\n')

NeVNeII_corr_Sy2_bib = df_Ne_Sy2_bib.corr()['log_[NeV]/[NeII]']
print(NeVNeII_corr_Sy2_bib)