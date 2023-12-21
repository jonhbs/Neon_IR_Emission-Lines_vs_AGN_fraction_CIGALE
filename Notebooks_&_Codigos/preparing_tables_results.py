#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 02:36:36 2023

@author: Jonhatan Bernal Salinas
jobernals@unal.edu.co
"""
import numpy as np
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')
from astropy.table import Table, join

path_FINAL = '../Datos/Final/'
path_CIGALE_out = '../Datos/Work/CIGALE_Output/'

#Reading the data with Ne line ratios, results from CIGALE, and phot data:
Sy_Ne = Table.read(path_FINAL+'Ne_ratios_-_CIGALE_results.tbl', format='ascii')
results = Table.read(path_CIGALE_out+'results0a90_All.fits',format='fits')
Phot = Table.read(path_FINAL+'CIGPhotSy_EnergyBal_All.tbl', format='ascii')

#Rename id column in Phot, to join with other tables:
Phot.rename_column('id', 'Main_id') 
#Adding information about number of phot points for each galaxy:
Phot['CountPoints'] = [sum(~np.isnan(list(j)))/
                       2 for j in Phot[Phot.colnames[2:]].as_array()]
#Join Phot table and Ne information
Ne_Sy = join(Phot, Sy_Ne)
print('Total number of galaxies with Ne ratios data: ',len(Ne_Sy))

#-----------------------------------------------------------------------------
#Convertimg the luminosities form Watts (W) to Solar Luminosities
for cigpar in ['bayes.agn.disk_luminosity']:
    Ne_Sy[cigpar] = (Ne_Sy[cigpar]*u.W).to(u.solLum)
    Ne_Sy[cigpar+'_err'] = (Ne_Sy[cigpar+'_err']*u.W).to(u.solLum)
for cigpar2 in ['bayes.agn.luminosity']:
    Ne_Sy[cigpar2] = (Ne_Sy[cigpar2]*u.W).to(u.solLum)
    Ne_Sy[cigpar2+'_err'] = (Ne_Sy[cigpar2+'_err']*u.W).to(u.solLum)
#-----------------------------------------------------------------------------
#Logaritmic values and Error Propagation
#For Neon lines:
Ne_Sy[r'log_[NeV]*'] = np.log10(Ne_Sy[r'[NeV]*'])
Ne_Sy[r'log_[NeII]*'] = np.log10(Ne_Sy[r'[NeII]*'])
Ne_Sy[r'log_[NeV]$^+$'] = np.log10(Ne_Sy[r'[NeV]$^+$'])
Ne_Sy[r'log_[NeIII]$^+$'] = np.log10(Ne_Sy[r'[NeIII]$^+$'])
Ne_Sy[r'log_[NeIII]$^-$'] = np.log10(Ne_Sy[r'[NeIII]$^-$'])
Ne_Sy[r'log_[NeII]$^-$'] = np.log10(Ne_Sy[r'[NeII]$^-$'])

Ne_Sy[r'log_[NeV]_err*'] = Ne_Sy[r'[NeV]_err*']/(Ne_Sy[r'[NeV]*']*np.log(10))
Ne_Sy[r'log_[NeII]_err*'] = Ne_Sy[r'[NeII]_err*']/\
    (Ne_Sy[r'[NeII]*']*np.log(10))
Ne_Sy[r'log_[NeV]_err$^+$'] = Ne_Sy[r'[NeV]_err$^+$']/\
    (Ne_Sy[r'[NeV]$^+$']*np.log(10))
Ne_Sy[r'log_[NeIII]_err$^+$'] = Ne_Sy[r'[NeIII]_err$^+$']/\
    (Ne_Sy[r'[NeIII]$^+$']*np.log(10))
Ne_Sy[r'log_[NeIII]_err$^-$'] = Ne_Sy[r'[NeIII]_err$^-$']/\
    (Ne_Sy[r'[NeIII]$^-$']*np.log(10))
Ne_Sy[r'log_[NeII]_err$^-$'] = Ne_Sy[r'[NeII]_err$^-$']/\
    (Ne_Sy[r'[NeII]$^-$']*np.log(10))

#For Neon line ratios:
Ne_Sy['log_[NeV]/[NeII]'] = np.log10(Ne_Sy[r'[NeV]/[NeII]*'])
Ne_Sy['log_[NeV]/[NeII]_err'] = Ne_Sy[r'[NeV]/[NeII]_err*']/\
    (Ne_Sy[r'[NeV]/[NeII]*']*np.log(10))
Ne_Sy['log_[NeV]/[NeIII]'] = np.log10(Ne_Sy[r'[NeV]/[NeIII]$^+$'])
Ne_Sy['log_[NeV]/[NeIII]_err'] = Ne_Sy[r'[NeV]/[NeIII]_err$^+$']/\
    (Ne_Sy[r'[NeV]/[NeIII]$^+$']*np.log(10))
Ne_Sy['log_[NeIII]/[NeII]'] = np.log10(Ne_Sy[r'[NeIII]/[NeII]$^-$'])
Ne_Sy['log_[NeIII]/[NeII]_err'] = Ne_Sy[r'[NeIII]/[NeII]_err$^-$']/\
    (Ne_Sy[r'[NeIII]/[NeII]$^-$']*np.log(10))
    
#For physical properties estimated with CIGALE:
Ne_Sy['log_bayes.sfh.sfr'] = np.log10(Ne_Sy['bayes.sfh.sfr'])
Ne_Sy['log_bayes.sfh.sfr_err'] = Ne_Sy['bayes.sfh.sfr_err']/\
    (Ne_Sy['bayes.sfh.sfr']*np.log(10))
Ne_Sy['log_bayes.stellar.m_star'] = np.log10(Ne_Sy['bayes.stellar.m_star'])
Ne_Sy['log_bayes.stellar.m_star_err'] = Ne_Sy['bayes.stellar.m_star_err']/\
    (Ne_Sy['bayes.stellar.m_star']*np.log(10))
Ne_Sy['log_bayes.agn.disk_luminosity'] = np.log10(
    Ne_Sy['bayes.agn.disk_luminosity'])
Ne_Sy['log_bayes.agn.disk_luminosity_err'] = \
    Ne_Sy['bayes.agn.disk_luminosity_err']/\
        (Ne_Sy['bayes.agn.disk_luminosity']*np.log(10))
Ne_Sy['log_bayes.agn.luminosity'] = np.log10(Ne_Sy['bayes.agn.luminosity'])
Ne_Sy['log_bayes.agn.luminosity_err'] = Ne_Sy['bayes.agn.luminosity_err']/\
    (Ne_Sy['bayes.agn.luminosity']*np.log(10))
#-----------------------------------------------------------------------------    
#Saving the table
Ne_Sy.write(path_FINAL+'final_results.tbl',format='ascii',overwrite=True)
Ne_Sy.write(path_FINAL+'final_results.fits',format='fits',overwrite=True)