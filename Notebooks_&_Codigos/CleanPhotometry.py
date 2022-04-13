#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 10:37:03 2022

The next class, clean the Photometry and ignore empty values in the tables.
This code was developed by Andrés Ramos Padilla (https://github.com/aframosp/AGNView).
"""

import astropy.units as u
import numpy as np
from numpy import unique as uniq
from astropy.table import Table, vstack


class CleanPhotometry:
    """Clean the photometry in CDS and NED tables"""
    def __init__(self, table_cds, table_ned):
        print('Cleaning')
        self.cds_table = table_cds
        self.ned_table = table_ned
        self.bib_codes = np.array([])
        self.join_tables()

    def get_bibcodes(self, bibcode):
        """Get the bibcodes for the tables to analyse them later."""
        bib = np.array(bibcode)
        self.bib_codes = np.unique(np.concatenate((self.bib_codes, bib)))

    def clean_cds_to(self, table):
        """Remove rows with masked (empty) and null (0) values, then average the values."""
        table.remove_rows(np.where(table['sed_eflux'].mask)[0])
        table.remove_rows(np.where(table['sed_eflux'] == 0.0)[0])
        self.get_bibcodes(table['Bibcode'])
        values, counts = uniq(table['sed_flux'], return_counts=True)
        if len(values) > 0:
            avg = np.average(values, weights=counts)
            std = np.sqrt(np.sum(table['sed_eflux']**2))
            return(avg, std)
        return(np.nan, np.nan)

    def clean_ned_to(self, table):
        """Remove rows with masked (empty) and string values, then average the values."""
        table.remove_rows(np.where(table['Flux_Density'].mask)[0])
        table.remove_rows(np.where(table['NED_Uncertainty'] == '')[0])
        table.remove_rows(np.where(table['NED_Uncertainty'] == '+/-...')[0])
        self.get_bibcodes(table['Refcode'])
        values, counts = uniq(table['Flux_Density'], return_counts=True)
        table['NED_Uncertainty'] = [float(j.split('+/-')[-1]) for j in table['NED_Uncertainty']]
        if len(values) > 0:
            avg = np.average(values, weights=counts)
            std = np.sqrt(np.sum(table['NED_Uncertainty']**2))
            return(avg, std)
        return(np.nan, np.nan)

    def filters_ned(self, table_filter):
        """ Use only the information that we need from NED """
        data_tab = []
        band = [np.where(self.ned_table['Observed_Passband'] == filt)[0] for filt in table_filter]
        for bandinx, bandif in enumerate(band):
            avg, std = self.clean_ned_to(self.ned_table[bandif])
            if np.isnan(avg):
                continue
            to_um = self.ned_table[bandif]['Frequency'].to(u.micron,
                                                           equivalencies=u.spectral())
            data_tab.append([table_filter[bandinx], np.mean(to_um), avg*u.Jy, std*u.Jy])
        return(Table(np.array(data_tab),
                     names=['Filter', 'Wave', 'Flux', 'F_er'],
                     dtype=('U32', 'float64', 'float64', 'float64')))

    def filters_cds(self, table_filter):
        """Use only the information that we need from CDS"""
        data_tab = []
        band = [np.where(self.cds_table['sed_filter'] == filt)[0] for filt in table_filter]
        for bandinx, bandif in enumerate(band):
            avg, std = self.clean_cds_to(self.cds_table[bandif])
            if np.isnan(avg):
                continue
            to_um = self.cds_table[bandif]['sed_freq'].to(u.micron,
                                                          equivalencies=u.spectral())
            data_tab.append([table_filter[bandinx], np.mean(to_um), avg*u.Jy, std*u.Jy])
        return(Table(np.array(data_tab),
                     names=['Filter', 'Wave', 'Flux', 'F_er'],
                     dtype=('U32', 'float64', 'float64', 'float64')))

    def join_tables(self):
        """Finally, we join the tables and remove those whose errors are too high"""
        self.final_tab = vstack([self.filters_cds(CDSFilters), self.filters_ned(NEDFilters)])
        self.final_tab.remove_rows(np.where(self.final_tab['F_er']/self.final_tab['Flux'] >= 1))
        
CDSFilters = ['GALEX:FUV', 'GALEX:NUV',
              "SDSS:u'", "SDSS:g'", "SDSS:r'", "SDSS:i'", "SDSS:z'",
              'SDSS:u', 'SDSS:g', 'SDSS:r', 'SDSS:i', 'SDSS:z',
              '2MASS:J', '2MASS:H', '2MASS:Ks',
              ':=3.6um', ':=4.5um', ':=5.8um', ':=8um',
              'WISE:W1', 'WISE:W2', 'WISE:W3', 'WISE:W4',
              'IRAS:12', 'IRAS:25', 'IRAS:60', 'IRAS:100',
              'Spitzer/MIPS:24', 'Spitzer/MIPS:70', 'Spitzer/MIPS:160',
              'Herschel/PACS:70', 'Herschel/PACS:100', 'Herschel/PACS:160',
              'Herschel/SPIRE:250', 'Herschel/SPIRE:350', 'Herschel/SPIRE:500',
              ':=250um', ':=350um', ':=500um',
              ':=1.4GHz', ':=21cm', ':=1.5GHz', ':=20cm', ':=5GHz', ':=6cm']

NEDFilters = ['2-10 keV (XMM)', '0.5-2 keV (XMM)',
              'FUV (GALEX)', 'NUV (GALEX)',
              'u (SDSS) AB', 'g (SDSS) AB', 'r (SDSS) AB', 'i (SDSS) AB', 'z (SDSS) AB',
              'J (2MASS) AB', 'H (2MASS) AB', 'Ks (2MASS) AB',
              'W1 (WISE)', 'W2 (WISE)', 'W3 (WISE)', 'W4 (WISE)',
              '3.6 microns (IRAC)', '4.5 microns (IRAC)',
              '5.8 microns (IRAC)', '8.0 microns (IRAC)',
              '12 microns (IRAS)', '25 microns (IRAS)',
              '60 microns (IRAS)', '100 microns (IRAS)',
              '24 microns (MIPS)', '70 microns (MIPS)', '160 microns (MIPS)',
              '70 microns (PACS)', '100 microns (PACS)', '160 microns (PACS)',
              '250 microns (SPIRE)', '350 microns (SPIRE)', '500 microns (SPIRE)',
              '4.89 GHz (VLA)', '1.46 GHz (VLA)', '1.4GHz']