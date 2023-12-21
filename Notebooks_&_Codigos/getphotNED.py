#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 15:06:55 2022

@author: Jonhatan Bernal Salinas - jobernals@unal.edu.co
"""

from astroquery.ipac.ned import Ned
from astropy.table import QTable
import astropy.units as u

def photNED(name,oids,coord):
    
    '''
    This function is to get the photometry tables
    from NED. The input are the id, others ids 
    and coordinates of the object.
    '''
    
    phot_ned = QTable()
    id_used = []
    
    try:
        phot_ned = Ned.get_table(name, table='photometry') #We get the Photometry Table using first SIMBAD names
        id_used = name
        #print('Using main SMB id\n')
    except:
        other_ids = oids.split('|')
        #other_ids = Simbad.query_objectids(name) #We get other IDs for the object from SIMBAD using astroquery
        for other_id in other_ids:
            try:
                phot_ned = Ned.get_table(other_id, table='photometry') #Photometry Table using other SIMBAD names
                id_used = other_id
                #print('Trying other id')
                break
            except:
                pass
    if id_used == []:
        try:
            region = Ned.query_region(coord,radius= 5 * u.arcsec)
            id_ned = region["Object Name"]
            phot_ned = Ned.get_table(id_ned, table='photometry')
            id_used = id_ned
            #print('Using region\n')
        except:
            pass

    return phot_ned, id_used