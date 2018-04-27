#!/usr/bin/env python
'''
********************************************************************
Basic GIC calculation method following Lehtinen and Pirjola 1985.
Adapted from ComputeGIC.m from C. Beggan's (BGS) adaptation of
K. Turnbull's Fortran code.

The script requires no input right now and should be run as:
    $ python GIC_Model_Horton.py

NOTES:
-   Due to differences in distance calculations and Python Numpy
    methods, final values are off the exact Horton values by
    a few percent in most cases.
-   Lehtinen-Pirjola method needs more null points added in to 
    equate to nodal admission matrix method, hence the representation
    of lines between HV, LV, ground and switch nodes.
-   Note that in this case the transformer resistance is explicitly 
    given as a connection between nodes.

Created by R Bailey (ZAMG, Austria) on 2015-08-03.
#********************************************************************
'''

import os
import sys
import numpy as np
from scipy import interpolate
from math import radians, tan, atan, atan2, cos, sin, acos, asin
from math import sqrt, pi, log
try:
    import IPython
except:
    pass

#####################################################################
#                       FUNCTIONS                                   #
#####################################################################

def grc_azimuth(lonlat1, lonlat2):
    """
    Function to compute the geographic distance on a prolate ellipsiod such
    as a planet or moon. This computes the distance accurately - not that it
    makes much difference in most cases, as the location of the points of interest is
    generally relatively poorly known.
    
    This function is based on the formula from the Wikipedia page, verified
    against the Geoscience Australia website calculator
    
    Author: Ciaran Beggan
    Rewritten from matlab into python by R Bailey, ZAMG.
    
    Returns azimuth between two points.
    """

    a, b = 6378.137, 6356.752
    f = (a-b)/a

    u1 = atan((1.-f)*tan(pi/180.*(lonlat1[1])))
    u2 = atan((1.-f)*tan(pi/180.*(lonlat2[1])))

    L = pi/180.*(lonlat2[0] - lonlat1[0])
    Lambda, converge, iters = L, False, 0

    while not converge and iters < 20:
        sinsig = sqrt((cos(u2)*sin(Lambda))**2. + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(Lambda))**2.)
        cossig = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(Lambda)
        sig = atan2(sinsig, cossig)
    
        sinalpha = (cos(u1)*cos(u2)*sin(Lambda))/sinsig
        cossqalpha = 1. - sinalpha**2.
        cos2sigm = cossig - (2.*sin(u1)*sin(u2))/cossqalpha

        C = (f/16.) * cossqalpha*(4. + f*(4.-3.*cossqalpha))
    
        calclambda = L + (1.-C)*f*sinalpha*(sig + C*sinalpha*(cos2sigm + C*cossig*(-1. + 2.*cos2sigm) ))
    
        if (abs(Lambda - calclambda) < 10.**(-12.)):
            converge = True
            Lambda = calclambda
        else:
            iters = iters + 1
            Lambda = calclambda

    usq = cossqalpha * ((a**2. - b**2.)/b**2.)
    A = 1. + usq/16384. * (4096. + usq*(-768. + usq*(320. - 175.*usq)))
    B = usq/1024.* (256. + usq*(-128. + usq*(74. - 47.*usq)))
    delsig = B * sinsig * (cos2sigm + 0.25 * B *(cossig *(-1. + 2.*cos2sigm) -(1./6.)*B * cos2sigm*(-3. + 4.*sinalpha**2.)*(-3.+4.*cos2sigm**2.)   ))
    s = b*A*(sig - delsig)
    a1 = atan2(cos(u2)*sin(Lambda), cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(Lambda) )
    #a2 = atan2(cos(u1)*sin(Lambda), -sin(u1)*cos(u2) + cos(u1)*sin(u2)*cos(Lambda) )    
    #print usq, A, B, delsig, s, a1
 
    if np.isnan(a1):
        a1 = 0.

    #a1 = -a1
    if a1 < 0.:
        a1 = 2.*pi + a1

    return a1


def grc_distance(lat1, lon1, lat2, lon2, result='km'):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees) using the Haversine method.
    Combination of:
    http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
    &
    http://gis.stackexchange.com/questions/29239/calculate-bearing-between-two-decimal-gps-coordinates
    """

    # Convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula for distance:
    # Source: Wiki https://en.wikipedia.org/wiki/Haversine_formula
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    #c =  2.* asin( sqrt( sin((dlat)/2.)**2. + cos(lat1)*cos(lat2)* sin((dlon)/2.)**2. ) )
    r = 6371. # Radius of earth in kilometers. Use 3956 for miles
    
    if dlat == 0. and dlon == 0.:
        return 0.

    # Great circle distance:
    c = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon))

    if result == 'km':
        return c * r
    elif result == 'rad':
        return c


#####################################################################
#                       MAIN PROGRAM                                #
#####################################################################

if __name__ == '__main__':
    
    # ===============================================================
    # 1) DEFINE NETWORK CONFIGURATION
    # ===============================================================
    
    networkpath =       "network/HortonGrid.txt"
    connectionspath =   "network/HortonGrid_Connections.txt"
    
    # Read station and transformer data:
    network = open(networkpath, 'r')
    netdata = network.readlines()
    
    # Define number of nodes:
    nnodes = len(netdata)
    
    # Read data into arrays:
    npf = np.float32
    geolat, geolon = np.zeros(nnodes, dtype=npf), np.zeros(nnodes, dtype=npf)
    country = []
    sitecode, sitename = [], []
    sitenum = np.zeros(nnodes, dtype=np.int32)
    res_earth = np.zeros(nnodes, dtype=npf)
    res_trans = np.zeros(nnodes, dtype=npf)
    # Station-to-index and index-to-station dicts:
    s2i, i2s = {}, {}       
    
    for i in range(nnodes):
        # READ NODE DATA:
        # ---------------
        data = netdata[i].split("\t")
        # Dictionary for station to index:
        s2i[data[2]] = int(data[0]) - 1
        # Dictionary for index to station:
        i2s[int(data[0]) - 1] = data[2]
        # Site number:
        sitenum[i] = int(data[0]) - 1           # -1 to simplify python indices
        # Site readable names:
        sitename.append(data[1])
        # Site code names:
        sitecode.append(data[2])
        # Site country:
        country.append(data[3])
        # Geographic latitude of node:
        geolat[i] = float(data[4])
        # Geographic longitude of node:
        geolon[i] = float(data[5])
        # Earthing resistance of each node:
        res_earth[i] = float(data[6])
        # Transformer resistances (see notes):
        res_trans[i] = float(data[7])
        
    # READ CONNECTION DATA:
    # ---------------------
    connectionsfile = open(connectionspath, 'r')
    conndata = connectionsfile.readlines()

    # Number of connections:
    nconnects = len(conndata)

    nodefrom = np.zeros(nconnects, dtype=np.int32)
    nodeto = np.zeros(nconnects, dtype=np.int32)
    res_line = np.zeros(nconnects, dtype=npf)
    voltage_lines = np.zeros(nconnects, dtype=np.int32)

    for i in range(0, nconnects):
        conns = conndata[i].split("\t")
        # Connection starts at site:
        nodefrom[i] = int(s2i[conns[1]])
        # Connection ends at site:
        nodeto[i] = int(s2i[conns[2]])
        # Line resistance between connecting sites:
        # (Divide by 3 for full transformer representation as paper values are given per phase)
        res_line[i] = float(conns[5])/3.
        # Voltage level of line:
        voltage_lines[i] = float(conns[6])
        
    # Set inf values to very high resistance for computation:
    valInf = 1.e9
    res_trans[res_trans==np.inf] = valInf
    res_earth[res_earth==np.inf] = valInf
    res_line[res_line==np.inf] = valInf
                                         
    # ===============================================================
    # 2) DEFINE LOCATION VARIABLES
    # ===============================================================
    
    # Boundaries of geographic box:
    nbound = 35.0
    ebound = -80.0
    sbound = 32.0
    wbound = -88.0
    
    # Number of cells in box (roughly equidistant):
    x_inc, y_inc = 6./60., 8./60.
    
    # Latitude and longitude range:
    lat = np.arange(sbound, nbound, x_inc)
    lon = np.arange(wbound, ebound, y_inc)
    
    # TODO: Add in method to remove station

    # ===============================================================
    # 3) ESTABLISH MATRICES FOR LP1985 METHOD
    # ===============================================================
    
    # Define resistance matrices:
    resis = np.zeros((nnodes, nnodes), dtype=npf)
    connections = np.zeros((nnodes,nnodes), dtype=npf)
    for i in range(0,nconnects):
        x, y = int(nodefrom[i]), int(nodeto[i])
        if (resis[x,y] > 0. and res_line[i] > 0.):
            # If res already added to this line, add another line in parallel:
            resis[x,y] = 1./(1./resis[x,y] + 1./res_line[i])
        else:
            resis[x,y] = res_line[i]
        if resis[x,y] > 0.:
            connections[x,y] = 1./resis[x,y]
            connections[y,x] = 1./resis[x,y]

    # Calculate matrix of distance between each point:
    dists = np.zeros((nnodes, nnodes), dtype=npf)
    azi = np.zeros((nnodes, nnodes), dtype=npf)
    for i in range(0,nnodes):
        for j in range(0,nnodes):
            lati, loni = geolat[i], geolon[i]
            latj, lonj = geolat[j], geolon[j]
            if resis[i, j] != 0.:
                dists[i, j] = grc_distance(lati, loni, latj, lonj, result='km')
                dists[j, i] = grc_distance(lati, loni, latj, lonj, result='km')
                try:
                    azi[i, j] = grc_azimuth([loni, lati], [lonj, latj])
                    azi[j, i] = grc_azimuth([loni, lati], [lonj, latj])
                except ZeroDivisionError:
                    azi[i, j] = np.nan
                    azi[j, i] = np.nan

    # Create earth impedance matrix:
    # LP1984 eq. (3): **Z**
    # (assuming that nodes are spaced far enough apart to not affect one another)
    earthimp = np.diag(res_earth + res_trans)

    # Calculate network admittance matrix:
    # LP1984 eq. (10): **Y**
    netadmit = -1.*connections + np.diag(sum(connections))

    # Create system matrix (1+YZ), which needs to be inverted:
    # LP1984 eq. (part of 12): **1 + YZ**
    systemmat = np.dot(netadmit, earthimp) + np.identity(nnodes)

    # ===============================================================
    # 4) DEFINE GEOELECTRIC FIELD
    # ===============================================================
    
    testfield = True
    layermodel = None
    surfacemodel = None
           
    # Option 1: Use 1 V/km field:
    N1Vkm = True
    E1Vkm = True
    if testfield:
        en, ee = np.zeros((len(lat),len(lon)), dtype=npf), np.zeros((len(lat),len(lon)), dtype=npf)
        if N1Vkm:
            en.fill(1.)
        if E1Vkm:
            ee.fill(1.)
        en_int = interpolate.interp2d(lon, lat, en)
        ee_int = interpolate.interp2d(lon, lat, ee)
        
    # Option 2: Load in external field file:
    # TODO: Add this in!
            
    # ===============================================================
    # 5) INTEGRATE FIELD ALONG LINES
    # ===============================================================
        
    # Number of steps in path integration (more is slower but more exact):
    steps = 200
    pathlatsteps, pathlonsteps = np.zeros((nconnects,steps), dtype=npf), np.zeros((nconnects,steps), dtype=npf)
    Vn_tot, Ve_tot = np.zeros(nconnects, dtype=npf), np.zeros(nconnects, dtype=npf)
    E_e, E_n = np.zeros((nconnects,steps), dtype=npf), np.zeros((nconnects,steps), dtype=npf)

    for i in range(0,nconnects):
        slat = geolat[nodefrom[i]]
        slon = geolon[nodefrom[i]]
    
        flat = geolat[nodeto[i]]
        flon = geolon[nodeto[i]]

        # Set number of steps in path integration:
        steplat = (flat - slat)/float(steps)
        steplon = (flon - slon)/float(steps)

        if steplat != 0.:
            pathlatsteps[i,:] = np.arange(slat, flat-steplat*0.1, steplat)
        else:
            pathlatsteps[i,:] = slat

        if steplon != 0.:
            pathlonsteps[i,:] = np.arange(slon, flon-steplon*0.1, steplon)
        else:
            pathlonsteps[i,:] = slon

        for p in range(0,len(pathlonsteps[i,:])):
            E_n[i,p] = en_int(pathlonsteps[i,p], pathlatsteps[i,p])
            E_e[i,p] = ee_int(pathlonsteps[i,p], pathlatsteps[i,p])

        # Integrate to get V = int(E*dL), use cylindrical coordinates:
        intline = dists[nodefrom[i], nodeto[i]]
        intazi = azi[nodefrom[i], nodeto[i]]
        for j in range(steps-1):
            # For a North field:
            vnseg = (0.5 * ( E_n[i,j] + E_n[i,j+1] ) * cos(intazi)*(intline/steps) )
            # For an East field:
            veseg = (0.5 * ( E_e[i,j] + E_e[i,j+1] ) * sin(intazi)*(intline/steps) )
            
            Vn_tot[i] += vnseg
            Ve_tot[i] += veseg
            
    # ===============================================================
    # 6) USE V AND SYSTEM MATRIX TO DETERMINE CURRENT (J)
    # ===============================================================

    # Create Nvoltage matrix for northward field:
    Nvoltage = np.zeros((nnodes, nnodes), dtype=npf)
    for l in range(nconnects):
        Nvoltage[nodefrom[l],nodeto[l]] = Vn_tot[l]

    # Create Evoltage matrix for eastward field:
    Evoltage = np.zeros((nnodes, nnodes), dtype=npf)
    for l in range(nconnects):
        Evoltage[nodefrom[l],nodeto[l]] = Ve_tot[l]
    
    # Use Ohm's law to calculate the current along each pathlength:
    # LP1985 eq. (14): J = V/R
    Nlinecurr = np.zeros((nnodes, nnodes), dtype=npf)
    for m in range(nnodes):
        for n in range(nnodes):
            if resis[m,n] > 0.:
                newval = Nvoltage[m,n]/resis[m,n]
                if not np.isnan(newval):
                    Nlinecurr[m,n] = newval
                else:
                    Nlinecurr[m,n] = 0.

    Elinecurr = np.zeros((nnodes, nnodes), dtype=npf)
    for m in range(nnodes):
        for n in range(nnodes):
            if resis[m,n] > 0.:
                newval = Evoltage[m,n]/resis[m,n]
                if not np.isnan(newval):
                    Elinecurr[m,n] = newval
                else:
                    Elinecurr[m,n] = 0.
                    
    # Total line current from both components:
    Tlinecurr = Nlinecurr + Elinecurr

    # Total current at each node due to the E field:
    Nsourcevec = np.sum(Nlinecurr, axis=0) - np.transpose(np.sum(Nlinecurr, axis=1))
    Esourcevec = np.sum(Elinecurr, axis=0) - np.transpose(np.sum(Elinecurr, axis=1))

    # GIC formed at each node due to the E field:
    # LP1985 eq. (part of 12): **(1 + YZ)^(-1)*J**
    netconN, resid, rank, s = np.linalg.lstsq(systemmat, np.transpose(Nsourcevec))
    netconE, resid, rank, s = np.linalg.lstsq(systemmat, np.transpose(Esourcevec))
    netconT = netconN + netconE

    #****************************************************************
    results = [netconN, netconE]        # GIC PER NODE
    #****************************************************************
    
    # Calculate current flowing through transformers (NOT line current):
    Vn_trans = np.multiply(netconN, (res_earth + res_trans))
    Ve_trans = np.multiply(netconE, (res_earth + res_trans))
    Iline_N, Iline_E = np.zeros(nconnects), np.zeros(nconnects)
    for l in range(0,nconnects):
        Iline_N[l] = (Vn_trans[nodefrom[l]]-Vn_trans[nodeto[l]]) / res_line[l]
        Iline_E[l] = (Ve_trans[nodefrom[l]]-Ve_trans[nodeto[l]]) / res_line[l]
        
    #****************************************************************
    transresults = [Iline_N, Iline_E]   # GIC PER TRANSFORMER
    #****************************************************************
        
    # ===============================================================
    # 7) PRINT RESULTS
    # ===============================================================
    
    # Print GIC per node (only grounding currents are non-zero):
    print("")
    print("Network node\tGIC_N    \tGIC_E")
    for i in range(nnodes):
        print("{}  \t{:.2f}    \t{:.2f}".format(sitename[i], results[0][i], results[1][i]))
        
    # For convenience, convert connections to transformer codes:
    transID = {'Sub1_LV-Sub1_E': 'T1',
               'Sub2_LV-Sub2_E': 'T3 /T4',
               'Sub3_HV-Sub3_LV': 'T5 / T15 (series)', 'Sub3_LV-Sub3_E': 'T5 / T15 (common)',
               'Sub4_LV-Sub4_E': 'T2 / T13 (LV), T12 / T14 (common)', 'Sub4_HV-Sub4_LV': 'T12 / T14 (series)',
                    'Sub4_HV-Sub4_E': 'T2 / T13 (HV)',
               'Sub5_LV-Sub5_G': 'T8', 'Sub5_HV-Sub5_G': 'T9',
               'Sub6_HV-Sub6_E': 'T6 / T7',
               'Sub8_HV-Sub8_E': 'T10 / T11'}
        
    # Print GIC per transformer (transformers here being 'lines'):
    print("")
    maxslen = max(len(p) for p in transID.values())
    print("Transformer code\t\t\tGIC_N\tGIC_E")
    for i in range(nconnects):
        codes, codes_r = sitename[nodefrom[i]]+'-'+sitename[nodeto[i]], sitename[nodeto[i]]+'-'+sitename[nodefrom[i]]
        if codes in transID or codes_r in transID:
            try:
                print(("{}".ljust(1+maxslen-len(transID[codes]))+
                       "\t{:.2f}\t{:.2f}").format(transID[codes], Iline_N[i]/3., Iline_E[i]/3.))
            except:
                print(("{}".ljust(1+maxslen-len(transID[codes_r]))+
                       "\t{:.2f}\t{:.2f}").format(transID[codes_r], Iline_N[i]/3., Iline_E[i]/3.))
                       
    #IPython.embed()

