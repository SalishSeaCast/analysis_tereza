#! /usr/bin/env python3

from __future__ import print_function
from numpy import *
from scipy import *
from matplotlib.pyplot import *
#from oct2py import octave

###############################################
# 5.dec.2017 DI - looking at DE for Matt and Paul
# 14.dec.2017 DI - using /Users/dianson/work/data/westcoast/2010_36/dixon/profile_3.py as a template
################################################
# load the data - expects only one file

#infil = loadtxt('../../westcoast/2004_24/2004_24.raw')
infil = loadtxt('CPE.raw')
alefil = loadtxt('/Users/dianson/work/data/westcoast/ale/ale_mother.raw')
#sogfil = loadtxt('../../ale/ale_mother.raw')

###################################
# setting 'vectors' from the data file

crid= infil[:,0]
ln = infil[:,2]
stn = infil[:,3]
P = infil[:,8]
T = infil[:,9]
S = infil[:,10]
ox = infil[:,11]
ox_qf = infil[:,12]
dic = infil[:,13]
alk = infil[:,15]
dic_qf = infil[:,14]
alk_qf = infil[:,16]
no3 = infil[:,17]
no3_qf = infil[:,18]
si = infil[:,19]
si_qf = infil[:, 20]
po4 = infil[:,21]
po4_qf = infil[:, 22]

# alefil - for west coast properties
acrid = alefil[:,0]
astn = alefil[:,3]
aP = alefil[:,8]
aT = alefil[:,9]
aS = alefil[:,10]
aox = alefil[:,11]
aox_qf = alefil[:,12]
adic = alefil[:,13]
aalk = alefil[:,15]
adic_qf = alefil[:,14]
aalk_qf = alefil[:,16]
ano3 = alefil[:,17]
ano3_qf = alefil[:,18]
asi = alefil[:,19]
asi_qf = alefil[:, 20]
apo4 = alefil[:,21]
apo4_qf = alefil[:, 22]

##################################
# create masks for subsetting data, there are a lot here that aren't used in this code
# can you figure out the logic operators?

pres_a750 = (P <= 750.)
ox_qc = ((ox_qf == 2) | (ox_qf == 6)) &  (ox > 0.0) 
dic_qc = ((dic_qf == 2) | (dic_qf == 6)) & (dic > 0.0)
alk_qc = ((alk_qf == 2) | (alk_qf == 6)) & (alk > 0.0)
no3_qc = ((no3_qf == 2) | (no3_qf == 6)) & (no3 > 0.0)
si_qc = ((si_qf == 2) | (si_qf == 6)) & (si> 0.0)
po4_qc = ((po4_qf == 2) | (po4_qf == 6)) & (po4 > 0.0)
# rosette stations for G line 2, 4, 5, 7, 9
# rosette stations for DE line 0, 2, 4, 5, 9
#stn_in = ((stn == 0) | (stn == 2) | (stn == 4) | (stn == 5) | (stn == 9))
stn_in = (crid == 2004.24)
red = ((crid == 2010.36))

mskpo = pres_a750  & ox_qc
mskp = pres_a750
mskpoda = pres_a750  & ox_qc & dic_qc & alk_qc

###ale masks
aox_qc = ((aox_qf == 2) | (aox_qf == 6)) &  (aox > 0.0)
adic_qc = ((adic_qf == 2) | (adic_qf == 6)) & (adic > 0.0)
aalk_qc = ((aalk_qf == 2) | (aalk_qf == 6)) & (aalk > 0.0)
ano3_qc = ((ano3_qf == 2) | (ano3_qf == 6)) & (ano3 > 0.0)
asi_qc = ((asi_qf == 2) | (asi_qf == 6)) & (asi> 0.0)
apo4_qc = ((apo4_qf == 2) | (apo4_qf == 6)) & (apo4 > 0.0)

acrid_out = (acrid == 2007.00) # just to mask out NOAA units that are already umol/kg
acrid_in = ((acrid == 1998.23) | (acrid == 2004.24) | (acrid == 2010.36) | (acrid == 2010.01))

################################
#unit conversions raw data - all to umol/kg

#ox_kg = ox / 0.0223916 / dens * 1000.
#no3_kg = no3 / dens * 1000
#si_kg = si / dens * 1000
#po4_kg = po4 /dens * 1000

#aox_kg = aox / 0.0223916 / adens * 1000.
#ano3_kg = ano3 / adens * 1000
#asi_kg = asi / adens * 1000
#apo4_kg = apo4 /adens * 1000

###################################

###################################
# an example of printing out 
# stations in this example to help set masking - so you know what you have


print('lines', ln);
print('stations', stn);
print('how many have with good quality dic', stn[dic_qc]);


###################################################
# make some property property figs 
# specifically few of the ones I use to QC data)

f = figure(1)
a1 = f.add_subplot(311)
a2 = f.add_subplot(312)
a3 = f.add_subplot(313)

a1.plot(aS[acrid_in & aox_qc], aox[acrid_in & aox_qc], 'kx', markersize=3)
a1.plot(S[stn_in & ox_qc], ox[stn_in & ox_qc], 'bo', markersize=4)
a1.plot(S[ox_qc & red], ox[ox_qc & red], 'ro', markersize=5)
a1.grid()
a2.plot(aS[adic_qc], adic[adic_qc], 'kx', markersize=3)
a2.plot(S[stn_in & dic_qc], dic[stn_in & dic_qc], 'bo', markersize=4)
a2.plot(S[dic_qc & red], dic[dic_qc & red], 'ro', markersize=5)
a2.grid()
a3.plot(aS[aalk_qc], aalk[aalk_qc], 'kx', markersize=3)
a3.plot(S[stn_in & alk_qc], alk[stn_in & alk_qc], 'bo', markersize=4)
a3.plot(S[alk_qc & red], alk[alk_qc & red], 'ro', markersize=5)
a3.grid()

savefig('fig_prop3.pdf')
show()

###############################
# some profiles  
# shows one way to set axis limits and flip the y axis to be surface to deep 

P1 = 0;
P2 = 110;
n1 = 0;
n2 = 30; 
dic1 = 1800;
dic2 = 2300;
ox1 = 2;
ox2 = 8;

f = figure(2)
a4 = f.add_subplot(131)
a5 = f.add_subplot(132)
a6 = f.add_subplot(133)

# first profile
a4.plot(aS, aP, 'kx', markersize=1)
a4.plot(S[stn_in], P[stn_in], 'bo', markersize=3)
a4.plot(S[red], P[red], 'ro', markersize=5)
a4.grid()
a4.set_ylim([P1, P2])
a4.set_ylim(a4.get_ylim()[::-1])
#
# second profile
# looking at no3- profiles, next two lines
#a5.plot(no3[stn_in & no3_qc], P[stn_in & no3_qc], 'bo', markersize=3)
#a5.plot(no3[stn_in & no3_qc & red], P[stn_in & no3_qc & red], 'ro', markersize=5)
#a5.set_xlim([n1, n2])
#
# looking at dic, next two lines
a5.plot(adic[adic_qc], aP[adic_qc], 'kx', markersize=1)
a5.plot(dic[stn_in & dic_qc], P[stn_in & dic_qc], 'bo', markersize=3)
a5.plot(dic[dic_qc & red], P[dic_qc & red], 'ro', markersize=5)
a5.grid()
a5.set_xlim([dic1, dic2])
# the next two lines end of second profile common to both no3 and dic choice
a5.set_ylim([P1, P2])
a5.set_ylim(a5.get_ylim()[::-1])
#
# third profile
a6.plot(aox[aox_qc & acrid_in], aP[aox_qc & acrid_in], 'kx', markersize=1)
a6.plot(ox[stn_in & ox_qc], P[stn_in & ox_qc], 'bo', markersize=3)
a6.plot(ox[ox_qc & red], P[ox_qc & red], 'ro', markersize=5)
a6.grid()
a6.set_xlim([ox1, ox2])
a6.set_ylim([P1, P2])
a6.set_ylim(a6.get_ylim()[::-1])

savefig('fig_prof3.pdf')
show()
