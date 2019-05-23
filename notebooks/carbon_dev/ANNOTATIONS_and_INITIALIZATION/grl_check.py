#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from oct2py import octave
import seawater as sw


#============================================================================#
#                              LOAD DATA                                     #
#============================================================================#
infil = np.loadtxt('grl2016.raw')

crid = infil[:,0]
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
si_qf = infil[:,20]
po4 = infil[:,21]
po4_qf = infil[:,22]
pH_obs = infil[:,23]
pH_obs_qf = infil[:,24]


#============================================================================#
#                         CALCULATE DERIVED PROPERTIES                       #
#============================================================================#

# density vector and sigma theta
dens = sw.dens(S, T, P)
sigmaT = 1000. - dens  
ptmp = sw.ptmp(S, T, P, 0.)
dens_ze = sw.dens(S, ptmp, 0.)
sigma_theta = dens_ze - 1000.

# unit conversions raw data - all to umol/kg
ox_kg = ox / 0.0223916 / dens * 1000.
no3_kg = no3 / dens * 1000.
si_kg = si / dens * 1000.
po4_kg = po4 /dens * 1000.


#============================================================================#
#                         CREATE DATA MASKS                                  #
#============================================================================#

ox_qc = ((ox_qf == 2) | (ox_qf == 6)) &  (ox > 0.0)
dic_qc = ((dic_qf == 2) | (dic_qf == 6) | (dic_qf == 10)) & (dic > 0.0)
alk_qc = ((alk_qf == 2) | (alk_qf == 6) | (alk_qf == 10)) & (alk > 0.0)
no3_qc = ((no3_qf == 2) | (no3_qf == 6)) & (no3 > 0.0)
si_qc = ((si_qf == 2) | (si_qf == 6)) & (si> 0.0)
po4_qc = ((po4_qf == 2) | (po4_qf == 6)) & (po4 > 0.0)

carb_pair_qc = ( ((dic_qf == 2) & (alk_qf == 2)) | ((dic_qf == 6) & (alk_qf == 6)) | 
                  ((dic_qf == 12) & (alk_qf == 10)) | ((dic_qf == 10) & (alk_qf == 12)) |
                  ((dic_qf == 10) & (alk_qf == 10)) )

#

# mask for CO2SYS
#msk_C = dic_qc & alk_qc & po4_qc & si_qc
msk_C = carb_pair_qc & po4_qc & si_qc
#msk_C = po4_qc & si_qc 

#msk for plotting just alk and just dic
# all data
msk_dic = dic_qc 
msk_alk = alk_qc 
msk_ox = ox_qc 


########################################
#create the CO2SYS vectors - the max size for the omega and pH vectors

S_C   = S[msk_C]
T_C   = T[msk_C]
P_C   = P[msk_C]
ox_C  = ox[msk_C]
ox_qf_C = ox_qf[msk_C]
dic_C = dic[msk_C] 
alk_C = alk[msk_C]
#dic and alk qc will be all true by construction, same for po4 and si
#dic_qf_C = dic_qf[msk_C]
#alk_qf_C = alk_qf[msk_C]
no3_C =  no3[msk_C]
no3_qf_C = no3_qc[msk_C]
si_C = si[msk_C]
#si_qf_C = si_qf[msk_C]
po4_C = po4[msk_C]
#po4_qf_C = po4_qc[msk_C]
crid_C = crid[msk_C]
stn_C = stn[msk_C]
ln_C = ln[msk_C]

dens_C = dens[msk_C]
sigma_theta_C = sigma_theta[msk_C]

ln_in_CPE_C = (ln_C == -1)
ln_in_north_C = (ln_C == 1)
ln_in_south_C = (ln_C == 2)
ln_in_fraser_C = (ln_C == 3)
ln_in_haro_C = (ln_C == 4)
ln_in_jdf_C = (ln_C == 5)

###### and also the unit converted vectors

ox_C_kg = ox_kg[msk_C]
no3_C_kg = no3_kg[msk_C]
si_C_kg = si_kg[msk_C]
po4_C_kg = po4_kg[msk_C] 

#################### for deep water intrusion point

mskdwi = ((crid_C == 2012.57) & (stn_C == 39.5) & (dic_C == 2040.64) )

################################
# calculate carbonate parameters using co2sys

# with '4' constants which are Roberta's choice 
#Result, Headers = octave.CO2SYS(dic, alk, 2., 1., S, T, T, P, P, si_kg, po4_kg, 1., 4., 1.)

# with Leuker et al. 2010 constants ('10')
# with Millero et al. 2006 constants ('13')
# with Millero et al. 2010 constants ('14')
Result = octave.CO2SYS(dic_C, alk_C, 2., 1., S_C, T_C, T_C, P_C, P_C, si_C_kg, po4_C_kg, 1., 14., 1.)

#omega_arag = Result[:,30]  #checked against Roberta

#print(omega_arag)

# total scale
pH = Result[:,36]  #checked against Roberta

#print(pH)

# some more while I am at it
#co3 = Result[:,21]
#omega_calcite = Result[:,29]
pco2 = Result[:,18]
#print(min(pco2))

##################################
# calculating AOU 

#ptmp = octave.sw_ptmp(S, T, P, 0.)[0]
#surface_dens = octave.sw_dens(S, ptmp, 0.)[0] / 1000. # kg / L

#O2 = ox / 0.0223916 / surface_dens
#AOU = octave.O2sol(S, ptmp)[0] - O2

#AOU_C = AOU[msk_C]

# print(AOU) # - checked - matches Roberta

#============================================================================#
#                              CREATE PLOTS                                  #
#============================================================================#

# define some colours
harblue=(0.2, 0., 0.7)
purple=(0.3, 0., 0.5)
dpurple=(0.2, 0., 0.3)
lpurple=(0.9, 0., 0.9)
orange=(0.9, 0.6, 0.)

# set plot output parameters
params = {\
#  'figure.figsize'        : (7.48, 9.0551), #fullpg AGU 190mmX230mm
#  'figure.figsize'        : (3.74, 9.0551), #exactly halfpage AGU
#skinny is 95mm wide
  'figure.figsize'        : (3.74, 9.0551),
    'figure.subplot.left'   : 0.16,
# was 0.12 but need more space for font
    'figure.subplot.right'  : 0.92,
# was 0.88 but wanted plots a little less square
#    'font.family'           : 'serif',
#    'font.serif'            : 'Times',
#    'font.sans-serif'       : 'Helvetica',
#    'font.cursive'          : 'Zapf Chancery',
#    'font.monospace'        : 'Courier',
#     'text.usetex'           : True,
    'axes.labelsize'        : 8,
    'font.size'             : 8,
    'legend.fontsize'       : 8,
    'legend.numpoints'      : 1,
    'xtick.labelsize'       : 8,
    'ytick.labelsize'       : 8,
    'savefig.dpi'           : 600,
    }

plt.rcParams.update(params)

s1 = 26  # was 27 #was 24
s2 = 34.1 # was 34
dic1 = 1700 #1600
dic2 = 2300 #2400
alk1 = 1700
alk2 = 2400
ox1 = 50
ox2 = 500 # was 400
ph1 = 7.4
ph2 = 8.2
sigma_theta1 = 20
sigma_theta2 = 26.8 # was 26

f = plt.figure(1)

a = f.add_subplot(311)
a2 = f.add_subplot(312)
a3 = f.add_subplot(313)

a.plot([30.3, 30.3],[dic1, dic2], linewidth = 2, color = '0.75')
p1, = a.plot(S_C[ln_in_jdf_C], dic_C[ln_in_jdf_C], 'ko', markersize=4)
p2, = a.plot(S_C[ln_in_haro_C], dic_C[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
p4, = a.plot(S_C[ln_in_north_C], dic_C[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a.plot(S_C[mskdwi], dic_C[mskdwi], 'gx', markersize=8)
a.plot(S_C[mskdwi], dic_C[mskdwi], 'g+', markersize=8)
p3, = a.plot(S_C[ln_in_south_C], dic_C[ln_in_south_C], 'go', markersize=4)
#p2, = a.plot(S_C[ln_in_haro_C], dic_C[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
#a.plot(S_C[ln_in_jdf_C], dic_C[ln_in_jdf_C], 'ko', markersize=6)
a.set_xlim([s1, s2])
a.set_ylim([dic1, dic2])
#a.plot([30.3, 30.3],[dic1, dic2], linewidth = 2, color = '0.75')
#a.set_xlabel('S')
a.set_ylabel('DIC ($\mu$mol kg$^{-1}$)')
a.legend([p1, p2, p3, p4],["Juan de Fuca", "Haro", "South-SoG", "North-SoG"], loc = 4)
a.text(24.5, 2240, 'a)')
a.grid()

a2.plot([30.3, 30.3],[ph1, ph2], linewidth = 2, color = '0.75')
a2.plot(S_C[ln_in_jdf_C], pH[ln_in_jdf_C], 'ko', markersize=4)
a2.plot(S_C[ln_in_haro_C], pH[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a2.plot(S_C[ln_in_north_C], pH[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a2.plot(S_C[mskdwi], pH[mskdwi], 'gx', markersize=8)
a2.plot(S_C[mskdwi], pH[mskdwi], 'g+', markersize=8)
a2.plot(S_C[ln_in_south_C], pH[ln_in_south_C], 'go', markersize=4)
#a2.plot(S_C[ln_in_jdf_C], pH[ln_in_jdf_C], 'ko', markersize=4)
a2.set_xlim([s1, s2])
a2.set_ylim([ph1, ph2])
#a2.plot([30.3, 30.3],[ph1, ph2], linewidth = 2, color = '0.75')
#a2.set_xlabel('S')
a2.set_ylabel('pH')
a2.text(24.5, 8.11, 'b)')
a2.grid()

a3.plot([30.3, 30.3],[ox1, ox2], linewidth = 2, color = '0.75')
#a3.plot(S_C[ln_in_south_C], ox_C_kg[ln_in_south_C], 'go', markersize=4)
a3.plot(S_C[ln_in_jdf_C], ox_C_kg[ln_in_jdf_C], 'ko', markersize=4)
a3.plot(S_C[ln_in_haro_C], ox_C_kg[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a3.plot(S_C[ln_in_north_C], ox_C_kg[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a3.plot(S_C[mskdwi], ox_C_kg[mskdwi], 'gx', markersize=8)
a3.plot(S_C[mskdwi], ox_C_kg[mskdwi], 'g+', markersize=8)
a3.plot(S_C[ln_in_south_C], ox_C_kg[ln_in_south_C], 'go', markersize=4)
#a3.plot(S_C[ln_in_haro_C], ox_C_kg[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a3.set_xlim([s1, s2])
#a3.plot([30.3, 30.3],[ox1, ox2], linewidth = 2, color = '0.75')
a3.set_xlabel('Salinity')
a3.set_ylabel('O$_2$ ($\mu$mol kg$^{-1}$)')
a3.text(24.5, 455, 'c)')
a3.grid()

#savefig('fig.pdf')
#show()

f = plt.figure(13)

a4 = f.add_subplot(311)
a5 = f.add_subplot(312)
a6 = f.add_subplot(313)

a4.plot([23.5, 23.5],[dic1, dic2], linewidth = 2, color = '0.75')
p5, = a4.plot(sigma_theta_C[ln_in_jdf_C], dic_C[ln_in_jdf_C], 'ko', markersize=4)
p6, = a4.plot(sigma_theta_C[ln_in_haro_C], dic_C[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
p8, = a4.plot(sigma_theta_C[ln_in_north_C], dic_C[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a4.plot(sigma_theta_C[mskdwi], dic_C[mskdwi], 'gx', markersize=8)
a4.plot(sigma_theta_C[mskdwi], dic_C[mskdwi], 'g+', markersize=8)
p7, = a4.plot(sigma_theta_C[ln_in_south_C], dic_C[ln_in_south_C], 'go', markersize=4)
#p6, = a4.plot(sigma_theta_C[ln_in_haro_C], dic_C[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a4.set_xlim([sigma_theta1, sigma_theta2])
a4.set_ylim([dic1, dic2])
#a4.plot([23.5, 23.5],[dic1, dic2], linewidth = 2, color = '0.75')
#a4.set_xlabel('Sigma_${theta}$')
a4.set_ylabel('DIC ($\mu$mol kg$^{-1}$)')
a4.legend([p5, p6, p7, p8],["Juan de Fuca", "Haro", "South-SoG", "North-SoG"], loc = 4)
a4.text(18.688, 2240, 'a)')
a4.grid()

a5.plot([23.5, 23.5],[ph1, ph2], linewidth = 2, color = '0.75')
a5.plot(sigma_theta_C[ln_in_jdf_C], pH[ln_in_jdf_C], 'ko', markersize=4)
a5.plot(sigma_theta_C[ln_in_haro_C], pH[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a5.plot(sigma_theta_C[ln_in_north_C], pH[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a5.plot(sigma_theta_C[mskdwi], pH[mskdwi], 'gx', markersize=8)
a5.plot(sigma_theta_C[mskdwi], pH[mskdwi], 'g+', markersize=8)
a5.plot(sigma_theta_C[ln_in_south_C], pH[ln_in_south_C], 'go', markersize=4)
#a5.plot(sigma_theta_C[ln_in_haro_C], pH[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a5.set_xlim([sigma_theta1, sigma_theta2])
a5.set_ylim([ph1, ph2])
#a5.plot([23.5, 23.5],[ph1, ph2], linewidth = 2, color = '0.75')
#a5.set_xlabel('Sigma$_{theta}$')
a5.set_ylabel('pH')
a5.text(18.688, 8.11, 'b)')
a5.grid()

a6.plot([23.5, 23.5],[ox1, ox2], linewidth = 2, color = '0.75')
#a6.plot(sigma_theta_C[ln_in_south_C], ox_C_kg[ln_in_south_C], 'go', markersize=4)
a6.plot(sigma_theta_C[ln_in_jdf_C], ox_C_kg[ln_in_jdf_C], 'ko', markersize=4)
a6.plot(sigma_theta_C[ln_in_haro_C], ox_C_kg[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a6.plot(sigma_theta_C[ln_in_north_C], ox_C_kg[ln_in_north_C], c=orange, marker='o', markersize=4, linestyle='None')
a6.plot(sigma_theta_C[mskdwi], ox_C_kg[mskdwi], 'gx', markersize=8)
a6.plot(sigma_theta_C[mskdwi], ox_C_kg[mskdwi], 'g+', markersize=8)
a6.plot(sigma_theta_C[ln_in_south_C], ox_C_kg[ln_in_south_C], 'go', markersize=4)
#a6.plot(sigma_theta_C[ln_in_haro_C], ox_C_kg[ln_in_haro_C], c=lpurple, marker='o', markersize=4, linestyle='None')
a6.set_xlim([sigma_theta1, sigma_theta2])
#a6.plot([23.5, 23.5],[ox1, ox2], linewidth = 2, color = '0.75')
#a6.set_xlabel('$\sigma \Theta}$') #works
a6.set_xlabel('$\sigma_\Theta}$')
a6.set_ylabel('O$_2$ ($\mu$mol kg$^{-1}$)')
a6.text(18.688, 455, 'c)')
a6.grid()

#plt.savefig('figtheta.pdf')
plt.show()

print('how many northern data?', len(S_C[ln_in_north_C]))
