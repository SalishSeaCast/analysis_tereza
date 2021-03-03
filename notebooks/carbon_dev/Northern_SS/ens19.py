start ='2019-01-01'
end ='2019-12-31'
yrstr = '19'

#######
import numpy as np
import glob
import warnings
import pickle
import arrow
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import netCDF4 as nc
import gsw

#import LambertConformalTicks as lct

ns_x = 160; ns_y = 680; ns_col = 'olive'
ns2_x = 150; ns2_y = 640; ns2_col = 'yellowgreen'
ns3_x = 155; ns3_y = 710; ns3_col = 'palegoldenrod'
ns4_x = 145; ns4_y = 740; ns4_col = 'lawngreen'


############
def OmA_2D(grid,carp):
    tsal = grid['vosaline'][0,0,:,:]
    ttemp = grid['votemper'][0,0,:,:]
    tdic = carp['dissolved_inorganic_carbon'][0,0,:,:]
    tta = carp['total_alkalinity'][0,0,:,:]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(898,398)
    OmAr = OmegaA.reshape(898,398)
    OmCr = OmegaC.reshape(898,398)
    pco2r = pco2.reshape(898,398)
    
    return pHr, OmAr, OmCr, pco2r

####

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array1 = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array1.append(r)
dayslen = len(arrow_array1)

pco2_2019 = np.zeros([dayslen,4])
nit_10m_2019 = np.zeros([dayslen,4])
diat_10m_2019 = np.zeros([dayslen,4])
flag_10m_2019 = np.zeros([dayslen,4])
microzoo_10m_2019 = np.zeros([dayslen,4])

for i in range(0,dayslen):

    tdate = arrow_array1[i][0]
    ymd = tdate.format('YYYYMMDD')
    DD = tdate.format('DD')
    mon = tdate.format('MMM').lower()
    yd = f'{DD}{mon}{yrstr}'
    
    if i%10 == 0:
        print(ymd)
    tstr = glob.glob(f'/results2/SalishSea/nowcast-green.201905/{yd}/SalishSea_1d_*{ymd}*carp_T.nc')
    tnc = tstr[0]; tn_carp = nc.Dataset(tnc)
    
    tstr = glob.glob(f'/results2/SalishSea/nowcast-green.201905/{yd}/SalishSea_1d_*{ymd}*grid_T.nc')
    tnc = tstr[0]; tn_grid = nc.Dataset(tnc)
    pHr, OmAr, OmCr, pco2r = OmA_2D(tn_grid,tn_carp)
    #print(np.max(pco2r))
    pco2r[pco2r>1e10] = np.nan
    tstr = glob.glob(f'/results2/SalishSea/nowcast-green.201905/{yd}/SalishSea_1d_*{ymd}*ptrc_T.nc')
    tnc = tstr[0]; tn_ptrc = nc.Dataset(tnc)
    pco2_2019[i,0] = np.nanmean(pco2r[ns_y-10:ns_y+10,ns_x-10:ns_x+10])
    pco2_2019[i,1] = np.nanmean(pco2r[ns2_y-10:ns2_y+10,ns2_x-10:ns2_x+10])
    pco2_2019[i,2] = np.nanmean(pco2r[ns3_y-10:ns3_y+10,ns3_x-10:ns3_x+10])
    pco2_2019[i,3] = np.nanmean(pco2r[ns4_y-10:ns4_y+10,ns4_x-10:ns4_x+10])
    
    nit = tn_ptrc['nitrate'][:]; nit[nit == 0] = np.nan
    diat = tn_ptrc['diatoms'][:]; diat[diat == 0] = np.nan
    flag = tn_ptrc['flagellates'][:]; flag[flag == 0] = np.nan
    microzoo = tn_ptrc['microzooplankton'][:]; microzoo[microzoo == 0] = np.nan
    
    nit_10m_2019[i,0] = np.nanmean(nit[0,0:10,ns_y-10:ns_y+10,ns_x-10:ns_x+10])
    nit_10m_2019[i,1] = np.nanmean(nit[0,0:10,ns2_y-10:ns2_y+10,ns2_x-10:ns2_x+10])
    nit_10m_2019[i,2] = np.nanmean(nit[0,0:10,ns3_y-10:ns3_y+10,ns3_x-10:ns3_x+10])
    nit_10m_2019[i,3] = np.nanmean(nit[0,0:10,ns4_y-10:ns4_y+10,ns4_x-10:ns4_x+10])

    diat_10m_2019[i,1] = np.nanmean(diat[0,0:10,ns2_y-10:ns2_y+10,ns2_x-10:ns2_x+10])
    diat_10m_2019[i,0] = np.nanmean(diat[0,0:10,ns_y-10:ns_y+10,ns_x-10:ns_x+10])
    diat_10m_2019[i,2] = np.nanmean(diat[0,0:10,ns3_y-10:ns3_y+10,ns3_x-10:ns3_x+10])
    diat_10m_2019[i,3] = np.nanmean(diat[0,0:10,ns4_y-10:ns4_y+10,ns4_x-10:ns4_x+10])

    flag_10m_2019[i,1] = np.nanmean(flag[0,0:10,ns2_y-10:ns2_y+10,ns2_x-10:ns2_x+10])
    flag_10m_2019[i,0] = np.nanmean(flag[0,0:10,ns_y-10:ns_y+10,ns_x-10:ns_x+10])
    flag_10m_2019[i,2] = np.nanmean(flag[0,0:10,ns3_y-10:ns3_y+10,ns3_x-10:ns3_x+10])
    flag_10m_2019[i,3] = np.nanmean(flag[0,0:10,ns4_y-10:ns4_y+10,ns4_x-10:ns4_x+10])

    microzoo_10m_2019[i,1] = np.nanmean(microzoo[0,0:10,ns2_y-10:ns2_y+10,ns2_x-10:ns2_x+10])
    microzoo_10m_2019[i,0] = np.nanmean(microzoo[0,0:10,ns_y-10:ns_y+10,ns_x-10:ns_x+10])
    microzoo_10m_2019[i,2] = np.nanmean(microzoo[0,0:10,ns3_y-10:ns3_y+10,ns3_x-10:ns3_x+10])
    microzoo_10m_2019[i,3] = np.nanmean(microzoo[0,0:10,ns4_y-10:ns4_y+10,ns4_x-10:ns4_x+10])
            
pickle.dump(nit_10m_2019, open("./pkls/nit_10m_2019.pkl", 'wb'))
pickle.dump(diat_10m_2019, open("./pkls/diat_10m_2019.pkl", 'wb'))
pickle.dump(pco2_2019, open("./pkls/pco2_2019.pkl", 'wb'))
pickle.dump(microzoo_10m_2019, open("./pkls/microzoo_10m_2019.pkl", 'wb'))
pickle.dump(flag_10m_2019, open("./pkls/flag_10m_2019.pkl", 'wb'))



