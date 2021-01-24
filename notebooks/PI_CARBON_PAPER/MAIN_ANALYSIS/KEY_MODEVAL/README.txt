TJSJ-Jan 24, 2021
Description of folder KEY_MODEVAL contents

-folder Datasets:
   either raw data (see readme there) or csvs made from pandas dataframes in KEY_MODEVAL/Dataset_dfmaker.ipynb

-folder ncs:
   for the following observational datasets (COMB, GRL, 2014, 2015), ncs either extracting corresponding model point in depthspace (day, month, depth) or salspace (day, month, wherever model salinity corresponds to obs salinity (interpolated values of T,DIC,TA) - note, made with only data corresponding to the following mask: 
        mask = (((dic_qf==2) | (dic_qf==6)) & \
        ((alk_qf==2) | (alk_qf==6)) & \
        (dic > 0) & (alk >0))

   made using Modeval_DICTA-expandedanalysis_by_depth_EXTRACTER.ipynb
   Modeval_DICTA-expandedanalysis_by_sal_EXTRACTER.ipynb

