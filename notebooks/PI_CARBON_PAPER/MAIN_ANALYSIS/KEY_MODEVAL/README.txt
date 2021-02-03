TJSJ-Jan 24, 2021
Description of folder KEY_MODEVAL contents

-folder Datasets:
	>> either raw data (see readme there) or csvs made from pandas dataframes in KEY_MODEVAL/Dataset_dfmaker.ipynb

-folder extracted_calculated
	>> salspace depthspace model output corresponding to datasets D14, D15, COMB, or GRL (YR...etc) > made in Calculate_Omega_for_obsdata_and_corresponding_model_points.ipynb
	>> also MASKED datasets > made in Calculate_Omega_for_obsdata_and_corresponding_model_points.ipynb
	>>also stats in depth and sal space > made in calculate_statistics	

-folder ncs:
   for the following observational datasets (COMB, GRL, 2014, 2015), ncs either extracting corresponding model point in depthspace (day, month, depth) or salspace (day, month, wherever model salinity corresponds to obs salinity (interpolated values of T,DIC,TA) - note, made with only data corresponding to the following mask: 
        mask = (((dic_qf==2) | (dic_qf==6)) & \
        ((alk_qf==2) | (alk_qf==6)) & \
        (dic > 0) & (alk >0))

   made using either Modeval_DICTA-expandedanalysis_by_depth_EXTRACTER.ipynb
   Modeval_DICTA-expandedanalysis_by_sal_EXTRACTER.ipynb

-individual ipynbs:
	>> Calculate_Omega_for_obsdata_and_corresponding_model_points.ipynb
		- observational omega, model omega (in depth and salinity space), also add locational tags > goes in extracted_calculated
		- also convert DIC and TA to model units and  
	>> calculate_statistics.ipynb
		- calculate stats bias, RMSE, WSS, stdrat_MtoO, also bias with depth
	>> Dataset_dfmaker 
		- take raw data in .txt files and make dataframes, save as .csv (still in Datasets)
	>> Modeval_DICTA-expandedanalysis_by_depth_EXTRACTER.ipynb
		- find data corresponding to obs by depth and save in .nc
	>> Plot_statistics_DIC (or OmegaA, or TA) 
		- plot s 
