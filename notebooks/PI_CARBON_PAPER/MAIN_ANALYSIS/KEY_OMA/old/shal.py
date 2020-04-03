def find_depth_shallowalg(dp,prof,water_depth):
    #finds saturation horizon given a profile and corresponding depths
    first_proper_undersat = np.nan
    depth_undersat = np.nan    
    dummy_var = 0
    #tot masked
    if prof.mask.all():
        dummy_var = 0
        depth_undersat = np.nan
        #print('masks all around!')
    elif np.ma.min(prof) >=1:
        depth_undersat = water_depth
        #print('saturated column')
    elif np.ma.max(prof) <1:
        depth_undersat = 0
        #print('undersat to surface')        
    else:
        t_ind = np.where(prof<1)
        t_indar = t_ind[0][0]
        t_indss = np.where(prof>=1)
        t_indsssar = t_indss[0][0]
        if t_indar.size == 0:
            dummy_var = 0
        else:
            if (t_indar.size != 0) & (t_indsssar.size == 0):
                depth_undersat = 0
                first_proper_undersat = 0
                dummy_var = 0
                #print('undersat to surface!')
                max_supsat = np.nan
            else:    
                max_supsat = np.max(t_indsssar)    
                try:
                    first_proper_undersat = np.min(t_indar)
                except:
                    dummy_var = 0
                    #print("An exception occurred")
                if first_proper_undersat == 0:
                    depth_undersat = dp[0]
                if np.isnan(first_proper_undersat):
                    dummy_var = 0
                    #print('saturated watercolumn!')
                else:
                    depth_undersat = (dp[first_proper_undersat]+dp[first_proper_undersat-1])/2
    return depth_undersat
