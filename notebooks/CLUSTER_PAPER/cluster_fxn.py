def cluster_patterns(orig_data,cluster_des,cluster_no,noday):
    import numpy as np
    
    '''For a given data matrix (S stations X 365 daily signals), list of clusters designation (S stations), 
    and specific cluster, return an C_S x 365 matrix of the C_S annual signals in that cluster
    for the C_S stations in the cluster'''
    #which stations are in the cluster we are looking for?
    where_cluster = np.where(cluster_des == cluster_no)
    where_cluster = np.squeeze(where_cluster)
    #print(where_cluster.shape)
    #print(where_cluster.size)

    
    no_stns_in_cluster = where_cluster.size
    if no_stns_in_cluster == 1:
        this_stn = np.squeeze(where_cluster)
        where_cluster = this_stn
        signalmat = orig_data[this_stn,:]
       
    else:    
    
        signalmat = np.zeros([no_stns_in_cluster,noday])

        for stn in range(0,no_stns_in_cluster):

            this_stn = where_cluster[stn]
            signalmat[stn,:] = orig_data[this_stn,:]
    return where_cluster, signalmat

def find_mean_patterns(signal_mat,no_clusters,cluster_list,noday):
    
    '''Extract mean patterns for each of C clusters in the signal matrix of size S (total no. stns) x 365 (days)
        given a list of cluster numbers
        developed for the salinity matrix
        also compute standard deviation of cluster
        Return 2 matrices of size C x 365 - mean patterns, std_dev patterns
        sample use case: mean_patterns, std_dev_patterns = find_mean_patterns(signal_mat,no_clusters,cluster_list)'''
    import numpy as np
    mean_patterns = np.zeros([no_clusters,noday])
    std_dev_patterns = np.zeros([no_clusters,noday])
    for c in range(0,no_clusters):
        cl = c+1
        
        where_cluster, signalmat = cluster_patterns(signal_mat,cluster_list,cl,noday)
#         print(signalmat)
#         print(signalmat.shape)
        if signalmat.size == noday:
            mean_signal = signalmat
            std_dev_signal = np.zeros_like(signalmat)
        else:
            mean_signal = np.nanmean(signalmat, axis = 0)
            std_dev_signal = np.nanstd(signalmat, axis = 0)
        mean_patterns[c,:] = mean_signal
        std_dev_patterns[c,:] = std_dev_signal
    return mean_patterns, std_dev_patterns


def find_mean_patterns_mv(signal_mat,no_clusters,cluster_list,var_to_include, noday):
    
    '''Extract mean patterns for each of C clusters in the signal matrix of size S (total no. stns) x 365 (days) X V variables
        given a list of cluster numbers, and a list of indices in the signal mat of variables to include
        mv means multivariate
        developed for the biology matrix
        also compute standard deviation of cluster
        Return 2 matrices of size C x 365 - mean patterns, std_dev patterns
        sample use case: mean_patterns, std_dev_patterns = find_mean_patterns(signal_mat,no_clusters,cluster_list)'''
    import cluster_fxn as cf
    import numpy as np
    no_var = len(var_to_include)
    
    mean_patterns = np.zeros([no_clusters,noday,no_var])
    #print(mean_patterns.shape)
    std_dev_patterns = np.zeros([no_clusters,noday,no_var])
    for i in range(0,no_var):
    
        this_var = var_to_include[i]
        #print(this_var)
        signal_mat_thisvar = signal_mat[:,:,this_var]
        #print(signal_mat_thisvar.shape)
        for c in range(0,no_clusters):
            cl = c+1
            where_cluster, signalmat = cf.cluster_patterns(signal_mat_thisvar,cluster_list,cl,noday)
            mean_signal = np.nanmean(signalmat, axis = 0)
            std_dev_signal = np.nanstd(signalmat, axis = 0)
            mean_patterns[c,:,i] = mean_signal
            std_dev_patterns[c,:,i] = std_dev_signal
    return mean_patterns, std_dev_patterns

def map_clusters(tit,no_clusters,cl,fsx,fsy,markersize,titfontsize,legfontsize,fname,colors):
    """argmuents: map_clusters(tit,no_clusters,cl,fsx,fsy,markersize,titfontsize,legfontsize)
    cl is list of clusters"""
    import map_fxn as mf
    import matplotlib.pyplot as plt
    import numpy as np
    import cmocean
    from salishsea_tools import (
        nc_tools,
        viz_tools,
        geo_tools,
        tidetools,
        visualisations,
    )
    bath = '/data/tjarniko/MEOPAR/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    spacing = 10
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    no_stns = len(d_stn_x)
    
    fig = plt.figure(figsize=(fsx,fsy))

    plt.rcParams['image.cmap'] = 'YlGnBu'


    
    for i in range(1,2):
        ax = fig.add_subplot(1,1,i)
        viz_tools.set_aspect(ax)   
        fmask = (grid.fmask[0,0,:,:])    
        mesh = ax.pcolormesh(fmask, vmin=0, vmax=1)
        out=ax.set(title='Domain Extent')

        ax.set_ylim([0,898])
        ax.set_xlim([0,398])
        if i ==1:
            #colors = plt.cm.hsv(np.linspace(0, 1, no_clusters))
            for j in range(1,no_clusters+1):
                cluster = np.where(cl == j)
                cluster = np.squeeze(cluster)
                #find the xs and ys of the stations in a given cluster
                c1_x = np.take(d_stn_x,cluster)
                c1_y = np.take(d_stn_y,cluster)
                print(j)
                print(colors[j-1])
                pts = ax.scatter(c1_x,c1_y,s=markersize,c=colors[j-1], label=j ,marker='o')
                #pts = ax.scatter(c1_x,c1_y,s=markersize,c='xkcd:dust', label=j ,marker='o')
            ax.set_title(tit,fontsize= titfontsize)
            ax.set_xticklabels( () ) 
            ax.set_yticklabels( () ) 
            plt.legend(loc=1, fontsize = legfontsize)

    plt.savefig(fname,  transparent=True)
    plt.show()

    
def plot_signals(no_clusters,mean_patterns,std_dev_patterns,ymin,ymax,ylab,fsx,fsy,noday,fname):
    from matplotlib import pyplot as plt
    import numpy as np
    start = 1
    stop = noday+1
    step = 1
    days = np.arange(start,stop,step)

    fig = plt.figure(figsize=(fsx,fsy))

    for i in range(1,no_clusters+1):
        ax = fig.add_subplot(no_clusters,1,i)

        pattern = mean_patterns[i-1,:]
        st_dev = std_dev_patterns[i-1,:]
        ax.fill_between(days,pattern+st_dev,pattern-st_dev, color = 'k', alpha = 0.2)

        ax.plot(days,pattern)

        ax.set_ylim([ymin,ymax])
        ax.set_xlim([1,noday])
        plt.xlabel('days', fontsize=16)
        plt.ylabel(ylab, fontsize=16)
        titstring = ' cluster ' + str(i)
        ax.set_title(titstring, fontsize = 20)
        plt.legend()


    plt.tight_layout()
    plt.savefig(fname)
    plt.show()