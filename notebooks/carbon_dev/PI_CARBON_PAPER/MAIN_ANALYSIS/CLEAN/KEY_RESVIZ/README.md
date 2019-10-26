The Jupyter Notebooks in this directory are for sharing of Python code
techniques and notes about model results analysis code.
They were developed by Tereza Jarnikova.

The links below are to static renderings of the notebooks via
[nbviewer.jupyter.org](http://nbviewer.jupyter.org/).
Descriptions under the links below are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[OmA_horizon_and_shoaling_maps_and_regional_avgs_using_2algorithms_shallowdeep.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//OmA_horizon_and_shoaling_maps_and_regional_avgs_using_2algorithms_shallowdeep.ipynb)  
    
    **NOTEBOOK CONTENTS**  
      
        -plots of aragonite saturation horizon for summer and winter seasons for BR and PI runs calculated using 2 methods:  
            shallow: arag. hor. is shallowest undersaturat. water  
            deep: arag. hor. is shallowest undersaturat. water beneath all supersaturated water  
        -regional averages as bar plot  
      
    **Leading questions:**  
      
    **What *are* the saturation horizons? **  
      
        Present day values:  
        -Jdf: ~45m summer / 90m winter  
        -Nsog ~30m s / ~20m w  
        -Csog ~40m s / ~25m w  
        -Haro 60m s / 65m w  
          
    **Where do they shoal the most?**  
      
        The inner strait (~10m summer, 20-30 m winter)  
        vs Haro/Jdf (<5 m both summer & winter)  
      
    **Does it matter which algorithm you use?**  
      
        Once you fix the logic structure, not really.  

* ##[MONTHLY_AVG_mapped_DIC_difference_and_lineplots.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//MONTHLY_AVG_mapped_DIC_difference_and_lineplots.ipynb)  
    
    **NOTEBOOK CONTENTS**  
      
        -plots of DIC difference (BR-PI run) by month for averaged fields at 0, 20, 50, 100, 200 m  
        -monthly averaged timeseries of DIC values and their differences for the 4 broad subregions  
        -extra DIC by region and depth on one plot  
      
    **Leading questions:**  
      
    **Does the magnitude extra surface DIC vary much by month, depth, & region?**  
          
    **By month:**  
          
        -no interesting pattern in surface, possibly increase in N./C. strait at depth in summer  
          
    **By depth:**  
        -in general, accross all regions, deeper regions have less extra DIC  
        -surface extra DIC is very blotchy ('stochastic')  
       
    **By region:**  
        -N./C. strait have marginally more extra DIC (and more variable) DIC than Haro,   
        - JdF has considerably less extra DIC than the above 3  
          
    **Do we reproduce coastal high-carbon signal, seasonal surface drawdown?**  
      
        -YES! seasonal surface drawdown more pronounced in inner strait than JdF  
        - at depth (50m, 200m, see beautiful upwelling signal)  

* ##[Monthly_avg_JAN_JUL_PIBR_oma_Sensitivity.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//Monthly_avg_JAN_JUL_PIBR_oma_Sensitivity.ipynb)  
    
    **NOTEBOOK CONTENTS**  
      
        For January and July (averaged fields):  
            plots of averaged: DIC, DIC-TA, & OmA sens. to changes in DIC  
            2 cases: blanket 10 umol DIC perturbation, and BR-PI perturbation  
      
    **Leading questions:**  
      
    ** Which regions/ times of year have more DIC difference (PI vs BR)?**  
          
        - time of year doesn't matter.  
        - inner strait has more of a DIC difference than JdF (though blotchy (plume))  
          
      
    **Are different regions of the SS/ times of year more sensitive to perturbations in DIC? **  
      
    **Does this track with a large TA-DIC?**  
      
    **When doing a uniform 10 $\mu$ mol perturbation, where (and what time of year) is the change in Omega most pronounced?**  
       

* ##[MONTHLY_AVG_CO2_surface_flux.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//MONTHLY_AVG_CO2_surface_flux.ipynb)  
    
    **NOTEBOOK CONTENTS**  
      
        -mean monthly maps of CO2 air-sea flux for preindustrial (PI)  
        and present-day (BR) run  
        -line plot of monthly means by region  
      
    **Leading questions:**  
      
    **When in the year do we switch from ingassing (positive numbers) to outgassing (CO2 out of system, negative numbers)?**  
      
    **Which regions are important?**  



##License

These notebooks and files are copyright 2013-2019
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
