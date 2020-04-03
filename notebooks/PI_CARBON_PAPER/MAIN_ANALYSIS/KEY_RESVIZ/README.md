The Jupyter Notebooks in this directory are for sharing of Python code
techniques and notes about model results analysis code.
They were developed by Tereza Jarnikova.

The links below are to static renderings of the notebooks via
[nbviewer.jupyter.org](http://nbviewer.jupyter.org/).
Descriptions under the links below are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[Carbon_preindustrial_paper_overview.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//Carbon_preindustrial_paper_overview.ipynb)  
    
    **The inorganic carbon balance and fate of anthropogenic carbon in a temperate fjord system**  
      
    We are trying to write a paper investigating the fate of anthropogenic carbon in a fjord system, using the SalishSeaCast Model, with custom-built carbonate chemistry module, as our tool. The basic idea is to run the SKØG model in 3 configurations - present-day co2 forcing; preindustrial co2 forcing + present-day boundary conditions; and preindustrial co2 forcing + preindustrial boundary conditions.   
      
    Then we try to see:   
      
    a) where the extra carbon is   
      
    b) how it gets transported thru lateral (b1) and air-sea (b2) boundary and how that balance  changes in the 3 scenarios.  
      
    c) the effect this has on ecologically-meaningful quantities (eg $\Omega_A$, pH)    
      
    This notebook tries to keep track of where the model code is, the run files, the results, and all the different bits of analysis are for a) Tereza, who could stand to be more organized, and b) advisors/collaborators - Debby, Susan, et al. The idea is to put links to relatively tidy notebooks that distill what has been done. The paper is being written on Overleaf.     
      
    Outline of doc:  
      
    * **1a** Model Code/ Config  
    * **1b** Init atmospheric CO2  
    * **1c** PI boundary conditions  
    * **2** Model Runs Description  
    * **3** Model Evaluations  
    * **4** Results  
      
    QQ means to fill in  
    TD means... todo  

* ##[ALL3runs_MASSBAL_FATE_anthropogenic_carbon.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//ALL3runs_MASSBAL_FATE_anthropogenic_carbon.ipynb)  
    
* ##[ALL3runs_MASSBAL_BR_LA_PI.ipynb](http://nbviewer.jupyter.org/urls/bitbucket.org/tjarnikova/analysis-tereza/raw/tip/notebooks/carbon_dev/PI_CARBON_PAPER/MAIN_ANALYSIS/CLEAN/KEY_RESVIZ//ALL3runs_MASSBAL_BR_LA_PI.ipynb)  
    

##License

These notebooks and files are copyright 2013-2020
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
