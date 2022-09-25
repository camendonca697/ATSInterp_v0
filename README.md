# ATSInterp_v0
Programs for manuscript "THIN SHEET AUTOMATIC INTERPRETATION USING REGULARIZED SECOND DERIVATIVES"

There are two main files to reproduce all figures in the manuscript:  
    RunMe_Synthetics tests the proposed procedure to a set of 10 models with multiple prismatic bodies
    RunMe_RealData applies the proposed method to a profile interpreted by Cavalcante et al. Journal of Hydrology, 2020, DOI:10.1016/j.jhydrol.2020.125079.


TESTS WITH SYNTHETIC DATA   
RunMe_Synthetics.m
aux01_TestingModels.m configure a set of 10 prismatic models for testing
aux02_PrismPlotting.m graphical representation of the testing models
multiprism.m model response evaluation

REAL DARA APPLICATION
RunMe_RealData.m    
dd.dat interval between stations
x0.dat vector with distance along the profile
ft.dat magnetic anomaly along the profile
tt.dat ampltude of the magnetic anomaly 
JH_model.dat model parameters as in Cavalcante et al (2020)       
JH_DataFitting.dat  data fitting as evaluated by Cavalcante et al. (2020)    
aux01_ModelJH.m make the picture of the prismatic models from Cavalcante et al. (2020)  
aux03_ZaPlotting.m picture with automatic solutions
ZaResults.xls table with model solutions 
          
rdiff.m  evaluate Tikhonov regularized derivatives 
from Wagner, J., 2021, Regularized numerical 427 differentiation. Accessed on July 17.
https://www.mathworks.com/matlabcentral/fileexchange/74165       
     
              
               
