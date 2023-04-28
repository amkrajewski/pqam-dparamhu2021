install.packages("locfit")  # You only need to install the locfit package once to run the rest parts of the code.

dir=getwd()
setwd(dir)

source("HEA_pred.R")

composition=c(2,0,0,0,1,0,1,1,0,0)
###composition###########################################################################################
#Above you input the alloy composition of interest. The input is enssentially a composition array in 
#the order of c(Ti,Zr,Hf,V,Nb,Ta,Mo,W,Re,Ru). For example, for an composition of Ti2NbMoW, its corresponding
#input should be composition=c(2,0,0,0,1,0,1,1,0,0). Please do not change the length of this composition
#array. 
########################################################################################################

prediction=HEA_pred(composition)
###prediction###########################################################################################
#Once the composition anrray is correctly set up, please execute every line in this file. The predicted
#GSF energy, Surface energy, and D parameter will show in the Console window below. The unit of GSF and
#Surface energy is J/m^2. 
#
#Please note that the predictions made here are based on the SL models trained with all the binary, 
#ternary, and quateranry data listed in Table 1 in the manuscript. 
########################################################################################################
