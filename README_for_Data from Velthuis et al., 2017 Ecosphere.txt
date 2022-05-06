
Datafiles:

biweekly_data.csv: Data on Chl-a (Fig. 1), DIN & DIP (Fig. 4a-b), NO3:NH4 ratio (Supplement 3), phytoplankton size classes (Fig. 2) and chytrid prevalence (Fig. 5a)
biweekly_data metadata.csv: metadata file belonging to biweekly_data.csv

weekly_data.csv: Data on seston C:N, C:P and N:P ratios (Fig. 2), light availability (I50) and dissolved CO2 (Supplement 3)
weekly_data metadata.csv: metadata file belonging to weekly_data.csv

LME_zooplankton_data.csv: Data on abundances of rotifers, cladocera and copepods (Fig. 5b-d)
LME_zooplankton_data metadata.csv: metadata file belonging to LME_zooplankton_data.csv

Si_data.csv: Data on dissolved silica concentrations (Fig. 4c)
Si_data metadata.csv: metadata file belonging to Si_data.csv

Temperature_data : Data on the temperature over the course of the experiment (Supplement 1)
Temperature_data metadata.csv: metadata file belonging to Temperature_data.csv

Light_data: Data on the hours of light over the course of the experiment (Supplement 1)
Light_data metadata.csv: metadata file belonging to Light_data.csv

Phytoplankton_microscopy_data.xls: Data on the relative abundance of phytoplankton groups (microscopic determination) (Supplement 2)
Phytoplankton_microscopy_data metadata.csv: metadata file belonging to Phytoplankton_microscopy_data.csv

periphyton_data.csv: Data on periphyton abundance (Supplement 4)
periphyton_data metadata.csv: metadata file belonging to periphyton_data.csv

epipelon_data.csv: Data on the abundance of benthic algae (Supplement 4)
epipelon_data metadata.csv: metadata file belonging to epipelon_data.csv

Floating_and_filamentousalgae_data.csv: Data on the percent volume infested (PVI) of filamentous algae and the cover of floating algae (Supplement 4)
Floating_and_filamentousalgae_data metadata.csv: metadata file belonging to Floating_and_filamentousalgae_data.csv


--- More info on all of the above datafiles can be found in their specific metadata file --- The datafiles below contain the same data but reorganized for their respective R-analysis ---


weibull_data_chla: Chla data from the biweekly_data file, reorganized in such a way that it can be used as input for the weibull analysis (Table 2)

weibull_data_chytrids: Chytrid prevalence data from the biweekly_data file, reorganized in such a way that it can be used as input for the weibull analysis (Table 2)

weibull_data_cladocera: Cladocera abundance data from the weekly_data file, reorganized in such a way that it can be used as input for the weibull analysis (Table 2)

weibull_data_copepods: Copepod abundance data from the weekly_data file, reorganized in such a way that it can be used as input for the weibull analysis (Table 2)

weibull_data_rotifers: Rotifer abundance data from the weekly_data file, reorganized in such a way that it can be used as input for the weibull analysis (Table 2)

data_variable_contribution_to_chla: datafile used as input for the analysis presented in Table 3. Is a combination of the weekly_data and biweekly_data file 


R scripts:


R script LME: produces the output of the linear mixed effect models for time serie analysis (Table 1)

R script Weibull fits: produces the output for the weibull analysis (Table 2)

R script variable contribution to chlorophyll a dynamics: produces the output presented in Table 3


Folders:


Weibull area under the curve calculations: Contains calculations for total abundance as presented in Table 2. Specific files in this folder match to the results of the different variables in Table 2.


