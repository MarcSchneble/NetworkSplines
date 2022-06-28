# NetworkSplines
Code and data for reproducing the article

Marc Schneble and Göran Kauermann (2022): "Intensity Estimation on Geometric Networks with Penalized Splines"

This folder contains the following sub-folders and R-files:

Data: Contains all the external data
	- CarCrashesMontgomery_Incidents.csv: Original data about car crashes in Montgomery County downloaded from https://data.montgomerycountymd.gov/Public-Safety/Crash-Reporting-Incidents-Data/bhju-22kf (date of download: 17th November 2020).
	- data_2019_clean.rds: Preprocessed data about parking event in Melbourne. Orginal data downloaded from https://data.melbourne.vic.gov.au/Transport/On-street-Car-Parking-Sensor-Data-2019/7pgd-bdf2.
	- EdgesMelbourne.xlsx: Edges of the Melbourne network (vertex IDs of start and end point, own creation)
	- EdgesMontgomery.xlsx: Edges of the Montgomery network (vertex IDs of start and end point) 
	  including internal covariates (own creation)
	- Parking_Lots.rds: Preprocessed data about on-street parking lots in Melbourne. Original data downloaded from https://data.melbourne.vic.gov.au/Transport/On-street-Parking-Bay-Sensors/vh2v-4nfs.
	- VerticesMelbourne.xlsx: Coordinates of the vertices of the Melbourne network (own creation)
	- VerticesMontgomery.xlsx: Coordinates of the vertices of the Montgomery network (own creation)

Plots: In this folder, all plots will be saved when the following R-scripts are sourced.

Chicago.R: Runs the analyses in Section 5 and the simulation study in Section 8
Figure1.R: Produces Figure 1
Functions.R: Contains function which are used in other scripts. 
MelbourneParking.R: Runs the analyses in Section 6 
MontegomeryCarCrashes.R: Runs the analyses in Section 7