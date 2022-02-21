bootstrap_distpara_obtain_documented_Jan2018.m is both a script and a documentation of the bootstrap exercise. 
The goal is that with this file, and the dependencies (made copies of their source in folder Masao and saved in "copies of dependencies that has not been updated"), one could redo the bootstrap as long as he knows how to run matlab on a server. 

The beginning sessions of bootstrap_distpara_obtain_documented_Jan2018.m also achieves code that translate optimization results to parameters in SummaryXXXXXX.xlsx.

As is mentioned in bootstrap_distpara_obtain_documented_Jan2018.m, it first optimize some parameters, which are inputs to the likelihood function, and then use the likelihood function to compute some other parameters. Currently, there is the issue of duplicated book index as discussed in email thread "Bootstrap" in Nov/Dec 2017. Hongkai has updated the output of the second step using correct indexing and old outputs of the first step. The hope is that the issue mainly affects the second step. 

Reflecting this partial fix, files in this folder updates some previous files:

bootstrap_distpara_obtain_documented_Jan2018.m replaces bootstrap_books.m (mode 1) and bootstrap_dispara_obtain.m (mode 2)

bootstrap_welfare_Dec2017.csv replaces bootstrap_welfare.csv (bootstrap estimates of parameters in format of likelihood function input and bootstrap estimates of welfare)

bootstrap_estimates_Dec2017.csv replaces bootstrap_estimates2017_1105.csv (bootstrap estimates of parameters in same format as SummaryXXXXX.xlsx)

bootstrap_welfare_estimates_Dec2017.csv replaces bootstrap_welfare_estimates2017_1105.csv (bootstrap estimates of welfare in same format as SummaryXXXXX.xlsx)

SummaryDec2017.xlsx replaces Summary201609 (true estimates, s.e. and c.i. of bootstrap)

bootstrap_distpara_obtain_documented_Jan2018.m documentation still refers to names of previous files in case of replacement to make it easier to follow previous emails communications.