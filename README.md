# hedgerow_metacom
Analysis code for "Proximity of restored hedgerows interacts with local floral diversity and species traits to shape  long-term pollinator metacommunity dynamics". 
Please refer to witeup (writeup/writeup.pdf) for a description of the organization of the repo and instructions on how to use the scripts. The data accompanies the code for reproducibility, but we ask to be contacted (lponisio@gmail.com and ckremen@zoology.ubc.ca) if authors wish to use the data for a publication.

In our study we examine the metacommunity dynamics of plant-pollinator communities using variety of different methods including 1) occupancy modeling and 2) network analyses.  We are committed to reproducible science and all analytical code will be maintained on github, along with this write up.

The entire analysis is executable from the main.sh file. All of the packages needed to run the analyses are listed in the packages.sh file. All analyses were run using R (version $3.5.1$) and nimble (0.6-12).

Note, there is an incompatibility with Nimble version 7.0 released in February 2019. If this effects you can you either install the older version of Nimble, or install from the branch "avoid-protect-stack-overflow" directly from github. See post at \url{https://groups.google.com/forum/#!topic/nimble-users/k6VMapOfxOk}

Hopefully this bug will be fixed soon and the models will run on any version of NIMBLE.

Navigate to the analysis folder within the github repo (hedgerow\_metacommunity) then the main.sh file can be selected and run (a warning, the occupancy analyses each take several hours on my 2.3 GHz imac pro, so all together they will take quite a while), you could run all of the analyses in the study by running this line in BASH.

This will somewhat helpfully print the results of each analysis and re-create any accompanying figures.

I walk through each the main script for each analysis individually in the write-up.pdf

