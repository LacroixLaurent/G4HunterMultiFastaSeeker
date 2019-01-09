# G4Hunter Multifasta Seeker
#### See [G4HunterApps, Lacroix. Bioinformatics 2018](https://doi.org/10.1093/bioinformatics/bty951).

#### Shiny App related to G4Hunter published in [Bedrat _et al._ NAR 2016][paper ref].  
Supplementary Data can be downloaded from [NCBI](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770238/bin/supp_44_4_1746__index.html) or [Github](https://github.com/LacroixLaurent/G4HunterPaperGit).  


> ##### The top part of this page allows you to compute the G4Hunter score of a single sequence. This main App (G4Hunter Seeker) identifies DNA or RNA regions in a longer sequence for which the G4Hunter score is above the chosen threshold in windows of the selected size. Please cite _"Bedrat A, Lacroix L, & Mergny JL (2016) Re-evaluation of G-quadruplex propensity with G4Hunter. Nucleic Acids Res 44(4):1746-1759."_, when reporting results obtained with this App.

##### The app requires the following packages:
* Biostrings
* GenomicRanges
* shiny

##### To run the app
1- download the project from github.  
2- set the directory where you extract the project as your working directory in R by using the command **setwd('PATH_TO_THE_PROJECT')**.  
3- install the required packages by running the **install-packages.r** script.  
4- in the R-console, type **runApp()**.  
5- a browser page should open with the app.  
6- for the next time you want to run the app, you can just go your R-console and type **shiny::runApp('PATH_TO_THE_PROJECT')**. Thus if you have installed this app in a directory named work under your home directory, you should type **shiny::runApp('~/work/G4HunterMultiFastaSeeker/')**.  


### Quick G4Hunter score
Just type your sequence in the box and you get the G4Hunter score below.  
Spaces are automaticaly removed, lower or upper cases are accepted.  
Characters orther than **G** or **C** are kept and counted as **A**,**T** or **U**.

###  G4Hunter Multifasta Seeker
Choose your multifasta file  (**Fasta File entry**).

The file to download has to be a **DNA fasta** file that does not exceed the size limit imposed by **Shiny** (default is 5Mb, see below how to change it).  
Invalide DNA letters are skipped and give rise to a warning.  
Please note that if the number of entries in your multifasta file is big (more than 1000), the app will become **very slow**. You might need to use **R** outiside of the **Shiny** interface. To avoid unwanted long computation time, you need to click on the button **Please click here to start the computation** to start the G4Hunter process after uploading your file and checking the number of fasta entries.  

The **Threshold** and **Window size** determine the parameters for the sequence search as described in the [publication][paper ref].  
The higher the threshold, the more stringent the search: fewer G4 motifs wil be found, but these will be the most stable/likely ones.  

The **Report sequences** option adds the nucleotide sequences in the output.  
The **Report G_sequences** option changes sequences with a negative score (C-rich sequences) into their reverse complement. Thus the output reports only G-rich sequences.

The **hits** report the number of sequences retrieved that match the settings.  


The output table contains a unique name for each hit (**hitnames**), the sequence name (**seqnames**) corresponding to the name of the fasta entry where this hit has been found, the **start**, **end** and **width** of the _refined_ sequences that meet the search criteria.  
The **strand** is **+** if the proposed G4 forming sequence is in the Input Sequence and this is set to **-** if the G4 forming sequence in on the reverse complementary strand.  
The **score** is the G4Hunter score of the _refined_ sequence and **max_score** is the highest score in absolute value in a window of the chosen **window size** for the sequence.  
**Threshold** and **window** are respectively the **Threshold** and **Window size** used for the search.  
The **sequence** corresponds to the _refined_ sequence in the **Input sequence**. This field is sensitive to the **Report G-sequences** option. This field is not present if the **Report sequences** option is not selected.  

> ##### Please note that the procedure extracts sequences that have a G4Hunter score above the threshold (in absolute value) in a window, fuses the overlapping sequences and then _refines_ theses sequences by removing bases at the extremities that are not G for sequences with a positive score (or C the negative ones). It also looks at the first neigboring base and adds it to the sequence if it is a G for sequences with a positive score (C for sequences with a negative score).  
> ##### Please see the [publication][paper ref] and Figure S1B for more details

The output can be exported to a text file that can be directly opened with _Microsoft Excel_.

--------------------------------------------------------------------------
As this is a simple **Shiny** app, there is a limit to the size you can upload (5Mb).  
If necessary (but might not be good for the host), the size limit can be increased  to XX Mb by adding the following code at the top of the server.R file (_adapted [from stackoverflow](http://stackoverflow.com/questions/18037737/how-to-change-maximum-upload-size-exceeded-restriction-in-shiny-and-save-user)_).  
```{r}
options(shiny.maxRequestSize=XX*1024^2)
```

--------------------------------------------------------------------------
_This project is independent of the [G4-Hunter python repository](https://github.com/AnimaTardeb/G4-Hunter) even if both project are related to the same [publication][paper ref]._


[paper ref]:http://doi.org/10.1093/nar/gkw006
