# Codes for solving spatial and non-spatial positive plant-soil feedback model
## Manusciprt titled "Niche theory for positive plant-soil feedbacks" has the details of the model and results
## Authors: Athmanathan Senthilnathan, Rafael D'Andrea

## Files
* *parGen.sh*: Bash script to generate parameter files
* *main.m*: MATLAB function to numerically solve for a single parameter combination
* *fig1.nb*: Mathematica notebook to produce figures 1 and S1
* *fig2a.m*: MATLAB script to produce figure 2a
* *fig2b.m*: MATLAB script to produce figure 2b
* *fig3a.m*: MATLAB script to produce figure 3a
* *fig3b.m*: MATLAB script to produce figure 3b
* *myTxtFmt.m* : MATLAB script to change font size and interpreter
* *printPdf.m*: MATLAB script to print plot as PDF

## Instructions to find numerical solutions
1. Create *output* directory
2. Edit *parGen.sh* and run the bash script to create batch sub-directory under *output* and generate parameter files
> ./parGen.sh batch\_name
3. Run *main.m* in MATLAB to numerically solve for one parameter combination. Code requires *batch_name*, *run_number* and *visualization status* (0 or 1 to indicate whether output should be plotted in runtime.
