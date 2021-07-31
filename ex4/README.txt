The folder contains:

"mulitplication.f" is the Fortran script
"mult_proc.py" is the script used for the simulation
"run.py" was the previous script, before writing the parallelization.
"multiplicationN.txt" are the executables called by "mult_proc.py", where N=1,2,3,4.
"dimensionsN.txt" are the txt file in which the dimensions of the matrix are saved by "mult_proc.py" and read by the executables.
"std_mult.txt", "clmn_mult.txt" and "built_in.txtx" are the files containing the simulations, for the three multiplicaation methods respectively.
"fit.plt" is the gnuplot script for fitting and plotting. It must be load inside this directory with either the command "gnuplot 'fit.plt'" from terminal, or "load 'fit.plt'" from gnuplot.
The .png files are the graphs of simulation.
THe folder old_file contains a previous simulation
logfile.log is the log file for the gnuplot fits.
Ex4_Bassi_Report.pdf is the report.

From the next time maybe I will directly build a program to be installed, along the instructions to install it.
