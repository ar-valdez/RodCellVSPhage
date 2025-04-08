# RodCellVSPhage

Forked from https://github.com/huiprobable/CellsMD3D

# Compilation
A makefile has been added. You must create directory for compilation objects.
The executable file is located in "/run_sim/"

# Execution
Inside "/run_sim/" open a terminal and type the following statement

./Main T7Plaque.txt 32 ResultsFolder 3 monoLayer/Starter.dat > ResultsFolder/Simul.log

# Explanation
1. The command line argument T7Plaque.txt is the input txt file, contains ALL the custom input parameters. 
2. The command line argument 32 means we run this code using 32 Threads.
3. The command line argument "ResultsFolder" is the directory to output the data.
4. The command line argument "monoLayer/Starter.dat" (optional) initializes the test with a previous known solution.


---
# Contact

Andres R. Valdez "arvaldez@psu.edu"

---
# How to cite
If you find this repository useful, cite the following paper:

@article{Valdez2025,
title={Biomechanical modeling of spatiotemporal bacteria-phage competition},
author={Valdez, Andr\'es R and Sun, Paul and Weiss, Howie and Aranson, Igor S },
journal={Nature Communications Physics},
doi={https://doi.org/10.1038/s42005-025-02078-1},
year={2025}}

