# CellSimulation

## Summary ##
This project is used to simulate the life cycle of _N_ cells over _T_ transitions. Each cell is represented as a collection of (27634) CpG sites. Each CpG site represents two binary bits (taking on values of _0_ or _1_), which is concretely stored as an integer with value _0_, _1_, or _2_.

At each discrete transition, a variety of alterations can be made to each cell. Cells can die and become randomly replaced by a copy of another cell or a brand-new cell can pop up in its old place. Cells can mutate and only change at a few of their CpG sites. See the attached input and output sections for more information on the configurability of the tool.

Information about the cells is stored in `csv` format. Each row represents a single cell. Each column represents a unique CpG site. (There are _N_ rows and 27634 columns) The final state of the cell "colony" is displayed via two CSVs (as not to exceed the Excel column limit). The first file contains the first 27634/2 = 13817 columns and _N_ rows. The second file follows the same dimensions, but the data from the latter half of columns. We collect intermediate and final statistics about the cells. These are contained within a separate statistics output file. 


## Instructions
- Download or clone this project! Unzip it and navigate inside the directory via your Terminal.
- Verify that your machine is able to compile C++ code. You need `g++`.
- Run `make` to compile the project.
- Debug/prerequisite instructions (for mac):
    - Verify that `xcode` is downloaded on your machine. Verify that `command line developer tools` are downloaded on your machine.
    - Are you seeing `xcrun: error: invalid active developer path`? Follow instructions linked [here](https://apple.stackexchange.com/questions/254380/why-am-i-getting-an-invalid-active-developer-path-when-attempting-to-use-git-a).
    - Are you seeing `fatal error: 'filesystem' file not found`? Follow instructions linked [here](https://stackoverflow.com/questions/39231363/fatal-error-filesystem-no-such-file-or-directory/39231488).
- Run `./simulation` (after successfully compiling) to actually run the cell simulator. Expect the simulation to take a number of hours if transitions or cells number in the tens of thousands. 
- Nuclear Option if nothing works... (this will require a bit more time but has consistent results):
    - Follow these [Virtual Machine Download Instructions](http://bits.usc.edu/cs104/installing-course-vm.html). This creates a virtual machine and a consistent environment in which to compile/run C++.
    - Launch the VM and follow this repo's (not the vm download's) instructions from inside it.

## Cleaning Up
- Run `make clean` to purge the executable and associated files.

## Input
- param.h 
    - N: Number of Cells.
    - T: Number of Cell Transitions.
    - B: Mode of CpG Site Flip Rates [Select 0 for an Exponential Relationship (peak at Bin 0.5) and 1 for Manual Entry.
    - SMin: Lower bound of CpG Site Flip Rate.
    - SMax: Upper bound of CpG Site Flip Rate.
    - R: Random Cell Replacement Rate. [i.e. At every transition, each cell has R prob. of regenerating with age = 0]
    - OR: Orderly Replacement Rate. [i.e. At every transition, OR*N cells will regenerate with age = 0]
    - X: Neoplasia Cycle Number. The cycle number at which a random cell becomes cancerous/neoplastic.
    - E: Expansion Rate. At every transition, each neoplastic cell has an E prob. of turning a random cell into a neoplastic one.
    - M: Max Expansion Rate. The number of neoplastic cells will never exceed M*N.
    - C: Transition Survival Rate. C*N cells will survive a transition. (1-C) * N cells will be replaced by random living cells.
    - A: Random Seed. Keep this the same to make experiments repeatable. Set it to -1 if you want a consistently random experiment.
    - P: Statistics Return Regularity. Set P to frequency at which statistics should be returned. Set P to -1 if statistics should just print at the end. [e.g. P=100 will return stats_X.csv for X=100,200,300,etc.]
    - init_table_fp: Gives file path to read from to obtain bin sizes.
    - matrix_output_fp: Prefix to matrix output filepaths. [e.g. Given "matrix", outputs will be "matrix1.csv" and "matrix2.csv"]
    - stats_output_fp: Prefix to simulation statistics filepath. [e.g. Given "stats", output will be "stats.csv"]
- init_table_fp
    - Pulls information from the Integer Bins column. (Doesn't use any other columns)

## Output
- matrix_output_fp
    - Outputs final state of simulation. This includes cell age AND CpG values for all ~27k sites. Output is split across two files to support MS Excel viewing.
- stats_output_fp
    - Outputs average values across all cells for each CpG site. ~27k dimensional vector.
    - Outputs average values across all neoplastic cells for each CpG site. ~27k dimensional vector.
    - Outputs a variety of single value statistics.
        - Mean CpG value across all sites and cells. 
        - Variance of CpG value across aforementioned average vector.
        - Mean Cell Age.
        - Number of Neoplastic Cells.
        - Mean CpG value across all sites in neoplastic cells.
        - Variance of CpG value across average vector for neoplastic cells.
