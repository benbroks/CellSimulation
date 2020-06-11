# CellSimulation

## Input
- param.h 
    - N: Number of Cells.
    - T: Number of Cell Transitions.
    - S: CpG Flip Rate. Used to calculate Methy and Demethy Rates. [e.g. for Bin 0.02 -> Methy = S*0.02 & Demethy = S*(1-0.02)]
    - R: Random Cell Replacement Rate. [i.e. At every transition, each cell has R prob. of regenerating with age = 0]
    - OR: Orderly Replacement Rate. [i.e. At every transition, OR*N cells will regenerate with age = 0]
    - X: Neoplasia Cycle Number. The cycle number at which a random cell becomes cancerous/neoplastic.
    - E: Expansion Rate. At every transition, each neoplastic cell has an E prob. of turning a random cell into a neoplastic one.
    - M: Max Expansion Rate. The number of neoplastic cells will never exceed M*N.
    - A: Random Seed. Keep this the same to make experiments repeatable. Set it to -1 if you want a consistently random experiment.
    - P: Statistics Return Regularity. Set P to frequency at which statistics should be returned. Set P to -1 if statistics should just print at the end. [e.g. P=100 will return statsX.csv,matrix_X_1.csv,matrix_X_2.csv for X=100,200,300,etc.]
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
