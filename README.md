# Genetic Algorithm for the Flow-Shop Scheduling Problem
Author: Andrii Morozov (andrii.morozov23@gmail.com)

# Overview
This project implements a genetic algorithm to optimize the flow-shop scheduling problem. The algorithm processes tasks through multiple machines while minimizing the makespan (total completion time). The solution supports customization of key parameters, allowing users to fine-tune performance based on their needs.

How to Use
1. Input Data:
Provide a .csv file as input.
Rows represent processing times for each task on a specific machine.
Columns represent different machines.
Each task must be processed sequentially from Machine 1 to Machine N.
A task can only be executed if the assigned machine is free.


2. Configurable Parameters:
Population Size: Adjust via the POP variable (line 11).
Number of Iterations: Set using the k variable (line 12).
Mutation & Recombination Probability: Modify these values on lines 13-14.
Dynamic Mutation & Recombination:
Enable/disable probability adjustments at runtime (lines 204-207).
Comment out this section if you prefer a fixed probability.


3. Selection & Crossover Methods:
Roulette Selection:
Enable by uncommenting lines 67-95.
Ensure you modify the Genetic function (lines 150-170) by replacing TournamentSelection with RouletteSelection.
Crossover Methods:
Use OX crossover instead of PMX by uncommenting lines 158-160 and commenting out lines 155-157.


4. Output
After execution, the results will be saved in an output.txt file. The repository also includes sample .csv datasets and example output files for reference.
