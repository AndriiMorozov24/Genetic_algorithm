# Genetic_algorithm
 Genetic algorithm for the flow-shop scheduling problem
andrii.moroozv23@gmail.com

Enter name of your input file. (3rd line) To start calculations, provide data (4th line), it must be .csv format. Rows must contain units of time for particular task on particular machine. Columns contain machine number. Each task should go from machine 1 to machine N (depends how many you set). Task on particular machine can be done only if machine is free. Fitness function (27th line) take care of it (makespan).

You can regulate speed of the program changing number of persons in the population by changing variable "POP". (11th line) Also you can change number of iterations by changing "k" variable. (12th line) and set probability of mutation and recombination (13-14 line). In (204-207 lines) you can set at which point probability of mutation / recombination will increase / decrease or you can even deny such option by commenting this block of code, feel free to manipulate with this via your needs.

Also Roulette selection method was implemented, feel free to use it by uncommenting (67-95 lines). If you will do this be sure to change Genetic function (150-70 lines) by changing "TournamentSelection" into "RouletteSeletion".

If you want to use an OX crossover instead of PMX, uncomment (158-160 lines) and comment (155-157 lines).

At the end will be created file called "output.txt" with results. Also repository contains couple .csv data files and example outputs.
