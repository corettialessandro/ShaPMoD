SHAKE POLARIZATION MOLECULAR DYNAMICS (ShaPMoD) - VER. 0.2

The code simulates a system of N-interacting ions taking into account the polarization effects. The polarization is modeled through a Shell-Model and the positions of the shells are computed via Mass-Zero Constrained Dynamics.


*Installation and executing:

It will be assumed that all the files of the program have been moved in a folder "/ShaPMoD". The source of the code (*.c files) is located in the folder "/ShaPMoD/source". The corresponding headers are located in the folder "/ShaPMoD/source/headers".

To compile the program just open a Terminal window, move in the folder "/ShaPMoD" and type "make". If the procedure succeeds without errors the executable "ShaPMoD" should appear in the folder "/ShaPMoD". To launch the program, type in the Terminal window "./ShaPMoD".


*Options for launching:

The following options are available for execution:

-i: allows to specify the path of the input file. If not given it is implicitly assumed the file "/ShaPMoD/input/input.txt" as input file.

-o: allows to specify the path of the folder for the output files. If not given it is implicitly assumed the folder "/ShaPMoD/output/" as output folder. If this folder does not exist it is automatically created. WARNING: If the folder already exists it will be cleared by the program and all the file in it will be removed. Anyway a confirmation is asked by the program during execution.

-r: specifies the name of the run. This means that all the input and output folder will be named as <runname>_input.txt <runname>_output.txt 

-d: debug mode.

-v: verbose mode. To be used in polarization mode to monitor the convergence of the algorithm.


*Input file:

The input file should be located in the folder "/ShaPMoD/input/". The path can anyway be specified using the option -i. An example of input file is given in the folder "/ShaPMoD/input/". WARNING: In this version only binary mixtures are implemented. 

It is possible to specify to modes for the execution: Standard mode 'S' and Polarization mode 'P'. The number of molecules should always be half of a perfect cube (4 = 2^3/2, 108 = 6^3/2). This version only implements a simple cubic starting configuration. It is not possible to restart a simulation from a previous final configuration (coming soon) so each time the program restarts from a simple cubic configuration.

The units are reduced units (RU). The parameters B and LAMBDA, set to one in the input file, have been used as the unit of energy and length respectively (see Jacucci et al.). Masses and charges are expressed in amu and units of electron charge.

There is the possibility to specify a starting mass for the Shell, but this feature is not used in this version.