# Project 25: Wang-Landau simulations and Lee-Kosterlitz method to study discontinuous phase transitions

This project and all the code was implemented by João Pedro Ferreira Biu, for the final project of the course SI2530 Computational Physics at KTH.

The code of the project is strtuctured in the following way:
- The folder `WangLandau/` contains the code needed to run the Wang-Landau simulations and also a Jupyter Notebook with the implementation of Lee-Kosterlitz method to study phase transitions.
- The folder `WangLandau/` is structured in the following manner:
- `main/`: contains the main `.cpp` code to run the Wang-Landau simulations `FirstOrderPT.cpp`.
- `src/`: contains the header (`.h`) and source files (`.cpp`) of the lattice objects and Wang-Landau simulations.
- `lib/`: contains the library files `.a` of the project.
- `bin/`: contains the executable files of the project.
- `results/`: contains the results of the simulations.
- `LeeKosterlitz.ipynb`: Jupyter Notebook with the implementation of the Lee-Kosterlitz method to study phase transitions.

To run this code:
- run `make` in the `WangLandau/` folder to compile the code and the libray.
- run `./bin/FirstOrderPT` to run the Wang-Landau simulations, the results are save in the folder `results/`.
- run the Jupyter Notebook `LeeKosterlitz.ipynb` to study the phase transitions, changing manually the value of the system dimension L.
- run `make clean` to clean the project.

Finally, the folder `images/` contains the images used in the report, `joaobiu_finalproject_si2530.pdf`.
