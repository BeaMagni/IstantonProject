# Instantons and Monte Carlo Methods in Quantum Mechanics

In this repository, you will find the Python codes to recreate some of the results presented in the "Instantons and Monte Carlo Methods in Quantum Mechanics" paper by T. Schäfer (https://arxiv.org/abs/hep-lat/0411010).

What you will find:
- Main.py: this is the main program, you should run this one to obtain all the possible results, following the instructions;
- General_functions.py: this file contains most of the functions used in the various programs;
- Plots.py: this file contains all the functions related to the plots;
- lastly, there are all the other files called in the Main.py file, containing specific functions and/or the main part, giving the corresponding results.

To run the codes, you simply need to download the files and use a Python compiler. 
It's important that all of the .py files are in the same folder on your computer.

We included a function that, if it doesn't already exists, creates the directory (in the same location where the file .py are stored) containing the various .txt files with the results of the computations. The main folder, called instanton_project, will be generated in the created directory: inside this, there will be created all the other folders, divided by the specific results you want to obtain (for example, you will have different folders for the exact solution and the Monte Carlo simulations results). 
