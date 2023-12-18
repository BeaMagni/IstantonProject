import Exact_solution
import Monte_Carlo
import MC_cooling
import Free_energy_as
import Instanton_hist
import Instanton_action
import Plots
from tabulate import tabulate

headers = ["NUMBER", "CODE"]
rows = [[1, "Exact solution"], [2, "Monte Carlo solution"], [3, "Monte Carlo Cooling"], [4, "Free energy adiabatic switching"], [5, "Instantons histogram"], [6, "Instantons action"], [7, "Plots"]]

print(tabulate(rows, headers=headers))

num = input('Select the code to run: ')
if num == 1:
    Exact_solution.
if num == 2:
    Monte_Carlo.
if num == 3:
    MC_cooling.    
if num == 4:
    Free_energy_as.
if num == 5:
    Instanton_hist.
if num == 6:
    Instanton_action.
if num == 7:
    head = ["NUMBER", "PLOT"]
    lines = [[1, "Anharmonic potential"], [2, "Energy eigenvalues"], [3, "Ground state wavefunction"], [4, "Monte Carlo correlators"], [5, "Monte Carlo log correlatos"], [6, "Euclidean path"], 
             [7, "Monte Carlo cooled correlators"], [8, "Monte Carlo cooled log correlators"], [9, "Free energy"], [10, "Instantons distribution"], [11, "Instantons density"], [12, "Instantons action"]]

    print(tabulate(lines, headers=head))
    plot = imput('Select the plot to show: ')
    if plot == 1:
      Plots.
    if plot == 2:
      Plots.
    if plot == 3:
      Plots.
    if plot == 4:
      Plots.
    if plot == 5:
      Plots.
    if plot == 6:
      Plots.
    if plot == 7:
      Plots.
    if plot == 8:
      Plots.
    if plot == 9:
      Plots.
    if plot == 10:
      Plots.
    if plot == 11:
      Plots.
    if plot == 12:
      Plots.



           
