from tabulate import tabulate
import Exact_solution 
import Monte_Carlo
import MC_cooling
import Free_energy_as
import Instanton_hist
import Instanton_action_density
import Plots

stop = False
while stop == False: #the table reappears until the exit condition is chosen
    headers = ["NUMBER", "CODE"]
    rows = [[1, "Exact solution"], [2, "Monte Carlo solution"], [3, "Monte Carlo Cooling"], [4, "Free energy adiabatic switching"], [5, "Instantons histogram"], [6, "Instantons action"], [7, "Plots"], [8, "Stop running"]]

    print(tabulate(rows, headers=headers)) #creation of the table showing the possible programs one can run

    num = int(input('Select the code to run: '))

    if num == 1:
        Exact_solution.main()
    elif num == 2:
        Monte_Carlo.main()
    elif num == 3:
        MC_cooling.main()
    elif num == 4:
        Free_energy_as.main()
    elif num == 5:
        Instanton_hist.main()
    elif num == 6:
        Instanton_action_density.main()
    elif num == 7:
        stop_plot = False
        while stop_plot == False:
            head = ["NUMBER", "PLOT"]
            lines = [[1, "Anharmonic potential"], [2, "Energy eigenvalues"], [3, "Ground state wavefunction"], [4, "Monte Carlo correlators"], [5, "Monte Carlo log correlatos"], [6, "Euclidean path"], 
                     [7, "Monte Carlo cooled correlators"], [8, "Monte Carlo cooled log correlators"], [9, "Free energy"], [10, "Instantons distribution"], [11, "Instantons density"], [12, "Instantons action"],
                     [13, "Stop running plots"]]

            print(tabulate(lines, headers=head)) #creation of the table showing the possible plots one can obtain
            
            plot = int(input('Select the plot to show: '))
            if plot == 1:
                Plots.potential_eigenvalues()
            elif plot == 2:
                Plots.energy_eigenvalues_variation()
            elif plot == 3:
                Plots.ground_state_hist()
            elif plot == 4:
                Plots.correlation_functions()
            elif plot == 5:
                Plots.log_correlation_functions()
            elif plot == 6:
                Plots.euclidean_path()
            elif plot == 7:
                Plots.correlation_functions_cooling()
            elif plot == 8:
                Plots.log_correlation_functions_cooling()
            elif plot == 9:
                Plots.free_energy()
            elif plot == 10:
                Plots.instanton_distribution()
            elif plot == 11:
                Plots.instanton_density()
            elif plot == 12:
                Plots.instanton_action()
            elif plot == 13:
                stop_plot = True

    elif num == 8:
        stop = True
