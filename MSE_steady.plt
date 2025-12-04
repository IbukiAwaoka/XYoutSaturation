# ====================================================================
# gnuplot script to plot MSE(t) from simulation data
# ====================================================================

# Set plot options
set logscale y 
set yrange [1e-9:1]
set xlabel font "Arial,16" "S2/S1"
set ylabel font "Arial,16" "Steady State MSE"
set tics font "Arial,16"
set key right top
set size ratio 0.75
set grid

# Set the terminal to output a high-quality PNG file (optional)
# set terminal pngcairo enhanced font "Arial,12" size 800,600
# set output "MSE_steady_state.png"

# Plot the simulation data from output folder
plot [0:3] \
    'output/dyhouwa1intS2loop, tEnd=300, mu=0.1, S1=0.3, S2=, xi=0.1.txt' using 1:2 with lines linetype 1 lw 1 title "Theory S1=0.3", \
    'output/dyhouwa1intS2loop, tEnd=300, mu=0.1, S1=0.5, S2=, xi=0.1.txt' using 1:2 with lines linetype 2 lw 1 title "S1=0.5", \
    'output/dyhouwa1intS2loop, tEnd=300, mu=0.1, S1=1, S2=, xi=0.1.txt' using 1:2 with lines linetype 3 lw 1 title "S1=1", \
    'output/dyhouwa1intS2loop, tEnd=300, mu=0.1, S1=2, S2=, xi=0.1.txt' using 1:2 with lines linetype 4 lw 1 title "S1=2", \
    'output/dyhouwaerrorbar(N=200, S1=0.3, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 1 pt 7 ps 0.5 title "Simulation S1=0.3", \
    'output/dyhouwaerrorbar(N=200, S1=0.5, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 2 pt 7 ps 0.5 title "S1=0.5", \
    'output/dyhouwaerrorbar(N=200, S1=1, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 3 pt 7 ps 0.5 title "S1=1", \
    'output/dyhouwaerrorbar(N=200, S1=2, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 4 pt 7 ps 0.5 title "S1=2"

# Uncomment to save as PNG
# set output