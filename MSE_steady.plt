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
    'output/Steady State MSE(N=200,T=5000,S1=0.3,mu=0.1,ens=100,xi=0).txt' using 1:2 with linespoints linetype 2 pt 7 title "Simulation (S1=0.3, μ=0.1)"

# Uncomment to save as PNG
# set output