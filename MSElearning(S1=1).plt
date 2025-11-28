# ====================================================================
# gnuplot script to plot MSE(t) from simulation data
# ====================================================================

# Set plot options
#set logscale y 
set yrange [0:1]
set xlabel font "Arial,16" "t = n/N"
set ylabel font "Arial,16" "MSE"
set tics font "Arial,16"
set key right top
set size ratio 0.75
set grid

# Set the terminal to output a high-quality PNG file (optional)
# set terminal pngcairo enhanced font "Arial,12" size 800,600
# set output "MSE_learning_curves.png"

# Plot the simulation data from output folder
plot [0:100] \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=0.3,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with linespoints linetype 2 pt 2 title "S2=0.3", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=0.5,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with linespoints linetype 3 pt 2 title "S2=0.5", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=1,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with linespoints linetype 4 pt 2 title "S2=1", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=2,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with linespoints linetype 5 pt 2 title "S2=2"

# Uncomment to save as PNG
# set output