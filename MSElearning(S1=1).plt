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
    'output/dyhouwa2MSE, tEnd=100, mu=0.1, S1=1, S2=0.3.txt' using 1:4 with lines linetype 1 lw 1 title "Theory S2=0.3", \
    'output/dyhouwa2MSE, tEnd=100, mu=0.1, S1=1, S2=0.5.txt' using 1:4 with lines linetype 2 lw 1 title "S2=0.5", \
    'output/dyhouwa2MSE, tEnd=100, mu=0.1, S1=1, S2=1.txt' using 1:4 with lines linetype 3 lw 1 title "S2=1", \
    'output/dyhouwa2MSE, tEnd=100, mu=0.1, S1=1, S2=2.txt' using 1:4 with lines linetype 4 lw 1 title "S2=2", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=0.3,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with points linetype 1 pt 5 ps 0.5 title "Simulation S2=0.3", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=0.5,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with points linetype 2 pt 7 ps 0.5 title "S2=0.5", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=1,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with points linetype 3 pt 9 ps 0.5 title "S2=1", \
    'output/MSE Simulation(N=200,T=100,S1=1,S2=2,mu=0.1,ens=1000,xi=0,a=2).txt' using 1:2 with points linetype 4 pt 11 ps 0.5 title "S2=2"

# Uncomment to save as PNG
# set output