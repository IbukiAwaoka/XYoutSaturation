# ====================================================================
# gnuplot script to plot MSE(t) from simulation data
# ====================================================================

# Set plot options
#set title "Steady State MSE vs S2/S1 Ratio"
set logscale y 
set format y "10^{%T}"
set yrange [1e-9:1]
set xlabel font "Arial,20" "S_2/S_1"
set ylabel font "Meiryo,20" "準定常MSE"
set tics font "Arial,12"
set key right bottom font "Arial,12"
set size ratio 0.75
#set grid

# Set the terminal to output a high-quality PNG file (optional)
# set terminal pngcairo enhanced font "Arial,12" size 800,600
# set output "MSE_steady_state.png"

# Plot the simulation data from output folder
plot [0:3] \
    'output/multcore_MSE_S2sweep_parallel_t=300_mu=1_S1=0.3.txt' using 1:3 with lines linetype 1 lw 1 title "Theory S1=0.3", \
    'output/parfeval_MSE_S2sweep_t=300_mu=1_S1=0.5.txt' using 1:3 with lines linetype 2 lw 1 title "S1=0.5", \
    'output/multcore_MSE_S2sweep_parallel_t=300_mu=1_S1=1.txt' using 1:3 with lines linetype 3 lw 1 title "S1=1", \
    'output/multcore_MSE_S2sweep_parallel_t=300_mu=1_S1=2.txt' using 1:3 with lines linetype 4 lw 1 title "S1=2", \
    'output/dyhouwaerrorbar(N=200, S1=0.3, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 1 pt 7 ps 0.5 title "Simulation S1=0.3", \
    'output/dyhouwaerrorbar(N=200, S1=0.5, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 2 pt 7 ps 0.5 title "S1=0.5", \
    'output/dyhouwaerrorbar(N=200, S1=1, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 3 pt 7 ps 0.5 title "S1=1", \
    'output/dyhouwaerrorbar(N=200, S1=2, T=300, mu=0.1, ens=100, xi=0_errorbars).txt' using 1:2:3:4 with yerrorbars linetype 4 pt 7 ps 0.5 title "S1=2"


# PNG画像への出力
set terminal pngcairo transparent enhanced font "Arial,12" size 800,500
set output "output/準定常MSE_xi=0.png"
replot
set output
set terminal qt