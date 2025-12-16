# ====================================================================
# gnuplot script to plot MSE(t) from simulation data
# ====================================================================

# プレビュー（ウィンドウ）表示
set title "μ=1" font "Arial,20"
set yrange [0:1]
set xlabel font "Arial,20" "t = n/N"
set ylabel font "Arial,20" "MSE"
set tics font "Arial,12"
set key right top
set size ratio 0.75

# ----- ここで一度plotする（qt/wxtなどウィンドウ端末向け） -----
# （端末指定を省略すればデフォルトの画面描画端末になります）
plot [0:100] \
    'Ioutput/dyhouwa1int,tEnd=100,mu=1,S1=1,S2=0.3.txt' using 1:4 with lines linetype 1 lw 1 title "Theory S2=0.3", \
    'Ioutput/dyhouwa1int,tEnd=100,mu=1,S1=1,S2=0.5.txt' using 1:4 with lines linetype 2 lw 1 title "S2=0.5", \
    'Ioutput/dyhouwa1int,tEnd=100,mu=1,S1=1,S2=1.txt' using 1:4 with lines linetype 3 lw 1 title "S2=1", \
    'Ioutput/dyhouwa1int,tEnd=100,mu=1,S1=1,S2=2.txt' using 1:4 with lines linetype 4 lw 1 title "S2=2", \
    'Ioutput/dyhouwaMSE Simulation(N=200,T=5000,S1=1,S2=0.3,mu=1,ens=1000,xi=0).txt' using 1:2 with points linetype 1 pt 5 ps 0.5 title "Simulation S2=0.3", \
    'Ioutput/dyhouwaMSE Simulation(N=200,T=5000,S1=1,S2=0.5,mu=1,ens=1000,xi=0).txt' using 1:2 with points linetype 2 pt 7 ps 0.5 title "S2=0.5", \
    'Ioutput/dyhouwaMSE Simulation(N=200,T=5000,S1=1,S2=1,mu=1,ens=1000,xi=0).txt' using 1:2 with points linetype 3 pt 9 ps 0.5 title "S2=1", \
    'Ioutput/dyhouwaMSE Simulation(N=200,T=5000,S1=1,S2=2,mu=1,ens=1000,xi=0).txt' using 1:2 with points linetype 4 pt 11 ps 0.5 title "S2=2"

# PNG画像への出力
set terminal pngcairo transparent enhanced font "Arial,12" size 800,500
set output "output/MSE_S1=1_mu=1.png"
replot
set output
set terminal qt