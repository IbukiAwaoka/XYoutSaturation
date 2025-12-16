# ====================================================================
# gnuplot script to plot MSE heatmap from simulation data
# ====================================================================

# Set plot options
#set title "Sim Steady State MSE Heatmap"
set xlabel font "Arial,28" "S_1" offset 0,-1
set ylabel font "Arial,28" "S_2" offset -2.8,0
set cbrange [0:1]
set cbtics 0.1
set cblabel font "Arial,28" "MSE" offset 2.5,0
set format x "%.1f"
set format y "%.1f"
set format cb "%.1f"
set tics font "Arial,16"
set size ratio 1.0
set grid

# Set color palette (jet カラーマップ風)
set palette rgbformulae 33,13,10

# Set view for 2D heatmap (pm3d map)
set view map
set pm3d at b

# 画面にプレビュー表示のみ（保存しない）
set terminal wxt size 800,800
splot 'Ioutput/NEWmse_data_xi=0.00.txt' using 1:2:3 with pm3d notitle


