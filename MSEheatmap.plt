# ====================================================================
# gnuplot script to plot MSE heatmap from simulation data
# ====================================================================

# Set plot options
set xlabel font "Arial,14" "S1"
set ylabel font "Arial,14" "S2"
set cblabel font "Arial,14" "MSE"
set tics font "Arial,14"
set size ratio 1.0
set grid

# Set color palette (jet カラーマップ風)
set palette rgbformulae 33,13,10

# Set view for 2D heatmap (pm3d map)
set view map
set pm3d at b

# 画面にプレビュー表示のみ（保存しない）
set terminal wxt size 800,800
splot 'Ioutput/mse_data_xi=0.00.txt' using 1:2:3 with pm3d notitle


