set term tikz size 6.0cm, 5.0cm fontscale 0.5 dl 0.5 fulldoc
#set terminal pngcairo
set output "E2-tex.tex"
#set output "E2.png"
set tics nomirror
set tics scale 0.5
set decimalsign ","
#set xrange[0:2.0]

set grid
set title "Energia Média ao Quadrado vs Temperatura"
set xlabel "Temperatura $T$"
set ylabel "Energia média ao Quadrado $\\langle E \\rangle^2$"
set key at graph 0.90, 0.30

plot "../../data/output_5.dat" u 1:3 w lp pt 7 lw 2 ps 0.3 lc 6 t'L=05',"../../data/output_10.dat" u 1:3  w lp pt 7 lw 2 ps 0.3 lc 2 t'L=10',"../../data/output_20.dat" u 1:3 w  lp pt 7 lw 2 ps 0.3 lc 7 t'L=20', "../../data/output_40.dat" u 1:3 w  lp pt 7 lw 2 ps 0.3 lc 4 t'L=40'


unset output
