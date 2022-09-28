set terminal pngcairo enhanced size 600,600 font "Times New Roman, 18"
set output "rdf.png"

#set size square
#set size ratio 0.50

set tmargin at screen 0.85
set lmargin at screen 0.20
set rmargin at screen 0.80
set bmargin at screen 0.25

set title "RDF (POPE-POPE)" offset -0.5,0.0
set key below right font 'Times New Roman,15'
#set key left at graph 0.25,-0.25 font 'Times New Roman,15' 
#set key right top font 'Times New Roman,12'


set xlabel "r (nm)" offset 0.0,0.0 font 'Times New Roman, 20'
set ylabel "g(r)" offset 1.0,0.0 font 'Times New Roman, 20'
set y2label "coordinate number" offset -1.5,0.0 font 'Times New Roman, 20'

set label font "Times New Roman, 18"
#set xrange [0:24] 
#set yrange [0:40.0]
#set y2range [0:40.0]

set xtics nomirror 
set ytics nomirror 
set y2tics nomirror 

#set xtics 10
#set mxtics 10
#unset mxtics
#set ytics 0.2
#mxtics 2
#unset mytics
set tics font "Times New Roman, 14"
set format y "%2.1f"

plot 'GR2_GG_POPE--POPE_CG.dat' u 1:2 axis x1y1 with line lt rgbcolor "red" title "pSPICA",\
	 'GR2_GG_POPE--POPE_AA.dat' u 1:2 axis x1y1 with line lt rgbcolor "orange" title "AA",\
	 'GR2_GG_POPE--POPE_CG.dat' u 1:3 axis x1y2 with line lt rgbcolor "blue" title "pSPICA",\
	 'GR2_GG_POPE--POPE_AA.dat' u 1:3 axis x1y2 with line lt rgbcolor "green" title "AA"

#set terminal pop 
