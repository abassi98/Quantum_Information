#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 8    last modified 2019-12-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2019
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
set style line 1  linecolor rgb "#00ffff"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 2  linecolor rgb "#00bfff"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 3  linecolor rgb "#0000ff"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 4  linecolor rgb "#000080"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 5  linecolor rgb "#0000ff"  linewidth 1.000 dashtype solid pointtype 6 pointsize 1
set style line 6  linecolor rgb "#00008b"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.1
set style line 7  linecolor rgb "#ff0000"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 8  linecolor rgb "#ff0000"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.2
set style line 9  linecolor rgb "#6495ed"  linewidth 1.000 dashtype solid pointtype 6 pointsize 0.5
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels 5
set cntrparam levels auto
set cntrparam firstlinetype 0 unsorted
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" textcolor lt -1 norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" textcolor lt -1 norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "it_IT.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
s=0.5
omega=0.05
maxeigen=10
f(x)=omega*(x+s)
fit [0:maxeigen]f(x) 'data/eigenvalues.txt' u 1:2 via omega,s
omega=0.05
#define the hermite polynomials
p0(x)=(omega**(0.25))*exp(-omega*(x**2)/2)*0.7511255445
p1(x)=2*(x*(omega**(0.5)))*(omega**(0.25))*exp(-omega*(x**2)/2)*0.5311259660383879
p2(x)=(4*(x*(omega**(0.5)))**2-2)*(omega**(0.25))*exp(-omega*(x**2)/2)*0.26556298301919395
p3(x)=(8*(x*(omega**(0.5)))**3-12*(x*(omega**(0.5))))*(omega**(0.254))*exp(-omega*(x**2)/2)*0.1084156338280698
p4(x)=(16*(x*(omega**(0.5)))**4-48*(x*(omega**(0.5)))**2+12)*(omega**(0.25))*exp(-omega*(x**2)/2)*0.03833071493323291
p5(x)=(32*(x*(omega**(0.5)))**5-160*(x*(omega**(0.5)))**3+120*(x*(omega**(0.5))))*(omega**(0.25))*exp(-omega*(x**2)/2)*0.012121236353164491 #max 6 states


set terminal postscript enhanced
#cut from now on
set key left 

xmax=2

ymax=1.4
set term png
set output 'data/eigenvalues.png'
set xlabel "n"
set ylabel "E"
p [0:2*maxeigen]'data/eigenvalues.txt' u 1:2 w p ls 5 ti 'Data', f(x) ti 'Fit' ls 6
set xlabel "x"
set ylabel ""
set output 'data/eigenvector0.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$2) w p ls 9 ti 'State n=0', p0(x) ls 8 ti 'Theor'

set output 'data/eigenvector1.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$3) w p ls 9 ti 'State n=1', p1(x) ls 8 ti 'Theor'
 
set output 'data/eigenvector2.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:($4) w p ls 9 ti 'State n=2', p2(x) ls 8 ti 'Theor'

set output 'data/eigenvector3.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$5) w p ls 9 ti 'State n=3', p3(x) ls 8 ti 'Theor'
set output 'data/eigenvector4.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:($6) w p ls 9 ti 'State n=4', p4(x) ls 8 ti 'Theor'
set output 'data/eigenvector5.png' 
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$7) w p ls 9 ti 'State n=5', p5(x) ls 8 ti 'Theor' 

set output 'data/cumulative.png'
p [-xmax/2:xmax/2][-ymax:ymax]'data/eigenvectors.txt' u 1:(-$2) w p ls 1 ti '{/Symbol y}_0', 'data/eigenvectors.txt' u 1:(-$3) w p ls 2 ti '{/Symbol y}_1','data/eigenvectors.txt' u 1:(-$4) w p ls 3 ti '{/Symbol y}_2','data/eigenvectors.txt' u 1:($5) w p ls 4 ti '{/Symbol y}_3',((omega*x)**2)/2 ls 7 ti 'Potential'


set ylabel "Residuals"
set output 'data/res0.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$2-p0($1)) w p ls 3 notitle

set output 'data/res1.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$3-p1($1)) w p ls 3 notitle

set output 'data/res2.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:($4-p2($1)) w p ls 3 notitle

set output 'data/res3.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$5-p3($1)) w p ls 3 notitle

set output 'data/res4.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:($6-p4($1)) w p ls 3 notitle

set output 'data/res5.png'
p [-xmax/2:xmax/2]'data/eigenvectors.txt' u 1:(-$7-p5($1)) w p ls 3 notitle




