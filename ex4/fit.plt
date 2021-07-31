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
set key fixed left top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
set label 1 "" at 0.00000, 0.00000, 0.00000 left norotate back nopoint
unset arrow
set style increment default
unset style line
set style line 1  linecolor rgb "#ff007f"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
set style line 2  linecolor rgb "#33ffff"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
set style line 3  linecolor rgb "#3333ff"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
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
set xlabel "Size" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 0.00000 : 4000.00 ] noreverse writeback
set x2range [ 1.00000 : 4000.00 ] noreverse writeback
set ylabel "Time(s)" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 0.00000 : 450.000 ] noreverse writeback
set y2range [ 1.00001e-06 : 405.851 ] noreverse writeback
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
set fit logfile 'logfile.log' quiet errorvariables nocovariancevariables errorscaling prescale nowrap v5
f1(x)=a1+b1*x**c1
f2(x)=a2+b2*x**c2
f3(x)=a3+b3*x**c3
GNUTERM = "qt"
a1 = 9.99999930468876e-06
b1 = 1.46998410465353e-09
c1 = 3.17625306797556
a2 = 1.00062639892295e-05
b2 = 6.18134606081662e-09
c2 = 2.83605559273062
a3 = 9.99999999996983e-06
b3 = 2.33678521816486e-09
c3 = 2.54422520216536
x = 0.0
FIT_CONVERGED = 1
FIT_NDF = 66
FIT_STDFIT = 4.28478323707378
FIT_WSSR = 1211.71824765476
FIT_P = 0.0
FIT_NITER = 1
a1_err = 0.588301160457171
b1_err = 5.45965610873233e-10
c1_err = 0.044742831227472
a2_err = 0.0551435166924764
b2_err = 5.0716550456311e-10
c2_err = 0.00988271442547673
a3_err = 0.0148074062616047
b3_err = 1.28324903824407e-09
c3_err = 0.0663015406937224
fit f1(x) "std_mult.txt" u 1:2 via a1,b1,c1
fit f2(x) "clmn_mult.txt" u 1:2 via a2,b2,c2
fit f3(x) "built_in.txt" u 1:2 via a3,b3,c3
g1(x)=h1+k1*x
g2(x)=h2+k2*x
g3(x)=h3+k3*x
GNUTERM = "qt"
h1 = -17.5272321856737
k1 = 2.74639633630251
h2 = -16.0947260940355
k2 = 2.45702050809165
h3 = -20.5473824111374
k3 = 2.59161979957259
FIT_CONVERGED = 1
FIT_NDF = 54
FIT_STDFIT = 0.235370946825112
FIT_WSSR = 2.99157206090489
FIT_P = 1.0
FIT_NITER = 4
h1_err = 0.267398597891585
k1_err = 0.0433298532847444
h2_err = 0.104113483044413
k2_err = 0.0168707763647541
h3_err = 0.149977540871681
k3_err = 0.0243026885451825
fit   g1(x) "std_mult.txt" u (log($1)):(log($2)) via h1,k1
fit  g2(x) "clmn_mult.txt" u (log($1)):(log($2)) via h2,k2
fit  g3(x) "built_in.txt" u (log($1)):(log($2)) via h3,k3

set term png
set autoscale
set output "fit.png"
p "std_mult.txt" u 1:2 w p ls 1 ti "Standard-mult", f1(x) ti "" ls 1,"clmn_mult.txt" u 1:2 w p ls 2 ti "Col-par-col",  f2(x) ti "" ls 2, "built_in.txt" u 1:2 w p ls 3 ti "Built-in", f3(x) ti "" ls 3

set ylabel "Residual(s)"
set output "res_fit.png"
p "std_mult.txt" u 1:($2-f1($1)) w p ls 1 ti "Standard-mult", "clmn_mult.txt" u 1:($2-f2($1)) w p ls 2 ti "Col-par-col", "built_in.txt" u 1:($2-f3($1))w p ls 3 ti "Built-in"
set output "log.png"
set ylabel "log(Time)"
set xlabel "log(Size)"
## Last datafile plotted: "built_in.txt"
p "std_mult.txt" u (log($1)):(log($2)) w p ls 1 ti "Standard-mult",g1(x) ti "" ls 1, "clmn_mult.txt" u (log($1)):(log($2)) w p ls 2 ti "Col-par-col",g2(x) ti "" ls 2, "built_in.txt" u (log($1)):(log($2)) w p ls  3 ti  "Built-in", g3(x) ti "" ls 3
set ylabel "Residual(log(Time))"
set autoscale
set output "res_log.png"
set key right
p "std_mult.txt" u (log($1)):(log($2)-g1(log($1))) w p ls 1 ti "Standard-mult", "clmn_mult.txt" u (log($1)):(log($2)-g2(log($1))) w p ls 2 ti "Col-par-col", "built_in.txt" u (log($1)):(log($2)-g3(log($1))) w p ls 3 ti "Built-in"
#    EOF
j1(x)=exp(1)**h1*x**k1
j2(x)=exp(1)**h2*x**k2
j3(x)=exp(1)**h3*x**k3
set autoscale
set key left
set xlabel "Size"
set ylabel "Time(s)"
set output "Comparison.png"
p "std_mult.txt" u 1:2 w p ls 1 ti "Standard-mult", j1(x) ls 1 ti "", "clmn_mult.txt" u 1:2 w p ls 2 ti "Col-par-col", j2(x) ls 2 ti "", "built_in.txt" u 1:2 w p ls 3 ti "Built-in", j3(x) ls 3 ti ""
set ylabel "Residual(s)"
set output "Comparison_res.png"
p "std_mult.txt" u 1:($2-j1($1)) w p ls 1 ti "Standard-mult", "clmn_mult.txt" u 1:($2-j2($1)) w p ls 2 ti "Col-par-col", "built_in.txt" u 1:($2-j3($1)) w p ls 3 ti "Built-in"



