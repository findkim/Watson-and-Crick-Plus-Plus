#
# Creates stacked MinMax plots (5) with pValue plot below
# auto-layout capability to label organism and axis tics
#
#

set terminal postscript portrait
set output "Ypes.fasta.ps"

#set terminal latex
#set output "stacked_min_max_plots.tex"

# Top & bottom margin
set tmargin 0
set bmargin 1

# Left & right margin to have same size plots
set lmargin 5
set rmargin 5

# Data columns are separated with ,
set datafile separator ","


# Turn off xtics for all plots but the bottom one
unset xtics
unset ytics

# Displays origin of x-axis
set xzeroaxis lt 1 linecolor rgb "#000000"

set yrange[-100:100]
#set format y "rare	common"
set ytics (" Common" 50,  "Rare" -50)
set ytics axis in scale 0,0 nomirror rotate by 90 font "Helvetica,10"
#set ytics font "Helvetica,5"

set key autotitle columnhead
set key inside left bottom vertical nobox
set key samplen -1

set multiplot layout 9, 1 title "MinMax plots for Ypes.fasta"
set style histogram rowstacked
set style data histogram

# Histogram bars filled with solid color and black border
set style fill solid noborder


#
#	MinMax Plots
#
# plot 'datafile' using #column
plot 'Ypes.fasta.csv' using ($1 > 0 ? $1 : 0) linecolor rgb "#0000FF", "" using ($1 < 0 ? $1 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($2 > 0 ? $2 : 0) linecolor rgb "#0000FF", "" using ($2 < 0 ? $2 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($3 > 0 ? $3 : 0) linecolor rgb "#0000FF", "" using ($3 < 0 ? $3 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($4 > 0 ? $4 : 0) linecolor rgb "#0000FF", "" using ($4 < 0 ? $4 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($5 > 0 ? $5 : 0) linecolor rgb "#0000FF", "" using ($5 < 0 ? $5 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($6 > 0 ? $6 : 0) linecolor rgb "#0000FF", "" using ($6 < 0 ? $6 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($7 > 0 ? $7 : 0) linecolor rgb "#0000FF", "" using ($7 < 0 ? $7 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($8 > 0 ? $8 : 0) linecolor rgb "#0000FF", "" using ($8 < 0 ? $8 : -0.05) linecolor rgb "#FF00FF"

plot 'Ypes.fasta.csv' using ($9 > 0 ? $9 : 0) linecolor rgb "#0000FF", "" using ($9 < 0 ? $9 : -0.05) linecolor rgb "#FF00FF"



#
#	pValue plot
#
set style histogram

# Resets ytics and xtics preset from earlier
unset ytics
unset xtics

# Autoscales x and y axis according to data
set autoscale y
set autoscale x

# Nomirror (only on left y axis & bottom x axis)
# in/out draws tics inwards/outwards of plot
# Startrange, increment, {optional end range}
set ytics nomirror out 5 font "Helvetica,10"
set xtics nomirror out 0,20 font "Helvetica,10"

#plot 'sample_min_max2.txt' using 15 linecolor rgb "#FF0000"

#
unset multiplot
#pause -1
#
#
#
# maybe two pages
# unset multiplot
# unset title
# set multiplot
