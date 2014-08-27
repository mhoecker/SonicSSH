set xdata time
set timefmt "%s"
set format x "%m/%d %H:%M"
set xtics rotate by -30
set xtics 
set grid
set rmargin 10
set xtics 900
set yrange [0:1024]
set ylabel "Distance Measured (cm)"
set term pdf enhanced color size 9in,6in font "Helvetica,12" linewidth 2
#
set output "10kts.pdf"
set title "Test at 10 kts"
set xlabel "UTC"
set label "Distance to sensor (cm)"
plot "LOGGER11.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5  w lines lt -1 notitle
#
set output "CTD5.pdf"
set title "CTD Station 5"
plot "LOGGER12.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5  w lines lt -1 notitle
#
set output "VMP4.pdf"
set title "VMP Station 4"
plot "LOGGER13.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5  w lines lt -1 notitle
#
set output "VMP6.pdf"
set title "VMP Station 6"
plot "LOGGER14.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5  w lines lt -1 notitle
#
set output "CTD51.pdf"
set title "CTD Station 51"
plot "LOGGER15.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5 w lines lt -1 notitle
#
set output "VMP12.pdf"
set title "VMP Station 12"
plot "LOGGER15.CSV" u ($2+(floor($1)%1000)/1000+7*3600):5 w lines lt -1 notitle
