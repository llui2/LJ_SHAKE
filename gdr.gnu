set terminal pdf 
set output "gdr.pdf"

set title " "
set xlabel "r"
set ylabel "g(r)"

plot 'gdr.dat' using 1:2 w l title ''