gnuplot -p << EOF

set dgrid3d 30,30
set hidden3d
splot "H" u 1:2:3 with lines

EOF
