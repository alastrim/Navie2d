if [ -d H ]; then
	rm -rf H
	rm -rf V1
	rm -rf V2
fi

mkdir H
mkdir V1
mkdir V2
make
./navie $1 $2 $3 $4 $5 $6

if [ ! -d H ]; then
	exit
fi

#1 - T
#2 - t_step_count
#3 - X
#4 - x_step_count
#5 - Y
#6 - y_step_count

gnuplot -p << EOF

f_name = "H"
t_step_count = ($2)
x_point_count = ($4+1)
y_point_count = ($6+1)
T = ($1)
X = ($3)
Y = ($5)

set terminal gif animate delay 5
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
    set hidden3d
    set dgrid3d y_point_count,x_point_count
    set xrange[0:X]
    set yrange[0:Y]
    set zrange[-5:5]
    F = sprintf ("%s/%d", f_name, i)
    splot F with lines title f_name
}

EOF

gnuplot -p << EOF

f_name = "V1"
t_step_count = ($2)
x_point_count = ($4+1)
y_point_count = ($6+1)
T = ($1)
X = ($3)
Y = ($5)

set terminal gif animate delay 5
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
    set hidden3d
    set dgrid3d y_point_count,x_point_count
    set xrange[0:X]
    set yrange[0:Y]
    set zrange[-5:5]
    F = sprintf ("%s/%d", f_name, i)
    splot F with lines title f_name
}

EOF

gnuplot -p << EOF

f_name = "V2"
t_step_count = ($2)
x_point_count = ($4+1)
y_point_count = ($6+1)
T = ($1)
X = ($3)
Y = ($5)

set terminal gif animate delay 5
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
    set hidden3d
    set dgrid3d y_point_count,x_point_count
    set xrange[0:X]
    set yrange[0:Y]
    set zrange[-5:5]
    F = sprintf ("%s/%d", f_name, i)
    splot F with lines title f_name
}

EOF


eog H.gif &
eog V1.gif &
eog V2.gif &
