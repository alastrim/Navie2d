if [ -d H ]; then
	rm -rf H
	rm -rf V1
	rm -rf V2
	rm -rf realH
	rm -rf realV1
	rm -rf realV2
fi

mkdir H
mkdir V1
mkdir V2
mkdir realH
mkdir realV1
mkdir realV2
make
./navie $1 $2 $3 $4 $5 $6 $7

if [ ! -d H ]; then
	exit
fi

#1 - T
#2 - t_step_count
#3 - X
#4 - x_step_count
#5 - Y
#6 - y_step_count
#7 - anim delay

gnuplot -p << EOF

f_name = "H"
real_name = "realH"
t_step_count = ($2)
t_point_count = (t_step_count+1)
x_step_count = ($4)
x_point_count = (x_step_count+1)
y_step_count = ($6)
y_point_count = (y_step_count+1)
T = ($1)
X = ($3)
Y = ($5)
x_step = (X/x_step_count)
y_step = (Y/y_step_count)
x_start = (0+x_step/2-1)
x_end = (X-x_step/2+1)
y_start = (0+y_step/2-1)
y_end = (Y-y_step/2+1)
grid_row_count = (y_step_count)
grid_column_count = (x_step_count)

set terminal gif animate delay ($7)
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
#    set hidden3d
    set dgrid3d grid_row_count,grid_column_count
    set xrange[x_start:x_end]
    set yrange[y_start:y_end]
    set zrange[-5:5]
    set xlabel "x"
    set ylabel "y"
    set zlabel "z"
    F = sprintf ("%s/%d", f_name, i)
    R = sprintf ("%s/%d", real_name, i)
    splot F with lines title f_name, R with lines title real_name
}

EOF

gnuplot -p << EOF

f_name = "V1"
real_name = "realV1"
t_step_count = ($2)
t_point_count = (t_step_count+1)
x_step_count = ($4)
x_point_count = (x_step_count+1)
y_step_count = ($6)
y_point_count = (y_step_count+1)
T = ($1)
X = ($3)
Y = ($5)
x_step = (X/x_step_count)
y_step = (Y/y_step_count)
x_start = (0-1)
x_end = (X+1)
y_start = (0-1)
y_end = (Y+1)
grid_row_count = (y_step_count+1)
grid_column_count = (x_step_count+1)

set terminal gif animate delay ($7)
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
#    set hidden3d
    set dgrid3d grid_row_count,grid_column_count
    set xrange[x_start:x_end]
    set yrange[y_start:y_end]
    set zrange[-5:5]
    set xlabel "x"
    set ylabel "y"
    set zlabel "z"
    F = sprintf ("%s/%d", f_name, i)
    R = sprintf ("%s/%d", real_name, i)
    splot F with lines title f_name, R with lines title real_name
}

EOF

gnuplot -p << EOF

f_name = "V2"
real_name = "realV2"
t_step_count = ($2)
t_point_count = (t_step_count+1)
x_step_count = ($4)
x_point_count = (x_step_count+1)
y_step_count = ($6)
y_point_count = (y_step_count+1)
T = ($1)
X = ($3)
Y = ($5)
x_step = (X/x_step_count)
y_step = (Y/y_step_count)
x_start = (0-1)
x_end = (X+1)
y_start = (0-1)
y_end = (Y+1)
grid_row_count = (y_step_count+1)
grid_column_count = (x_step_count+1)

set terminal gif animate delay ($7)
set output sprintf ("%s.gif", f_name)

do for [i=0:t_step_count] {
#    set hidden3d
    set dgrid3d grid_row_count,grid_column_count
    set xrange[x_start:x_end]
    set yrange[y_start:y_end]
    set zrange[-5:5]
    set xlabel "x"
    set ylabel "y"
    set zlabel "z"
    F = sprintf ("%s/%d", f_name, i)
    R = sprintf ("%s/%d", real_name, i)
    splot F with lines title f_name, R with lines title real_name
}

EOF

eog V1.gif &
eog V2.gif &
eog H.gif &
