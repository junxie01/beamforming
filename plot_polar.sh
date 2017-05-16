#!/bin/sh
# created by wangkai
# 2012/10/28
###############################################################################
#
#-------- GMT script----plot the phase/group velocity
# 
#
#######################---parameters set---######################################
# The following variables are used for all the GMT commands.
# Modify the settings here, further down all references in the script
# such as "$region" are replaced by the shell.
gmtset ANOT_FONT_SIZE 14
gmtset HEADER_FONT_SIZE 16
gmtset BASEMAP_TYPE fancy
gmtset PAPER_MEDIA A4
#gmtset X_AXIS_LENGTH 105c
#gmtset Y_AXIS_LENGTH 15C
#gmtset X_ORIGIN 2.5c 
#gmtset Y_ORIGIN 4.5c 
gmtset TICK_LENGTH -0.2c	
dir=/home/junxie/beamforming
cpt=$dir/mymatlabjet.cpt


# set  -P for portrait mode, else landscape
# set  -V for verbose GMT output, else leave blank for quite mode
if [ $# -ne 3 ];then
       echo usage: plot_polar.sh input_file 0 0.5
       exit
fi

#set the input file name
inputfile=$1
#set the output file name
tomo=$1.ps
# set the projection
projection=-JPa5
# set the region (west/east/south/north boundaries)
region=-R0/360/$2/$3
# set tick
tick=-Ba30f10N
# set title
title=:.${1}Hz:
# Set the X and Y offset on the postscript plot (in inches)

################################################################################
# Use gridimage to plot topography
# results.d---three column data need to tomography
# results.m---
# results.grd---grid file by using surface
# mydata.cpt--- color palette table from grid file

cat $inputfile | awk '{print $1, $2, $3}' > results.d

# Get the maximum and minimum value of surface wave group velocity, and this value will be used in grdimage
max=$(cat results.d | awk '{print $3}'|sort -n |sed -n '$p')
min=$(cat results.d | awk '{print $3}'|sort -n -r |sed -n '$p')

echo The max Z value is $max
echo The min Z value is $min
dz=`pwd |awk -v a=$min -v b=$max '{printf"%.1f",(b-a)/6}'`
echo $dz
# Get the maximum and minimum coordinates of seismic station and use this information to clip the suitable tomography region
# first column of input file is latitude and second column is longitude
ymax=$(cat $inputfile | awk '{print $2}'|sort -n |sed -n '$p')
ymin=$(cat $inputfile | awk '{print $2}'|sort -n -r |sed -n '$p')
xmax=$(cat $inputfile | awk '{print $1}'|sort -n |sed -n '$p')
xmin=$(cat $inputfile | awk '{print $1}'|sort -n -r |sed -n '$p')

echo xmin $xmin
echo xmax $xmax
echo ymin $ymin
echo ymax $ymax

#min=`bc -l <<EOF 
#$min + 0
#EOF`

#-------------
blockmedian $inputfile $region -I1/0.005 > results.m
surface results.m -Gresult.grd -I1/0.005 $region -C0.01 -f
#xyz2grd $inputfile -Gresult.grd -I5/0.01 $region 
grd2cpt result.grd -C$cpt -S$min/$max/0.005 -D -Z > mydata.cpt
grdimage result.grd $region $offsets $projection ${tick}${title}  -Cmydata.cpt -K -X2i -Y4i -P >$tomo

if [ 1 -eq 1 ];then
pstext $projection $region -K -O -N -P -G255/0/0 >>$tomo <<EOF
80 0.1 12 0 0 1 0.1
80 0.2 12 0 0 1 0.2
80 0.3 12 0 0 1 0.3
80 0.4 12 0 0 1 0.4
80 0.5 12 0 0 1 0.5
EOF

psxy -JP -R -Sc  -K -N -O -m -W0.950/0/0t5_10:0 >>$tomo <<EOF
0 0 1
0 0 2
0 0 3
0 0 4
EOF
fi
psxy -JP -R0/360/0/$3  -K -N -O -m -W0.950/0/0t5_10:0 >>$tomo <<EOF
>
0 0
0 $3
>
0 0
30 $3
>
0 0
60 $3
>
0 0
90 $3
>
0 0
120 $3
>
0 0
150 $3
>
0 0
180 $3
>
0 0
210 $3
>
0 0
240 $3
>
0 0
270 $3
>
0 0
300 $3
>
0 0
330 $3
EOF



psscale -Cmydata.cpt  -Ba${dz}/:"beampower (db)":  -D12/5/9/0.4 -N -O  >> $tomo

############################################################################################
rm results.d results.m result.grd
#rm *.cpt
#gs $tomo
