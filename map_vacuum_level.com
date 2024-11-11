#! /bin/tcsh -f
#
#	try to find the "vacuum level" of the provided electron
#	density map using "robust statistics"
#
#
set map = "$1"

set tempfile = ${CCP4_SCR}/vacuum$$


# make one ASU
echo "xyzlim ASU" | mapmask mapin $map mapout ${tempfile}test.map > /dev/null
if($status) then
    echo "something is wrong with map: $map"
    exit 9
endif


# examine this map to get the grid and axis convention for this SG
echo "" | mapdump mapin ${tempfile}test.map >! ${tempfile}maphead.txt
set xyzgrid = `awk '/Grid sampling on x, y, z ../{print $8,$9,$10}' ${tempfile}maphead.txt  | tail -1`
set voxels = `awk '/Number of columns, rows, sections/{print $7*$8*$9}' ${tempfile}maphead.txt  | tail -1`
#set mean   = `awk '/Mean density/{print $NF}' ${tempfile}maphead.txt`

set CELL = `awk '/Cell dimensions ../{print $4,$5,$6,$7,$8,$9}' ${tempfile}maphead.txt`
echo $CELL |\
awk 'NF==6{DTR=atan2(1,1)/45; A=cos(DTR*$4); B=cos(DTR*$5); G=cos(DTR*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf "%.3f\n", $1*$2*$3*sqrt(skew)}' |\
cat >! ${tempfile}volume
set CELLvolume = `cat ${tempfile}volume`
rm -f ${tempfile}volume
rm -f ${tempfile}maphead.txt


# convert voxel values to text
set datasize = `echo $voxels | awk '{print 4*$1}'`
set filesize = `ls -l ${tempfile}test.map | awk '{print $5}'`
set head = `echo "$filesize $datasize" | awk '{print $1-$2}'`
echo "reading $map"
od -vf -w4 -j $head ${tempfile}test.map | awk 'NF==2{print $2}' > ! ${tempfile}map.txt
set mean = `awk '{++n;sum+=$1} END{print sum/n}' ${tempfile}map.txt`

echo "cell volume is : $CELLvolume    mean electron density: $mean electrons/A^3"

if ($?histograms) then
    ~jamesh/awk/histogram.awk -v bs=0.01 ${tempfile}map.txt | sort -g >! hist0.txt
endif


# start by looking only at negative side
set cut = `echo $mean | awk '$1<=0{$1=1e-99} {print}'`
set mad = 1e-6

awk -v cut=$cut '$1<cut' ${tempfile}map.txt |\
sort -g >! ${tempfile}sorted.txt

set median = $cut
set lastmedian = 1e99
echo "symmetric median rejection"
while ( $median != $lastmedian ) 
    cat ${tempfile}sorted.txt |\
    awk -v cut=$cut '{print} $1>cut{exit}' |\
    tee ${tempfile}new.txt |\
    awk 'NR==1{min=$1} {++n;v[n]=$1}\
     END{print min+0,v[int(n/2)]+0,n}' >! ${tempfile}median.txt

    set cut = `awk '{print $1+2*($2-$1)}' ${tempfile}median.txt`
    set n = `awk '{print $3}' ${tempfile}median.txt`
    set frac = `echo $n $voxels | awk '{printf "%d",$1/$2*100}'`
    if($frac < 1) goto finish
    set lastmedian = $median
    set median = `awk '{print $2}' ${tempfile}median.txt`
    set min = `head -n 1 ${tempfile}sorted.txt | awk '{print $1+0}'`
    echo "$map median= $median   cut= $cut   frac= ${frac}%   min= $min"
    set median = $min
    mv ${tempfile}new.txt ${tempfile}sorted.txt
end
# now use 3x median absolute deviation as upper rejection cutoff
foreach sig ( 4 3  )
    echo "+$sig x MAD rejection"
set lastmedian = 1e99
while ( $median != $lastmedian ) 
    cat ${tempfile}sorted.txt |\
    awk -v median=$median '{print sqrt(($1-median)^2)}' |\
    sort -g | awk '{++n;v[n]=$1}\
     END{print v[int(n/2)]}' >! ${tempfile}mad.txt
    set mad = `cat ${tempfile}mad.txt`

    cat ${tempfile}sorted.txt |\
    awk -v median=$median -v mad=$mad -v sig=$sig '\
	$1 <= median+sig*mad' |\
    tee ${tempfile}new.txt |\
    awk '{++n;v[n]=$1}\
     END{print v[int(n/2)]+0,n}' >! ${tempfile}median.txt

    set n = `awk '{print $2}' ${tempfile}median.txt`
    set frac = `echo $n $voxels | awk '{printf "%d",$1/$2*100}'`
    if($frac < 1) goto finish
    set min = `head -n 1 ${tempfile}sorted.txt | awk '{print $1+0}'`
    set lastmedian = $median
    set median = `awk '{print $1}' ${tempfile}median.txt`
    echo "$map median= $median   mad= $mad   frac= ${frac}%   min= $min"
    mv ${tempfile}new.txt ${tempfile}sorted.txt
end
end
# now use 2x median absolute deviation as upper and lower rejection cutoff
foreach sig ( 4 3  )
    echo "+/- $sig x MAD rejection"
set lastmedian = 1e99
while ( $median != $lastmedian ) 
    cat ${tempfile}sorted.txt |\
    awk -v median=$median '{print sqrt(($1-median)^2)}' |\
    sort -g | awk '{++n;v[n]=$1}\
     END{print v[int(n/2)]+0}' >! ${tempfile}mad.txt
    set mad = `cat ${tempfile}mad.txt`

    cat ${tempfile}sorted.txt |\
    awk -v median=$median -v mad=$mad -v sig=$sig '\
	sqrt(($1-median)^2) <= sig*mad' |\
    tee ${tempfile}new.txt |\
    awk '{++n;v[n]=$1}\
     END{print v[int(n/2)]+0,n}' >! ${tempfile}median.txt

    set n = `awk '{print $2}' ${tempfile}median.txt`
    set frac = `echo $n $voxels | awk '{printf "%d",$1/$2*100}'`
    if($frac < 1) goto finish
    set min = `head -n 1 ${tempfile}sorted.txt | awk '{print $1+0}'`
    set lastmedian = $median
    set median = `awk '{print $1}' ${tempfile}median.txt`
    echo "$map median= $median   mad= $mad   frac= ${frac}%   min= $min"
    mv ${tempfile}new.txt ${tempfile}sorted.txt
end
end
finish:
set vacuum = `awk '{++n;sum+=$1} END{print sum/n}' ${tempfile}sorted.txt`
echo "$map vacuum level: $vacuum +/- $mad  occupies ${frac}% of map"

if($?histograms) then
    set bs = `echo $mad | awk '$1==0{$1=1e-6} {print $1/10}'`
    awk -v vacuum=$vacuum '{print $1-vacuum}' ${tempfile}sorted.txt |\
    ~jamesh/awk/histogram.awk -v bs=$bs |\
      sort -g >! baseline_hist.txt
    awk -v vacuum=$vacuum '{print $1-vacuum}' ${tempfile}map.txt |\
    ~jamesh/awk/histogram.awk -v bs=$bs |\
      sort -g >! hist.txt
    cat ${tempfile}map.txt |\
    ~jamesh/awk/histogram.awk -v bs=$bs |\
      sort -g >! hist0.txt
endif

set mean = `awk '{++n;sum+=$1} END{print sum/n}' ${tempfile}map.txt`
set F000 = `echo $CELLvolume $mean $vacuum | awk '{print $1*($2-$3)}'`
echo "estimated F000 = $F000"

set offset = `echo $vacuum | awk '{print -$1}'`
mapmask mapin ${tempfile}test.map mapout vacuum_zero.map << EOF > /dev/null
scale factr 1 $offset
EOF
echo "vacuum_zero.map has a vacuum level of zero"


rm -f ${tempfile}test.map
rm -f ${tempfile}map.txt
rm -f ${tempfile}sorted.txt
rm -f ${tempfile}median.txt
rm -f ${tempfile}mad.txt


exit

#####################################################################3
####
#	notes and tests
#
#

~jamesh/Develop/randompdb.com 10 10 10 90 90 90 0.01
set F000 = `egrep "^ATOM" random.pdb | awk '$NF=="H"{e+=1} $NF=="C"{e+=6} $NF=="N"{e+=7} $NF=="O"{e+=8} $NF=="S"{e+=16} END{print e}'`
set CELLvolume = `awk 'BEGIN{print 10*10*10}'`

sfall xyzin random.pdb mapout test.map << EOF
mode atmmap
symm 1
vdwr 3
FORMFACTOR NGAUSS 5
EOF

./map_vacuum_level.com test.map

rm -f test.mtz
phenix.fmodel random.pdb k_sol=0 b_sol=0 high_res=1 \
  algorithm=direct file_name=test.mtz

fft hklin test.mtz mapout test.map << EOF
labin F1=FMODEL PHI=PHIFMODEL
VF000 $CELLvolume $F000
EOF

./map_vacuum_level.com test.map


egrep -v "ANISOU|WAT|HOH" confA.pdb >! dry.pdb
egrep "^ATOM" dry.pdb | awk '$NF=="H"{e+=1} $NF=="C"{e+=6} $NF=="N"{e+=7}\
 $NF=="O"{e+=8} $NF=="S"{e+=16} END{print 4*e}'
set F000 = 3460
set CELLvolume = 11350.989
sfall xyzin dry.pdb mapout test.map << EOF
mode atmmap
symm 19
vdwr 3
grid 128 128 128
EOF
./map_vacuum_level.com test.map

rm -f test.mtz
phenix.fmodel dry.pdb k_sol=0 b_sol=0 high_res=0.5 \
  algorithm=direct file_name=test.mtz

fft hklin test.mtz mapout test.map << EOF
labin F1=FMODEL PHI=PHIFMODEL
VF000 $CELLvolume $F000
EOF


./map_vacuum_level.com test.map




safexp(x) = (x>700 ? 1e300 : (x<-700 ? 0 : exp(x)))

lopsided(x,fwhm,n) = (x<0?safexp(-log(16)*(x/fwhm)**2):1./(1+(2*x/fwhm)**n))


