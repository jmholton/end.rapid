#! /bin/tcsh -f
#
#	change FP by RMS SIGFP				-James Holton  3-2-12
#
#


set mtzfile = refme.mtz
set F = ""
set SIGF = ""
set seed = `date +%N | awk '{print $1/1000}'`


set tempfile = ${CCP4_SCR}/kick_data$$


foreach arg ( $* )

    if( "$arg" =~ *.mtz ) set mtzfile  = "$arg"
    if(( "$arg" =~ *[0-9] )&&( "$arg" =~ [1-9]* )) set steps = "$arg"

    if( "$arg" =~ seed=* ) then
        set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
        if("$test" != "") set seed = $test
    endif
    if( "$arg" =~ F=* ) then
        set test = `echo $arg | awk -F "[=]" '$2+0>0{print $2+0}'`
        if("$test" != "") set user_F = $test
    endif
end

if(! -e "$mtzfile") then
    set BAD = "$mtzfile does not exist"
    goto exit
endif

# examine MTZ file
echo "go" | mtzdump hklin $mtzfile |\
awk '/OVERALL FILE STATISTICS/,/No. of reflections used/' |\
awk 'NF>5 && $(NF-1) ~ /^[FJDGKQLMPWABYIRUV]$/' |\
cat >! ${tempfile}mtzdmp

# use completeness, or F/sigF to pick default F
cat ${tempfile}mtzdmp |\
awk '$(NF-1) == "F"{F=$NF; meanF=$8; reso=$(NF-2); comp=substr($0,32)+0; \
  getline; if($(NF-1) != "Q") next; \
  S=$NF; if($8) meanF /= $8; print F, S, reso, comp, meanF;}' |\
sort -k3n,4 -k4nr,5 -k5nr >! ${tempfile}F

# and extract all dataset types/labels
cat ${tempfile}mtzdmp |\
awk 'NF>2{print $(NF-1), $NF, " "}' |\
cat >! ${tempfile}cards

#clean up
rm -f ${tempfile}mtzdmp

if("$F" == "" || "$SIGF" == "") then
    # pick F with best resolution, or F/sigma
    if($?user_F) then
	set F = `awk -v F=$user_F '$1==F{print}' ${tempfile}F`
    endif
    if($#F < 2) set F = `head -1 ${tempfile}F`
    if($#F > 2) then
	set SIGF = $F[2]
	set F    = $F[1]
	echo "selected F=$F SIGF=$SIGF "
    endif
    rm -f ${tempfile}F
endif
# capture remaining columns for imort later...
set otherstuff = `awk -v F=$F -v SIGF=$SIGF '$2!=F && $2!=SIGF{++n;print "E"n"="$2}' ${tempfile}cards`
rm -f ${tempfile}cards
#echo "$otherstuff"

# extract only F/SIGF so sftools does not get pissed off
cad hklin1 $mtzfile hklout F_SIGF.mtz << EOF > /dev/null
labin file 1 E1=$F E2=$SIGF
EOF

echo "using seed = $seed"
sftools << EOF >&! ${tempfile}sftools.log
read F_SIGF.mtz
calc seed $seed
calc col noise = col $SIGF RAN_G *
calc F col Fnew  = col $F col noise +
select col Fnew < 0
calc col Fnew = 0
select all
calc F col $F = col Fnew
delete  col Fnew noise
write new.mtz 
y
exit
y
EOF
if($status) then
    set BAD = "sftools failed"
    goto exit
endif
if("$otherstuff" != "") then
    cad hklin1 new.mtz hklin2 $mtzfile hklout kicked.mtz << EOF > /dev/null
labin file 1 all
labin file 2 $otherstuff
EOF
else
    cad hklin1 new.mtz hklout kicked.mtz << EOF > /dev/null
labin file 1 all
EOF
endif
if($status) then
    set BAD = "output mtz corrupted."
    goto exit
endif
rm -f F_SIGF.mtz
rm -f new.mtz

echo "kicked.mtz contains $F from $mtzfile modified by rms $SIGF"

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

exit


######################################################################3
#
#	notes and stuff
#

set seed = 1

foreach id ( `cat pdb_ids.txt` )

set ID = `echo $id | awk '{print toupper($1)}'`

foreach n ( `seq 1 20` )
    @ seed = ( $seed + 1 )
    echo "$seed" | tee -a ${ID}/seeds.txt
end

end



setenv PHENIX_OVERWRITE_ALL true

set pwd = `pwd`
set seed = 1

foreach id ( `cat pdb_ids.txt` )

set ID = `echo $id | awk '{print toupper($1)}'`
cd $ID 

if(-e busy || -e done) continue
touch busy

foreach seed ( `cat seeds.txt` )

../kick_data.com *_data.mtz seed=$seed

set ligands = `ls -1rt | awk '/cif$/'`
phenix.refine kicked.mtz *updated.pdb $ligands

mv *_map_coeffs.mtz seed_${seed}_map_coeffs.mtz

end

mv busy done
end



