#! /bin/tcsh -f
#
#
#	compute the average and variance of a series of maps
#
#	uses auxillary program float_mult				-James Holton 1-20-12
#
#
set tempfile = ${CCP4_SCR}/mapvar$$

set ref = ""
set maps = ""
foreach arg ( $* )
    
    if("$arg" =~ ref=*) then
	set ref = `echo $arg | awk '{print substr($0,5)}'`
	continue
    endif
    if("$arg" =~ *.map) then
	set maps = ( $maps $arg )
	continue
    endif
    if("$arg" =~ -norms* || "$arg" =~ norms*) then
	set NO_RMSD
	continue
    endif
end

if($#maps == 0) then
    echo "usage: $0 *.map"
    exit 9
endif
if(-e "$ref") goto got_ref

echo "scale factor 0 0" |\
  mapmask mapin $maps[1] mapout ${tempfile}sum.map > /dev/null

foreach map ( $maps )
    echo "$map"

    echo "maps add" |\
    mapmask mapin1 $map mapin2 ${tempfile}sum.map mapout ${tempfile}temp.map > /dev/null
    mv ${tempfile}temp.map ${tempfile}sum.map
end
set scale = `echo $#maps | awk '{print 1/$1}'`
echo "scale factor $scale 0" |\
  mapmask mapin ${tempfile}sum.map mapout avg.map > /dev/null
set ref = avg.map
echo "avg.map is the average of all maps"
if($?NO_RMSD) goto exit

got_ref:
echo "scale factor -1 0" |\
  mapmask mapin $ref mapout ${tempfile}negavg.map > /dev/null
echo "scale factor 0 0" |\
  mapmask mapin $ref mapout ${tempfile}sum.map > /dev/null

foreach map ( $maps )
    echo "$map - $ref"

    echo "maps add" |\
    mapmask mapin1 $map mapin2 ${tempfile}negavg.map mapout ${tempfile}temp.map > /dev/null
    echo "maps mult" |\
    mapmask mapin1 ${tempfile}temp.map mapin2 ${tempfile}temp.map mapout ${tempfile}sqrdiff.map > /dev/null
    echo "maps add" |\
    mapmask mapin1 ${tempfile}sum.map mapin2 ${tempfile}sqrdiff.map mapout ${tempfile}temp.map > /dev/null
    mv ${tempfile}temp.map ${tempfile}sum.map
end
if($#maps != 1) then
    set scale = `echo $#maps | awk '{print 1/($1-1)}'`
else
    echo "WARNING: only one map! output will simply be absolute difference."
    set scale = 1
endif
echo "scale factor $scale 0" |\
  mapmask mapin ${tempfile}sum.map mapout variance.map > /dev/null


rm -f ${tempfile}temp.map ${tempfile}sum.map ${tempfile}sqrdiff.map ${tempfile}negavg.map

echo "variance.map is the variance from $ref"

echo go | mapdump mapin variance.map | tee ${tempfile}mapdump.log |\
awk '/Grid sampling on x, y, z/{gx=$8;gy=$9;gz=$10}\
     /Maximum density /{max=$NF}\
     /Cell dimensions /{xs=$4/gx;ys=$5/gy;zs=$6/gz}\
     /Number of columns, rows, sections/{nc=$7;nr=$8;ns=$9}\
 END{print xs,ys,zs,nc,nr,ns,max}' >! ${tempfile}mapstuff.txt

set size = `awk '{print 4*$4*$5*$6}' ${tempfile}mapstuff.txt`
set head = `ls -l variance.map | awk -v size=$size '{print $5-size}'`
set skip = `echo $head | awk '{print $1+1}'`


# now somehow take the square root?
sqrt:
if(-e float_mult) set path = ( . $path )
if(-e `dirname $0`/float_mult ) set path = ( `dirname $0` $path )
rm -f sigma.map ${tempfile}output.bin
head -c $head variance.map >! ${tempfile}temp.map
tail -c +$skip variance.map >! ${tempfile}variance.bin
float_mult ${tempfile}variance.bin ${tempfile}variance.bin ${tempfile}output.bin -power1 0.5 -power2 0 > /dev/null
if(-e ${tempfile}output.bin) then
    cat ${tempfile}output.bin >> ${tempfile}temp.map
    echo "scale factor 1 0" | mapmask mapin ${tempfile}temp.map mapout sigma.map > /dev/null
    rm -f ${tempfile}output.bin
    echo "sigma.map is the rms deviation from $ref"
else
    if(! -e float_mult.c) goto compile_float_mult
endif
rm -f ${tempfile}temp.map
rm -f ${tempfile}variance.bin
rm -f ${tempfile}output.bin


echo "sigma.map :"
echo "go" | mapdump mapin sigma.map | egrep density


exit:

exit

#############################################################################
#############################################################################


compile_float_mult:
echo "attempting to generate float_mult utility ..."
cat << EOF >! float_mult.c

/* multiply two binary "float" files together                                           -James Holton           1-31-10

example:

gcc -O -O -o float_mult float_mult.c -lm
./float_mult file1.bin file2.bin output.bin 

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


char *infile1name = "";
char *infile2name = "";
FILE *infile1 = NULL;
FILE *infile2 = NULL;
char *outfilename = "output.bin\\0";
FILE *outfile = NULL;

int main(int argc, char** argv)
{
     
    int n,i,j,k,pixels;
    float *outimage;
    float *inimage1;
    float *inimage2;
    float power1=1.0,power2=1.0;
    float sum,sumd,sumsq,sumdsq,avg,rms,rmsd,min,max;
        
    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(strlen(argv[i]) > 4)
        {
            if(strstr(argv[i]+strlen(argv[i])-4,".bin"))
            {
                printf("filename: %s\\n",argv[i]);
                if(infile1 == NULL){
                    infile1 = fopen(argv[i],"r");
                    if(infile1 != NULL) infile1name = argv[i];
                }
                else
                {
                    if(infile2 == NULL){
                        infile2 = fopen(argv[i],"r");
                        if(infile2 != NULL) infile2name = argv[i];
                    }
                    else
                    {
                        outfilename = argv[i];
                    }
                }
            }
        }

        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-power1") && (argc >= (i+1)))
            {
                power1 = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-power2") && (argc >= (i+1)))
            {
                power2 = atof(argv[i+1]);
            }
        }
    }

    if(infile1 == NULL || infile2 == NULL){
        printf("usage: float_mult file1.bin file2.bin [outfile.bin] -power1 1.0 -power2 1.0\\n");
        printf("options:\\n");\\
//        printf("\\t-atom\\tnumber of atoms in the following file\\n");
//      printf("\\t-file filename.txt\\ttext file containing point scatterer coordinates in Angstrom relative to the origin.  The x axis is the x-ray beam and Y and Z are parallel to the detector Y and X coordinates, respectively\\n");
exit(9);
    }


    /* load first float-image */
    fseek(infile1,0,SEEK_END);
    n = ftell(infile1);
    rewind(infile1);
    inimage1 = calloc(n,1);
    inimage2 = calloc(n,1);
    fread(inimage1,n,1,infile1);
    fclose(infile1);
    fread(inimage2,n,1,infile2);
    fclose(infile2);

    pixels = n/sizeof(float);
    outfile = fopen(outfilename,"w");
    if(outfile == NULL)
    {
        printf("ERROR: unable to open %s\\n", outfilename);
        exit(9);
    }
    

    outimage = calloc(pixels,sizeof(float));
    sum = sumsq = sumd = sumdsq = 0.0;
    min = 1e99;max=-1e99;
    for(i=0;i<pixels;++i)
    {
        if(inimage1[i]<0.0 && power1 != ((int) power1) ) inimage1[i] = 0.0;
        if(inimage2[i]<0.0 && power2 != ((int) power2) ) inimage2[i] = 0.0;
        outimage[i] = powf(inimage1[i],power1) * powf(inimage2[i],power2);
        if(outimage[i]>max) max=outimage[i];
        if(outimage[i]<min) min=outimage[i];
        sum += outimage[i];
        sumsq += outimage[i]*outimage[i];
    }
    avg = sum/pixels;
    rms = sqrt(sumsq/pixels);
    for(i=0;i<pixels;++i)
    {
        sumd   += outimage[i] - avg;
        sumdsq += (outimage[i] - avg) * (outimage[i] - avg);
    }
    rmsd = sqrt(sumdsq/pixels);
    printf("max = %g min = %g\\n",max,min);
    printf("mean = %g rms = %g rmsd = %g\\n",avg,rms,rmsd);


    printf("writing %s as %d %d-byte floats\\n",outfilename,pixels,sizeof(float));
    outfile = fopen(outfilename,"w");
    fwrite(outimage,pixels,sizeof(float),outfile);
    fclose(outfile);


    return;
}

EOF
gcc -o float_mult float_mult.c -lm -static
set path = ( . $path )
goto sqrt



