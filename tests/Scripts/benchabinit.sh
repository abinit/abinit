#
# Provides a summary of the executing time, to be placed eventually on the
# Web page.
#
# Script to be called when executing all sequential tests (for benchmarks):
#
#    make tests_allseq
#
# Otherwise, can be executed from a tmp-* directory inside tests/cpu.
#
cd A3
echo "fourwf(pot)"
f20=`grep 'fourwf(pot)' report |head -1|awk '{print $NF}'`
f30=`grep 'fourwf(pot)' report |head -3|tail -1|awk '{print $NF}'`
f36=`grep 'fourwf(pot)' report |head -5|tail -1|awk '{print $NF}'`
echo "20 : $f20"
echo "30 : $f30"
echo  "36 : $f36"
cd ../B3
f48=`grep 'fourwf(pot)' report |head -1|awk '{print $NF}'`
f64=`grep 'fourwf(pot)' report |head -3|tail -1|awk '{print $NF}'`
f80=`grep 'fourwf(pot)' report |head -5|tail -1|awk '{print $NF}'`
f96=`grep 'fourwf(pot)' report |tail -1|awk '{print $NF}'`
echo "48 : $f48"
echo "64 : $f64"
echo "80 : $f80"
echo "96 : $f96"
#av=`echo $f20 $f30 $f36 $f48 $f64 $f80 $f96|awk '{print ($1+$2+$3+$4+$5+$6+$7)/7}'`
#echo 'fourwf(pot)average= '$av

cd ../A3
#echo "nonlop(apply)"
nl20=`grep 'nonlop(apply)' report |head -1|awk '{print $NF}'`
nl30=`grep 'nonlop(apply)' report |head -3|tail -1|awk '{print $NF}'`
nl36=`grep 'nonlop(apply)' report |head -5|tail -1|awk '{print $NF}'`
#echo "20 : $nl20"
#echo "30 : $nl30"
#echo "36 : $nl36"
cd ../B3
nl48=`grep 'nonlop(apply)' report |head -1|awk '{print $NF}'`
nl64=`grep 'nonlop(apply)' report |head -3|tail -1|awk '{print $NF}'`
nl80=`grep 'nonlop(apply)' report |head -5|tail -1|awk '{print $NF}'`
nl96=`grep 'nonlop(apply)' report |tail -1|awk '{print $NF}'`
#echo "48 : $nl48 "
#echo "64 : $nl64"
#echo "80 : $nl80"
#echo "96 : $nl96"
av=`echo $nl20 $nl30 $nl36 $nl48 $nl64 $nl80 $nl96|awk '{print ($1+$2+$3+$4+$5+$6+$7)/7}'`
echo 'nonlop(apply)average= '$av

cd ../A3
#echo "projbd"
pbd20=`grep 'projbd' report |head -1|awk '{print $NF}'`
pbd30=`grep 'projbd' report |head -3|tail -1|awk '{print $NF}'`
pbd36=`grep 'projbd' report |head -5|tail -1|awk '{print $NF}'`
#echo "20 : $pbd20"
#echo "30 : $pbd30"
#echo "36 : $pbd36"
cd ../B3
pbd48=`grep 'projbd' report |head -1|awk '{print $NF}'`
pbd64=`grep 'projbd' report |head -3|tail -1|awk '{print $NF}'`
pbd80=`grep 'projbd' report |head -5|tail -1|awk '{print $NF}'`
pbd96=`grep 'projbd' report |tail -1|awk '{print $NF}'`
#echo "48 : $pbd48 "
#echo "64 : $pbd64"
#echo "80 : $pbd80"
#echo "96 : $pbd96"
av=`echo $pbd20 $pbd30 $pbd36 $pbd48 $pbd64 $pbd80 $pbd96|awk '{print ($1+$2+$3+$4+$5+$6+$7)/7}'`
echo 'projbd average= '$av
