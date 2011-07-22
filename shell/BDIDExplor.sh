#!/bin/sh

cd /bgc/data/psych/dae
list=`ls -rt` # /bgc/data/psych/dae`

for f in $list
do
    #c=`grep -c or=90 $f`
    c=`grep -c bd= $f`
#    if [ "$c" -gt "0" ] ; then
#	if [ "$c" -lt "12" ]; then
#	    echo "file is $f and it has $c"
#	fi
#    fi
    if [ "$c" -gt "70" ]; then 
	echo "file is  $f and it has $c bd"
    fi
done

cd -
