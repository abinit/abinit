#!/bin/bash
#
version="1.0"
#
search="Unused variable"
cat /dev/null > rmstatus.log
#
fline="`egrep --line-number --max-count=1 '^make.*Entering.*src\/12_hide_mpi' make.log | awk -F':' '{print $1}'`"
#echo $fline
regexpr="1,${fline}d"
sed="sed -i.bak -e '${regexpr}' make.log"
if [ "${fline}" ];then
	echo "fline : $fline"
	eval $sed
fi
#
IFS=$'\n'
cmd="grep -B 4 '${search}' make.log"
i=0
for line in $(eval $cmd); do
	var[i]=$line
        let i=i+1
	status="echo \"$line\" | grep 'Warning: Unused variable'"
	status=`eval "${status}"`
	if [ ! "${status}" ]; then continue; fi
#        for ((a=0; a < $i ; a++)); do
#	  echo " $a: ${var[$a]}"
#	done
	let a=i-1
	let af=a-3
	let as=a-2
        i=0
	filename=`echo ${var[$af]} | awk -F'[:]' '{print $1}'`
	tmp=`echo ${var[$af]} | awk -F'[:]' '{print $2}'`
	linenumber=`echo ${tmp} | awk -F'.' '{print $1}'`
	charpos=`echo ${tmp} | awk -F'.' '{print $2}'`
	source=`echo ${var[$as]}`
	variable=`echo ${var[$a]} | awk -F"'" '{print $2}'`
	status="`grep ${variable}  src/*/${filename} | wc -l`"
	printf "$filename == ${source} == ${variable} == $status\n" >> rmstatus.log
	if [ "${status}" -gt "1" ]; then continue; fi
#
if [ "`echo \"$source\" | tr 'A-Z' 'a-z' | grep \"$variable(\"`" ]; then
        status="$source"
        cmd1="echo '$status' | sed -n -e 's/.*\($variable([^)]*)\)$/\1/Ip'"
        status2=`eval $cmd1`
        if [ -z "$status2" ]; then
                cmd1="echo '$status' | sed -n -e 's/.*\($variable([^)]*)\) *,.*/\1/Ip'"
                status2=`eval $cmd1`
        fi
        cmd1="echo '$status2' | sed -e 's/(/\\\(/g' | sed -e 's/)/\\\)/g'"
        findvar=`eval $cmd1`
else
#        echo "simple variable..."
        variable2="`echo $variable| tr 'a-z' 'A-Z'`"
        cmd1="echo '$source' | sed -n -e 's/.*[), :]*\($variable2\).*/\1/Ip'"
        findvar=`eval $cmd1`
        if [ -z "$findvar" ]; then
                cmd1="echo '$source' | sed -n -e 's/.*[), :]*\($variable\).*/\1/Ip'"
                findvar=`eval $cmd1`
        fi
     		
        cmd1="echo '$source' | sed -n -e 's/.*[), :]*\($variable=[^,]*\).*/\1/Ip'"
	test=`eval $cmd1`
	if [ -n "$test" ]; then
		findvar=$test
	fi
fi
        printf "$filename : $variable [$findvar] -> "
#
        source2=$source
while true; do

#        source2="`echo ${source2} | tr 'A-Z' 'a-z'`"
        cmd1="echo '${source2}' | egrep ', *${findvar}$'"
        position="`eval $cmd1`"
        if [ ! -z "$position" ]; then printf "'end' ";position="end";break; fi
#
        cmd2="echo '${source2}' | egrep ', *${findvar} *,'"
        posmid=$(eval $cmd2)
        if [ ! -z "$posmid" ]; then printf "'middle' ";position="middle";break; fi
#
        cmd3="echo '${source2}' | egrep '[) :]* ${findvar} *,'"
        position=$(eval $cmd3)
        if [ ! -z "$position" ]; then printf "'start'";position="start";break; fi
#
        cmd3="echo '${source2}' | egrep '[) :]* ${findvar}$'"
        position=$(eval $cmd3)
        if [ ! -z "$position" ]; then printf "'alone'";position="alone";break; fi
#
        cmd1="echo '${source2}' | egrep ',* *${findvar}=[0-9]*$'"
        position="`eval $cmd1`"
        if [ ! -z "$position" ]; then printf "'end' ";position="end";break; fi
#
        cmd1="echo '${source2}' | egrep ', *${findvar}=[0-9]* *,*'"
        position="`eval $cmd1`"
        if [ ! -z "$position" ]; then printf "'middle' ";position="middle";break; fi
#
        cmd1="echo '${source2}' | egrep '[) :]* ${findvar}=[0-9]* *,*'"
        position="`eval $cmd1`"
        if [ ! -z "$position" ]; then printf "'start' ";position="start";break; fi
#
        break
done
        printf " [${source}]\n"
#
	cmd="echo ${findvar} | sed -e 's/\\\//g'"
	findvar="`eval ${cmd}`"
	case "${position}" in
	 "start" | "middle" )
	 	regexpr="${linenumber}s/${findvar} *,//"
	 ;;
	 "end")
	 	regexpr="${linenumber}s/, *${findvar}.*//"
	 ;;
	 "alone")
	 	regexpr="${linenumber}s/^/!! removed/"
	 ;;
	 * )
		regexpr="none"
	 ;;
	esac
        sed="sed -i.bak -e '${regexpr}' src/*/${filename}"
	if [ "${regexpr}" != "none" ]; then eval $sed; fi
done
