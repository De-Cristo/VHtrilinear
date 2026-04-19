#!/bin/bash - 
dir=$1
#loop over all subdirectories
id=0
olpincfile=check_olp.inc
processfile=proc_ml
rm -f $olpincfile
rm -f $processfile
#the map that translate particles into pdg ids
#This should work in bash 4.0 and above, ksh, and zsh:
declare -A pdgmap
pdgmap["g"]="21"
pdgmap["d"]="1"
pdgmap["u"]="2"
pdgmap["s"]="3"
pdgmap["c"]="4"
pdgmap["b"]="5"
pdgmap["d~"]="-1"
pdgmap["u~"]="-2"
pdgmap["s~"]="-3"
pdgmap["c~"]="-4"
pdgmap["b~"]="-5"
#loop over all subprocesses
ifct=0
for i in `find "$dir"/SubProcesses -type d -name P0*`;do
#the process that we want to called MadLoop is 
    j=`cat $i/born.f|grep Process -m1|sed -e 's/^.*Process://g' -e 's/WEIGHTED.*$//g' -e 's/$/ [virt=QED]/g'`
    if test $id == 0;then
        j=`echo $j|sed -e 's/^/generate /g'`
    else
        j=`echo $j|sed -e 's/^/add process /g'`
    fi
#several partionic processes shared one subdirectory, let's loop over it
    for k in `cat $i/born.f|sed '/IMPLICIT/q'| grep Process|sed -e 's/^.*Process://g' -e 's/WEIGHTED.*$//g'|sed -e 's/ /*/g'`;do
#output code for choice of which matrix element of MadLoop
        k=`echo $k|sed -e 's/\*/ /g'`
        echo -n "      " >> $olpincfile
        if test $ifct != 0;then
            echo -n "else ">>$olpincfile
        fi
        let "ifct=ifct+1"
        echo  "if(" >> $olpincfile
        ctl=0
        for l in $k;do
            idl=${pdgmap[$l]}
            if test x$l != x">";then
                let "ctl=ctl+1"
            fi
            if test x$idl != x;then
                echo -n '     $ ' >> $olpincfile
                if test $ctl != 1;then
                    echo -n ' .and. '>>$olpincfile
                fi
                echo 'pdg('"$ctl"', i).eq.'$idl >>$olpincfile
            fi
        done
#        ./transpdg "$k" >> check_olp.inc
        echo "     $ ) then " >> $olpincfile
        echo "       CALL ML5_"$id"_SLOOPMATRIX_THRES(P,">> $olpincfile
        echo "     $     MATELEM,tolerance,PREC_FOUND,RETURNCODE)">>$olpincfile
    done
    echo $j >> $processfile
    let "id=id+1"
done
echo "       else " >> $olpincfile
echo "        write (*,*) 'Dont know what to do:'," >> $olpincfile
echo "     $  (pdg(j,i),j=1,nexternal-1)" >> $olpincfile
echo "        stop" >> $olpincfile
echo "       endif" >> $olpincfile
