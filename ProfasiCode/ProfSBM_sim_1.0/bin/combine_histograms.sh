#!/bin/bash

if [ -z "$profasi" ]; then
echo "This script needs the environment variable \"profasi\" to be set, "
echo "and point to where the ProFASi executables can be found. For instance "
echo "if you have ProFASi in your home directory under /home/you/PROFASI; "
echo "you can do this: "
echo 
echo "export profasi=/home/you/PROFASI/app/bin"
exit 1
fi

for i in $(for i in `ls n0/his_*`; do basename $i | tr '\n' ' '; done); do  
    $profasi/his1dmerge -o $i.dat n*/$i ; 
done

