#!/bin/bash

# inputs
in_fasta="$1"
out_dir="$2"

# resources
CPU="$3"
MEM="$4"
DB=$(echo $5 | sed "s|'||g") # -this needs to contain '-d path_to_db1 -d path_to_db2'

echo $CPU $MEM
echo $(pwd)
# setup hhblits
export HHLIB=<path/to/your/hhsuite>
export PATH=$HHLIB/bin:$PATH
#HHBLITS="hhblits -o /dev/null -n 1 -mact 0.35 -maxfilt 20000 -neffmax 20 -cpu $CPU -nodiff -realign_max 20000 -maxmem $MEM -n 4 $DB"
HHBLITS="hhblits -o /dev/null -n 1 -mact 0.35 -maxfilt 20000 -neffmax 20 -all -cpu $CPU -realign_max 20000 -maxmem $MEM -n 4 $DB"

mkdir -p $out_dir/hhblits
tmp_dir="$out_dir/hhblits"
out_prefix="$out_dir/t000_"

rm -f ${out_prefix}.msa0.a3m

# perform iterative searches
prev_a3m="$in_fasta"
#for e in 1e-80 1e-70 1e-60 1e-50 1e-40 1e-30 1e-20 1e-10 1e-8 1e-6 1e-4 1e-3 1e-1
#for e in 1e-30 1e-10 1e-6 1e-3 1e-1   # -dk speed up for cameo testing
for e in 1e-50 1e-40 1e-35 1e-30 1e-25 1e-20 1e-15 1e-10 1e-8 1e-6 1e-4 1e-3 # -minkyung speed up for cameo testing. it'll use representative MSA (msa0) only (Nov 12, 2020)
do
    echo $e >> $out_dir/msa0.a3m.log
    
    if [ -f "$tmp_dir/t000_.$e.a3m" ]; then
        echo "$tmp_dir/t000_.$e.a3m already exists. Skipping hhblits run."
    else
        $HHBLITS -i $prev_a3m -oa3m $tmp_dir/t000_.$e.a3m -e $e -v 0
    fi
    hhfilter -id 90 -cov 75 -maxseq 10000 -i $tmp_dir/t000_.$e.a3m -o $tmp_dir/t000_.$e.id90cov75.a3m
    hhfilter -id 90 -cov 50 -maxseq 10000 -i $tmp_dir/t000_.$e.a3m -o $tmp_dir/t000_.$e.id90cov50.a3m
    
    prev_a3m="$tmp_dir/t000_.$e.id90cov50.a3m"
    n75=`grep -c "^>" $tmp_dir/t000_.$e.id90cov75.a3m`
    n50=`grep -c "^>" $tmp_dir/t000_.$e.id90cov50.a3m`

    echo "At $e found $n75 sequences at cov 75 id 90"
    echo "At $e found $n50 sequences at cov 50 id 90"

    if ((n75>2000)) 
    then
        if [ ! -s ${out_prefix}.msa0.a3m ]
        then
	   echo "Copying the $e,75 alignment and breaking"
	   cp $tmp_dir/t000_.$e.id90cov75.a3m ${out_prefix}.msa0.a3m 
	   break # -minkyung add this (Nov 12, 2020)
        fi
    elif ((n50>5000)) 
    then
        if [ ! -s ${out_prefix}.msa0.a3m ]
        then
	    echo "Copying the $e,50 alignment and breaking"
	    cp $tmp_dir/t000_.$e.id90cov50.a3m ${out_prefix}.msa0.a3m
	    break # -minkyung add this (Nov 12, 2020)
        fi
    else
        continue
    fi

done

if [ ! -s ${out_prefix}.msa0.a3m ]
then
    echo "Copying the 1e-3 alignment and ending"
    cp $tmp_dir/t000_.1e-3.id90cov50.a3m ${out_prefix}.msa0.a3m
fi
#echo "Going to remove " $tmp_dir
#rm -r $tmp_dir

