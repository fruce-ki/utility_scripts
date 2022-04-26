bedtools unionbedg -i $1 $2 -filler 0 | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $0"\t"sum/(NF-4+1); }'
