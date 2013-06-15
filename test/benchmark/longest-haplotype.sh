if [[ $# -ne 1 ]]; then
  echo "usage: ./longest-haplotype.sh <in.vcf>"
  exit -1
fi

if ! [[ -e $1 ]]; then
  echo "Cannot read $1"
  exit -1
fi

STATS=`grep -v '^#' $1 | awk '{print $5}' | tr ',' '
' | awk '{print length}' | sort -n | awk 'BEGIN{max=0; sum=0;}
{ values[NR]=$1; sum += $1; if ( $1 > max ) { max = $1; } }
END{
  median = (NR % 2) ? values[(NR + 1) / 2] \
                    : (values[(NR / 2)] + values[(NR / 2) + 1]) / 2.0;
  print sum" "NR" "max" "sprintf("%.1f", sum/NR)" "median
}'`

sum=`echo $STATS | cut -d' ' -f1`
num=`echo $STATS | cut -d' ' -f2`
max=`echo $STATS | cut -d' ' -f3`
mean=`echo $STATS | cut -d' ' -f4`
median=`echo $STATS | cut -d' ' -f5`

echo "[Haplotype length (bp)] longest: $max; mean: $mean; median: $median  \
[$num alleles; $sum bp total]"
