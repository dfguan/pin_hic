bam2cool() {
	local bams=$1
	local sat=$2
	local ref=$3
	./build_hig $sat $bams > alignment.bed 2>/dev/null
	local bed=aligment.bed
	samtools faidx $ref
	local base=`basename $bed .bed`
	cut -f1,2 $ref.fai | sed 's/-/_/g' > $ref.chrom.sizes
	paste -d $'\t' - - < $bed | sed 's/-/_/g' | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '01'  | sort --parallel=8 -S10G -k3,3d -k7,7d > $base.pre.bed
	cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 $ref.chrom.sizes:1000 $base.pre.bed $base.pre.cool
	cooler zoomify -o $base.pre.mcool $base.pre.cool		
}


export -f bam2cool

if [ $# -lt 3 ]
then
	echo "bam2cool <bams> <sat> <ref>"
	exit 1
fi

bams="$1"
sat=$2
ref=$3

bam2cool "$bams"  $sat $ref  


