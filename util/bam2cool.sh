bam2cool() {
	set -x
	local sat=$1
	local bams=${@:2}
	local base=alignment
	local bed=$base.bed
	local satn=`basename $sat .sat`	
	./bam2hig_bed -s $sat $bams > $bed 2>/dev/null
	
	sort -k2,2nr chrom.sizes > chrom.sizes.srt	
	
	paste -d $'\t' - - < $bed | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '01'  | sort --parallel=8 -S10G -k3,3d -k7,7d > $base.pre.bed
	cat $base.pre.bed | cooler cload pairs -0 -c1 3 -p1 4 -c2 7 -p2 8 chrom.sizes.srt:1000 - $base.pre.cool
	cooler zoomify -o $base.pre.mcool $base.pre.cool		
}


export -f bam2cool

if [ $# -lt 2 ]
then
	echo "bam2cool <sat> <bams>"
	exit 1
fi

#sat=$1
#bams="${@:2}"

bam2cool $@ 


