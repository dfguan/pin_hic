bam2hic() 
{
	set -x
	local sat=$1
	local bams=${@:2}
	local base=alignment
	local bed=$base.bed
	local satn=`basename $sat .sat`	
	./bam2hig_bed -s $sat $bams > $bed 2>/dev/null
	
	bash ./sat2chrom_size.sh $sat > $satn.chrom.sizes
	paste -d $'\t' - - < $bed | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($12 == "-") $12 = "1"; else $12 = "0"; if ($6 == "-") $6 = "1"; else $6 = "0"; if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else { print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | sort --parallel=8 -S10G -k3,3d -k7,7d > $base.pre.bed
	java -Xmx70g -jar ./juicebox_tools.jar pre -q 1  $base.pre.bed $base.pre.hic $satn.chrom.sizes
}

export -f bam2hic

if [ $# -lt 2 ]
then
	echo "bam2hic <sat> <bams>"
	exit 1
fi


bam2hic $@ 


