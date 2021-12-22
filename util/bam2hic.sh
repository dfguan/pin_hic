bam2hic() 
{
	set -x
	local sat=$1
	local bams=${@:2}
	local base=alignment
	local bed=$base.bed
	local satn=`basename $sat .sat`	
			
	# generate .hic notice: part of the code is from the 3D-DNA visualization script
	./bam2hig_bed -s $sat $bams > $bed 2>/dev/null
	sort -k2,2nr chrom.sizes > chrom.sizes.srt	
	local totl=`awk '{sum += $2}END{print sum}' chrom.sizes.srt`		
	local scale=$(( 1 + $totl / 2100000000 ))
	
	paste -d $'\t' - - < $bed | awk -F$'\t' -vOFS=$'\t' -vscale=$scale 'NR==FNR{st[$1] = l; l += $2; next;}{if ($12 == "-") $12 = "1"; else $12 = "0"; if ($6 == "-") $6 = "1"; else $6 = "0"; if ($1 > $7) {print substr($4,1,length($4)-2),$12,"assembly",int(($8 + st[$7])/scale),"16",$6, "assembly", int(($2 + st[$1])/scale), "8", $11,$5} else { print substr($4,1,length($4)-2),$6, "assembly",int(($2 + st[$1])/scale),"8",$12,"assembly",int(($8 + st[$7])/scale),"16",$5,$11}}' chrom.sizes.srt /dev/stdin | sort --parallel=8 -S10G -k3,3d -k7,7d > $base.pre.bed
	
	local res_string="2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000"
	IFS=',' read -r -a res <<< "$res_string"
	local rLen=${#res[@]}
	local add_options=$(( res[0]/scale ))
	for (( i=1; i<$rLen; i++ ))
	do 
		    add_options=$add_options","$(( res[$i]/scale ))
	done
	add_options="-r "${add_options}

	awk -vscale=$scale '{sum += $2}END{print "assembly\t"int(sum/scale)}' chrom.sizes.srt > assembly.sizes
	
	java -Xmx70g -jar ./juicebox_tools.jar pre -q 1 ${add_options} $base.pre.bed $base.pre.hic assembly.sizes 

	# generate assembly
	satool agp $sat > $satn.agp
	cut -f1 chrom.sizes.srt | xargs -n1 -i grep -w {} $satn.agp > $satn.updated.agp
	# notice: agp2assembly.py from phasegenomics
	python3 ./agp2assembly.py $satn.updated.agp $satn.updated.assembly
}

export -f bam2hic

if [ $# -lt 2 ]
then
	echo "bam2hic <sat> <bams>"
	exit 1
fi


bam2hic $@ 


