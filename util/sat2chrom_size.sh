sat=$1
awk -vOFS=$'\t' '{if ($1=="S" || $1 == "P") s[$2] = $3; else if ($1 == "A") a[$2] = $NF; else if ($1 == "C") c=$2;}END{split(a[c], b, ","); for (i in b) print b[i],s[b[i]]}' $sat 
