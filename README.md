# Pins (拼接)

Scaffolding tool based on Hi-C reads. 

## Overview

Pin\_hic a scaffolder using Hi-C data. It applies a dual selection and local optimal strategy to bridge two contigs and output a SAT file for each iteration, the SAT format is the extension of GFA format which is able to record the scaffolding process, and can also be useful for further genomic analysis.


## Dependencies

1. zlib


## Installation
Run the following commands to install pins:

```
git clone https://github.com/dfguan/pins.git
cd pins/src && make

```

## Usage

### Scaffolding with Hi-C reads
##### Hi-C Read preprocessing
Given a list **hiclist** of Hi-C read files (suppose in fastq.gz format, paired files in a line) and the assembly **asm**, use the following code to generate Hi-C alignment files. 

```
bwa index $asm
while read -r r1 r2
do
	prefix=`basename $r1 .fastq.gz`
	bwa mem -SP -B10 -t12 $asm $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -b - > $prefix.bam
done < $hiclist
```

##### Hi-C scaffolding


Given Hi-C reads alignment **bams**, a draft assembly **asm** and a output directory **outdir**, if you want to build scaffols with Hi-C in **N** (default: 3) rounds, please try the following commands. The final assembly will be named as **scaffols_final.fa**.

```
samtools faidx $asm 
./bin/pin_hic_it -i $N -c $asm.fai -x $ref -O $outdir $bam1 $bam2 $bam3 ... 
```

Or you want to build scaffolds step by step:
##### Step 1. contact matrix calculation
From a draft assembly **asm**：

```
samtools faidx $asm
./bin/pin_hic link $bam1 $bam2 $bam3 ... > link.matrix  # this will calcuate contact numbers between any pairs of contigs.
```

From a **sat** file：

```
./bin/pin_hic link -s $sat $bam1 $bam2 $bam3 ... > link.matrix  # this will calcuate contact numbers between any pairs of contigs.
```

##### Step 2. Scaffolding graph construction
From a draft assembly **asm**:

```
/bin/pin_hic build -w100 -k3 -c $asm.fai link.matrix > scaffolds.sat # this will generate scaffolding paths. 
```

From a **sat** file:

```
/bin/pin_hic build -w100 -k3 -s $sat link.matrix > scaffolds.sat # this will generate scaffolding paths. 
```

##### Step 3. Mis-join detection
Given a **sat** file:

```
./bin/pin_hic break $sat $bam1 $bam2 $bam3 ... > scaffs.bk.sat
./bin/pin_hic gets -c $asm scaffs.bk.sat > scaffols_final.fa # get scaffold sequences.
```

A **scaffolding pipeline** of 3 iterations:

```
samtools faidx $asm
for i in `seq 1 3`
do
	if [ $i -eq 1 ]
	then 
		./bin/pin_hic link $bam1 $bam2 $bam3 ... > links_$i.matrix
		./bin/pin_hic build -w100 -k3 -c $asm.fai links_$i.matrix > scaffolds_$i.sat
	else
		./bin/pin_hic link -s scaffolds_$pi.sat $bam1 $bam2 $bam3 ... > links_$i.matrix
		./bin/pin_hic build -w100 -k3 -s scaffolds_$pi.sat links_$i.matrix > scaffolds_$i.sat 
	fi
	pi=i
done
./bin/pin_hic break -s scaffolds_$i.sat $bam1 $bam2 $bam3 ... > scaffolds_bk.sat 
./bin/pin_hic gets -c $asm scaffs.bk.sat > scaffols_final.fa 
```

### Output format: SAT (V 0.1)
SAT format is extended from the [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

#### Record types
| Tag | Description | Comment | 
| --- | --- | --- | 
| H   | Header | optional |
| S   | Sequence |  required | 
| L   | Link |  optional | 
| P   | Path |  optional | 
| A   | Scaffold set |  optional | 
| C   | Current scaffold set | optional | 


#### `H` Header
| Col | Field | Regexp | Description | Comment | 
| --- | --- | --- | --- | --- | 
| 1 | TAG | `H` | Tag | Required | 
| 2 | VER | `VN:Z:[0-9]\.[0-9]` | Version |  Required | 

#### `S` Sequence 
| Col | Field | Regexp | Description | Comment | 
| --- | --- | --- | --- | --- | 
| 1 | TAG | `S` | Tag | Required | 
| 2 | SNAME | `.+` | Sequence name | Required, primary key | 
| 3 | SLEN | `[0-9]+` | Sequence length | Required | 
| 4 | SEQ | `\*\|[A-Za-z]+` | Sequence | Required | 

#### `L` Link
| Col | Field | Regexp | Description | Comment | 
| --- | --- | --- | --- | --- | 
| 1 | TAG | `P` | Tag | Required | 
| 2 | SRCS | `.+` | Source sequence name | Required, foregin key S:SNAME | 
| 3 | SRCE | `[-+]` | Source end | Required, `+` for 5' end and `-` for 3' |
| 4 | TGTS | `.+` | Target sequence name | Required, foregin key S:SNAME | 
| 5 | TGTE |  `[-+]` | Target end | Required, `+` for 5' end and `-` for 3' |
| 6 | WGT | `wt:f:[0-9]*\.?[0-9]+` | Link weight | Optional |

#### `P` Path
| Col | Field | Regexp | Description | Comment | 
| --- | --- |  --- | --- | --- | 
| 1 | TAG | `P` | Tag | Required | 
| 2 | PNAME | `[cu][0-9]{9}` | Path name | Required, primary key | 
| 3 | PLEN | `[0-9]+` | Path length | Required | 
| 4 | NAMEL | `((.+[-+],)*(.+[-+]))\|((u[0-9]{9}[-+],)*u[0-9]{9}[-+])` | List of sequence names or path names | Required, foregin keys S:SNAME | 

#### `A` Scaffold set (or assembly set ?)
| Col | Field | Regexp | Description | Comment | 
| --- | --- |  --- | --- | --- | 
| 1 | TAG | `A` | Tag | Required | 
| 2 | ANAME | `a[0-9]{5}` | Scaffold set name | Required | 
| 3 | PNAMEL | `([cu][0-9]{9},)*[cu][0-9]{9}` | List of path names | Required, foregin keys P:PNAME| 

##### `C` Current scaffold set
| Col | Field | Regexp | Description | Comment | 
| --- | --- |  --- | --- | --- | 
| 1 | TAG | `C` | Tag | Required | 
| 2 | CNAME | `a[0-9]{5}` | Current scaffold set name | Required, foregin key A:ANAME

#### Example 

```
H	VN:Z:0.1
S	LR132056.1.4	138023	*
S	LR132056.1.5	1128790	*
S	LR132056.1.6	4496575	*
P	u000000004	662215	LR132053.1.4+,LR132053.1.5+,LR132053.1.6+
L	LR132051.1.5	+	LR132051.1.4	-	wt:f:0.028248
L	LR132051.1.6	+	LR132051.1.5	-	wt:f:0.009367
A	a00000	1	u000000004
C	a00000
```

## Limitation 


## FAQ


## Contact

Every one is Wellcomed to use and distribute the package. Bug report or any other suggestions, please use the github webpage or email me dfguan9@gmail.com. 
