#!/bin/sh
SOURCE=$1
THREAD_NUM=$2
BIN=$3
JROOT=$4
KMER_LENGTH=$5
COM_FLAG=$6
QUAL_FALG=$7
QUAL_CHAR=$8
COMPRESS_FLAG=$9
#4G -c 40 -c 6
LOG_TIME1=`date +%H:%M:%S`

#echo "f_flag:" $COMPRESS_FLAG $JROOT $COM_FLAG

if [ $COMPRESS_FLAG -eq 0 ] 
then

if [ $QUAL_FALG -eq 1 ]
then

if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
$JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM $SOURCE/*
else
$JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM $SOURCE
fi

else

if [ -d $SOURCE ]
then
$JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM $SOURCE/*
else
$JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM $SOURCE
fi

fi

else


if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
$JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM $SOURCE/*
else
$JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM $SOURCE
fi

else

if [ -d $SOURCE ]
then
$JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM $SOURCE/*
else
$JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM $SOURCE
fi

fi

fi




elif [ $COMPRESS_FLAG -eq 1 ] 
then




if [ $QUAL_FALG -eq 1 ]
then

if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
zcat $SOURCE/* | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
else
zcat $SOURCE | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
fi

else

if [ -d $SOURCE ]
then
zcat $SOURCE/* | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
else
zcat $SOURCE | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
fi

fi

else


if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
zcat $SOURCE/* | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM /dev/stdin
else
zcat $SOURCE | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM /dev/stdin
fi

else

if [ -d $SOURCE ]
then
zcat $SOURCE/* | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM /dev/stdin
else
zcat $SOURCE | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM /dev/stdin
fi

fi

fi





else







if [ $QUAL_FALG -eq 1 ]
then

if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
fastq-dump $SOURCE/* | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
else
fastq-dump $SOURCE | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
fi

else

if [ -d $SOURCE ]
then
fastq-dump $SOURCE/* | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
else
fastq-dump $SOURCE | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -Q $QUAL_CHAR -t $THREAD_NUM /dev/stdin
fi

fi

else


if [ $COM_FLAG -eq 1 ]
then

if [ -d $SOURCE ]
then
fastq-dump $SOURCE/* | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM /dev/stdin
else
fastq-dump $SOURCE | $JROOT/bin/jellyfish count -C -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 1G -t $THREAD_NUM /dev/stdin
fi

else

if [ -d $SOURCE ]
then
fastq-dump $SOURCE/* | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM /dev/stdin
else
fastq-dump $SOURCE | $JROOT/bin/jellyfish count -m $KMER_LENGTH -o $BIN/kmerInfo -c 7 -s 500M -t $THREAD_NUM /dev/stdin
fi

fi

fi



fi

#31.5% 21%(-m 31 -c 7 -s 1G)

BNUM=`ls -l $BIN/kmerInfo_* | wc -l` 
echo "jellyfish hash table number:" $BNUM


:<<BLOCK

if [ $BNUM -gt 1 ]
then
$JROOT/bin/jellyfish merge -o $BIN/kmerInfo $BIN/kmerInfo\_*
echo "Finish Jellyfish merge"
rm $BIN/kmerInfo_*
else
echo "No Jellyfish merge"
mv $BIN/kmerInfo_0 $BIN/kmerInfo
fi
BLOCK


LOG_TIME2=`date +%H:%M:%S`
echo "kmercounting time: "from $LOG_TIME1 to $LOG_TIME2
#$JROOT/bin/jellyfish dump -c -t -o $BIN/out $BIN/output
#LOG_TIME3=`date +%H:%M:%S`
#echo "dump time (txt transfer): "from $LOG_TIME2 to $LOG_TIME3
#rm $BIN/output
