#inputs:
# 1. Assembly gtf file
# 2. Output folder path
# 3. Unique subexon length threshold
# 4. ExonRefine binary


$4 -p $2subExon $1

grep subexon $2subExon.gtf > $2onlySubExons.gtf

awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NodeId") print $i,$(i+1) } }' $2onlySubExons.gtf | sed 's/\;//g' | grep -w NodeId | sort | uniq > $2subExonIds.txt

grep -w -o -f $2subExonIds.txt $2subExon.gtf | sort | uniq -c | awk '$1 == 1 {print $2,$3}' > $2uniqueSubExons.txt

grep -w -f $2uniqueSubExons.txt $2onlySubExons.gtf | awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "transcript_id") print $5-$4+1,$(i+1) } }' | sed 's/\"//g' | sed 's/\;//g' > $2transcripts_uniquesubExonLength.txt

awk -v threshold="$3" '$1 <= threshold {print $2}' $2transcripts_uniquesubExonLength.txt | sort | uniq > $2transcripts_uniqueSubexonLength_threshold.txt
