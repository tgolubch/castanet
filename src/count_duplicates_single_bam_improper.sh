BamPath="${1}"
BamName=$(basename "${BamPath%%.bam}")
BamName=$(sed s'/_dedup//' <<< ${BamName})

samtools view -F2048 -F4 ${BamPath} | python3 src/parse_bam_positions.py ${BamName} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g
