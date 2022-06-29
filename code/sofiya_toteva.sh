cut -d "	" -f1,8,10,12 gene_disease_opt.csv > rheumatoid_arthritis1.csv

#awk '{$4=="rheumatoid arthritis"}' rheumatoid_arthritis2.csv > rheumatoid_arthritis1.csv

sort -k 1 -r rheumatoid_arthritis1.csv > rheumatoid_arthritis.csv

echo “Number of the GPCRs:”

awk '($14 == "GPCR") {count++ } END { print count }' gene_disease_opt.csv