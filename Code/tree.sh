tools=/exports/mm-run/tools/
fasta=$1

echo " === Align"; ${tools}/mafft-7.505/bin/mafft --auto  --quiet ${fasta}.fasta  > ${fasta}_.fasta
echo " === Tree" ; ${tools}/FastTree/FastTreeMP   -gtr -nt -quiet ${fasta}_.fasta > ${fasta}_.tree
