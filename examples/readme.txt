./FastTree /Users/gc3045/PLUSplacer-taxtastic/examples2/ref.fasta > /Users/gc3045/PLUSplacer-taxtastic/examples2/bbtre.nwk

FastTree -nosupport -gtr -gamma -nt -log /Users/gc3045/PLUSplacer-taxtastic/examples2/bbtre.log -intree /Users/gc3045/PLUSplacer-taxtastic/examples2/bbtre.nwk < /Users/gc3045/PLUSplacer-taxtastic/examples2/ref.fasta > /Users/gc3045/PLUSplacer-taxtastic/examples2/reest_bbtre.nwk

python3 pplacer-tax-SCAMPP.py -i examples2/bbtre.log -t examples2/reest_bbtre.nwk -d examples2/ -o output.jplace -a examples2/aln.fasta -r examples2/ref.fasta

