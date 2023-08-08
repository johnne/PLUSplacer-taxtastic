FastTree -nosupport -gtr -gamma -nt -log /Users/gc3045/PLUSplacer-taxtastic/examples/bbtre.log -intree /Users/gc3045/PLUSplacer-taxtastic/examples/bbtre.nwk < /Users/gc3045/PLUSplacer-taxtastic/examples/ref.fasta > /Users/gc3045/PLUSplacer-taxtastic/examples/reest_bbtre.nwk

python3 pplacer-tax-SCAMPP.py -i examples/bbtre.log -t examples/reest_bbtre.nwk -d examples/ -o output.jplace -a examples/aln.fasta -r examples/ref.fasta

