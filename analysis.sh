
# input data file
fastafile="data/bla_5000.fa"

# generate results folder
mkdir -p results

# build pangraph
pangraph build \
    --len 20 \
    --alpha 0 \
    --beta 10 \
    --sensitivity 5 \
    --test \
    $fastafile > results/raw_pangraph.json

# polish pangraph and refine alignments
pangraph polish -c results/raw_pangraph.json > results/pangraph.json

# evaluate mash distance
mash triangle $fastafile > results/mash.tsv
python3 scripts/mash_triangle_to_csv.py \
    --mash_tri results/mash.tsv \
    --csv results/mash.csv

# generate figures
python3 analysis.py \

# extract alignments
python3 scripts/extract_alignments.py \
    --pangraph results/pangraph.json \
    --aln_folder results/alignments
