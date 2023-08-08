import argparse
import pathlib

import pypangraph as pp

from Bio import SeqIO, Seq


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Extract all of the block alignments in the specified folder"""
    )
    parser.add_argument("--pangraph", help="input pangraph", type=str)
    parser.add_argument("--aln_folder", help="output alignment folder", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # generate folder
    aln_fld = pathlib.Path(args.aln_folder)
    aln_fld.mkdir(parents=True, exist_ok=True)

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)

    for b in pan.blocks:
        bid = b.id
        aln, occs = b.alignment.generate_alignments()

        recs = []
        for a, o in zip(aln, occs):
            seq = Seq.Seq(a)
            fa_id, n_occ, strand = o
            strand = "+" if strand else "-"
            rec_id = f"{fa_id}|{n_occ}|{strand}"
            rec = SeqIO.SeqRecord(seq, id=rec_id, description="")
            recs.append(rec)

        fname = aln_fld / f"{bid}.fa"
        SeqIO.write(recs, fname, "fasta")
