import sys,os
from Bio.Seq import Seq
import pyfaidx
from pymsaviz import MsaViz
from subprocess import call
import argparse


def translate(seq):
    length = len(seq)
    s = length % 3
    proseq = ""
    for i in range(0, length - s, 3):
        aa = seq[i:i+3]
        aaseq = Seq(aa)
        if "-" in aa:
            aas = "-"
        else:
            aas = aaseq.translate()
        proseq += aas
    #print(proseq)
    # if s == 0:
    #     return proseq
    # else:
    #     proseq += "-"
    return proseq, seq


def parser_fasta(alignfile, outfile):
    out = open(outfile,"w")
    fasta = pyfaidx.Fasta(alignfile)
    namelist = fasta.keys()
    for name in namelist:
        seq = str(fasta[name]).strip()
        proseq, sss= translate(seq)
        proseq_format = ""
        for aa in proseq.strip():
            proseq_format += f"-{aa}-"
        out.write(f">{name}\n")
        out.write(f"{seq}\n")
        out.write(f">{name}__AA\n")
        out.write(f"{proseq_format}\n")
    out.close()


def plot_msa(infile, outpng):
    mv = MsaViz(infile, color_scheme="Nucleotide", wrap_length=81, show_grid=False, show_consensus=False, wrap_space_size=1)
    mv.savefig(outpng, dpi=300)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plot base and aa sequences")
    parser.add_argument("-i", "--input", type=str,help="the input fasta file", dest="infile")
    parser.add_argument("-o", "--output", type=str,help="the output png file", dest = "outfile")
    args = parser.parse_args()
    infile = args.infile
    outpng = args.outfile
    tmp = "tmp.fasta"
    parser_fasta(infile, tmp)
    plot_msa(tmp, outpng)
    # try:
    #     parser_fasta(infile, tmp)
    #     plot_msa(tmp, outpng)
    # except:
    #     print("Print usage: python plotbase -h")



        
        
            


            
    
    
    




