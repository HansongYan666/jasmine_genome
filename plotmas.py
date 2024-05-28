import sys,os
from Bio.Seq import Seq
import pyfaidx
from pymsaviz import MsaViz
from subprocess import call



def muscle_align(infasta, outfasta):
    fasta = pyfaidx.Fasta(infasta)
    out = open("tmp1.fasta", "w")
    for name in fasta.keys():
        seq = str(fasta[name])
        spe = seq[0:3].upper()
        spe = Seq(spe)
        spe = spe.reverse_complement()
        if spe in ["TAG", "TAA", "TGA"]:
            seq = Seq(seq).reverse_complement()
            seq = str(seq).strip()
            print(spe)
        out.write(f">{name}\n")
        out.write(f"{seq}\n")
    out.close()
    cmd = f"muscle -align tmp1.fasta -output {outfasta}\n"
    cmd += "#rm tmp1.fasta"
    call(cmd, shell=True)

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
    print(proseq)
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
    mv = MsaViz(infile, color_scheme="Taylor", wrap_length=85, show_grid=True, show_consensus=True)
    mv.savefig(outpng)



if __name__ == "__main__":
    infile = sys.argv[1]
    outpng = sys.argv[2]
    tmp = "tmp.fasta"
    alignfast = "align.fasta"
    muscle_align(infile, alignfast)
    parser_fasta(alignfast, tmp)
    plot_msa(tmp, outpng)



        
        
            


            
    
    
    




