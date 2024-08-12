from check_software import get_tool_paths
import subprocess
from Bio import SeqIO
import argparse
import os
get_tool_paths('mkvtree')
get_tool_paths('vmatch')

def reverse_sequence(fasta,prefix):
    with open(prefix+'reverse.fasta','w') as f1:
        with open(fasta,'r') as f:
            body=""
            for i in f:
                if i.startswith('>'):
                    if body:
                        f1.write('>'+header+'\n'+body+'\n')
                        body = body[::-1]
                        f1.write('>reversr_'+header+'\n'+body+'\n')
                        body=""
                    header = i.strip()[1:]
                else:
                    body = body + i.strip()
            f1.write('>'+header+'\n'+body+'\n')
            body = body[::-1]
            f1.write('>reversr_'+header+'\n'+body+'\n')
    return prefix+'reverse.fasta'

def make_db(db_fa):
    cmd_mkvtree=['mkvtree', '-db',db_fa, '-dna', '-pl', '-v', '-allout']
    subprocess.run(cmd_mkvtree)

def db_repeat_find(db_fa,sl,seedlen):
    cmd_vmatch=['vmatch',  '-v', '-l', sl , '-seedlength', seedlen ,'-h', '3', '-d', '-p',  db_fa]
    result1=subprocess.Popen(cmd_vmatch,stdout=subprocess.PIPE)
    table1=[]
    for i in result1.stdout:
        if not i.startswith(b'#'):
            tmp_table=i.strip().split()
            tmp_table=[x.decode('utf8') for x in tmp_table]
            table1.append(tmp_table)
    return table1

def remove_dup(final_tab):
    len1=len(final_tab)
    dup=[]
    for i in range(len1):
        for j in range(i+1,len1):
            if final_tab[i][2] == final_tab[j][2] and \
                    final_tab[i][3] == final_tab[j][6] and \
                    final_tab[i][4] == final_tab[j][7] and \
                    final_tab[j][3] == final_tab[i][6] and \
                    final_tab[j][4] == final_tab[i][7]:
                dup.append(j)
    #get item
    dup_item=[]
    for i in dup:
        dup_item.append(final_tab[i])
    #delet
    for i in dup_item:
        final_tab.remove(i)
    return final_tab

def generate_table(db_fa_tab,fasta):
    seqs=SeqIO.parse(fasta, "fasta")
    seq_ids,seq_seqs,seq_lens=[],[],[]
    for seq_record in seqs:
        seq_ids.append(seq_record.id)
        seq_seqs.append(seq_record.seq)
        seq_lens.append(len(seq_record))
    db_fa_tab=[x for x in db_fa_tab if int(x[1])%2==0]
    final_tab=[]
    for i in db_fa_tab:
        if (int(i[1])+int(i[5]))%2==0:
            Length= i[0]
            if i[3]=='D':
                Type='F'
            else:
                Type=i[3]
            Seq_A_id=seq_ids[int(i[1])]
            Repeat_A_start=int(i[2])+1
            Repeat_A_end=int(i[2])+int(i[0])
            Seq_B_id=seq_ids[int(i[5])]
            Repeat_B_start=int(i[6])+1
            Repeat_B_end=int(i[6])+int(i[4])
            Hamming_Distance=i[7]
            E_value=i[8]
            Repeat_A_seq=seq_seqs[int(i[1])][int(i[2]):int(i[2])+int(i[0])]
            Repeat_B_seq=seq_seqs[int(i[5])][int(i[6]):int(i[6])+int(i[4])]
        else:
            Length= i[0]
            if i[3]=='D':
                Type='R'
            else:
                Type='C'
            Seq_A_id=seq_ids[int(i[1])]
            Repeat_A_start=int(i[2])+1
            Repeat_A_end=int(i[2])+int(i[0])
            Seq_B_id=seq_ids[int(i[5])-1]
            Repeat_B_start=seq_lens[int(i[5])]-int(i[6])-int(i[4])+1
            Repeat_B_end=seq_lens[int(i[5])]-int(i[6])
            Hamming_Distance=i[7]
            E_value=i[8]
            Repeat_A_seq=seq_seqs[int(i[1])][int(i[2]):int(i[2])+int(i[0])]
            Repeat_B_seq=seq_seqs[int(i[5])-1][Repeat_B_start-1:Repeat_B_start-1+int(i[4])]
        tmp1=[Length,Type,Seq_A_id,Repeat_A_start,Repeat_A_end,Seq_B_id,Repeat_B_start,
                Repeat_B_end,Hamming_Distance,E_value,Repeat_A_seq,Repeat_B_seq]
        tmp1=[str(x) for x in tmp1]
        final_tab.append(tmp1)
    #final_tab=remove_dup(final_tab)
    final_tab=sorted(final_tab,key=(lambda x:int(x[0])),reverse=True)
    count=0
    final_tab1=[]
    for i in final_tab:
        count+=1
        x=['R'+str(count)]
        x.extend(i)
        final_tab1.append(x)
    return final_tab1

def generate_output(tab,prefix):
    with open(prefix+'_long_repeat.xls','w') as f:
        header=['ID','Length','Type','Seq_A_id','Repeat_A_start','Repeat_A_end','Seq_B_id','Repeat_B_start',
                'Repeat_B_end','Hamming_Distance','E-value','Repeat_A_seq','Repeat_B_seq']
        f.write('\t'.join(header)+'\n')
        for i in tab:
            f.write('\t'.join(i)+'\n')
    return prefix+'_long_repeat.xls'

def main(fasta_file,prefix,sl,seedlen):
    quy_fa=reverse_sequence(fasta_file,prefix)
    make_db(quy_fa)
    tab1=db_repeat_find(quy_fa,sl,seedlen)
    final_tab=generate_table(tab1,quy_fa)
    generate_output(final_tab,prefix)
    for root, dirs, files in os.walk('./'):
        for name in files:
            if name.startswith(prefix+'reverse.fasta'):
                os.remove(os.path.join(root, name))
                print("Delete File: " + os.path.join(root, name))
    return {prefix:final_tab}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for get chloroplast statistics")
    parser.add_argument("fasta", type=str, help="input fasta file")
    parser.add_argument("prefix", type=str, help="output prefix")
    parser.add_argument("-l", type=str,default="30",help="search minimal length. [default: 30]")
    parser.add_argument("-seedlength", type=str,default="7",help="Specify the seed length. [default: 7]")
    args = parser.parse_args()
    fasta_file = args.fasta
    prefix = args.prefix
    searchlen=args.l
    seedlen=args.seedlength
    main(fasta_file,prefix,searchlen,seedlen)
