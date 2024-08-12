import RepFind_chl_mito
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def count1(input_list):
    rep_types=['F','R','C','P']
    result=[]
    for i in range(4):
        coun=0
        for j in rep_types:
            tmp=len([x for x in input_list if -int(x[9]) == i and x[2]==j])
            result.append(tmp)
            coun+=tmp
        result.append(coun)
    return result
def count2(input_key,input_list):
    coun_F=len([x for x in input_list if x[2]=='F'])
    coun_R=len([x for x in input_list if x[2]=='R'])
    coun_C=len([x for x in input_list if x[2]=='C'])
    coun_P=len([x for x in input_list if x[2]=='P'])
    return {input_key:[coun_F,coun_R,coun_C,coun_P]}
def count3(input_key,input_list):
    input_list=[int(x[1]) for x in input_list]
    coun_30_35=len([x for x in input_list if 30<=x and x<35])
    coun_35_40=len([x for x in input_list if 35<=x and x<40])
    coun_40_45=len([x for x in input_list if 40<=x and x<45])
    coun_45_50=len([x for x in input_list if 45<=x and x<50])
    coun_50_55=len([x for x in input_list if 50<=x and x<55])
    coun_55_60=len([x for x in input_list if 55<=x and x<60])
    coun_60_65=len([x for x in input_list if 60<=x and x<65])
    coun_65_70=len([x for x in input_list if 65<=x and x<70])
    coun_70=len([x for x in input_list if 70<=x])
    return {input_key:[coun_30_35,coun_35_40,coun_40_45,coun_45_50,coun_50_55,coun_55_60,coun_60_65,
                       coun_65_70,coun_70]}
def plot1(input_list,prefix):
    df=pd.DataFrame(input_list,index=['F','R','C','P'])
    #df.loc['sum'] = df.apply(lambda x: x.sum())
    pdf = PdfPages(prefix+'_repeat_type.pdf')
    ax=df.plot(kind='bar', stacked=True, figsize=(12,8))
    #ax=df[4:].T.plot(kind='line',figsize=(10,5))
    #df[:4].T.plot(kind='bar',ax=ax)
    plt.xlabel('Type of repeats')
    plt.ylabel('Number of repeats')
    plt.legend(bbox_to_anchor=(1.04, 0), loc=3, frameon=False, markerscale=0.3, fontsize=10, borderaxespad=0)
    #plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    #plt.xticks(rotation=0)
    #plt.subplots_adjust(right = 0.88)
    plt.subplots_adjust(right = 0.8,bottom = 0.18)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pdf.savefig()
    plt.close()
    pdf.close()
    pdf = PdfPages(prefix+'_repeat_type_stack.pdf')
    df[:4].T.plot(kind='bar',stacked=True,figsize=(10,5))
    plt.ylabel('Number of repeats')
    plt.legend(bbox_to_anchor=(1.04, 0), loc=3, frameon=False, markerscale=0.3, fontsize=10, borderaxespad=0)
    #plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0)
    plt.subplots_adjust(right = 0.8)
    plt.xticks(rotation=0)
    pdf.savefig()
    plt.close()
    pdf.close()
def plot2(input_list,prefix):
    df=pd.DataFrame(input_list,index=['30~34','35~39','40~44','45~49','50~54','55~59','60~64','65~69','>=70'])
    pdf = PdfPages(prefix+'_repeat_length.pdf')
    ax=df.plot(kind='bar', stacked=True, figsize=(12,8))
    plt.xlabel('Length of repeats')
    plt.ylabel('Number of repeats')
    plt.legend(bbox_to_anchor=(1.04, 0), loc=3, frameon=False, markerscale=0.4, fontsize=10, borderaxespad=0)
    plt.subplots_adjust(right = 0.8,bottom = 0.18)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pdf.savefig()
    plt.close()
    pdf.close()
    pdf = PdfPages(prefix+'_repeat_length_spec.pdf')
    ax=df.T.plot(kind='bar',stacked=True,figsize=(10,5))
    plt.ylabel('Number of repeats')
    plt.xticks(rotation=0)
    plt.legend(bbox_to_anchor=(1.05, 0), loc=3, frameon=False, markerscale=0.4, fontsize=10, borderaxespad=0)
    plt.subplots_adjust(right = 0.8,bottom = 0.18)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    pdf.savefig()
    plt.close()
    pdf.close()

def output_tab(list1,prefix):
    with open(prefix+'_rep_stat.xls','w') as f:
        header1=[' ','Hamming Distance = 0',' ',' ',' ',' ','Hamming Distance = 1',' ',' ',' ',' ',
                'Hamming Distance = 2',' ',' ',' ',' ','Hamming Distance = 3',' ',' ',' ',' ']
        header2=['id','F','R','C','P','all','F','R','C','P','all',
                 'F','R','C','P','all','F','R','C','P','all',]
        f.write('\t'.join(header1)+'\n')
        f.write('\t'.join(header2)+'\n')
        df=pd.DataFrame(list1)
        df.set_index([0], inplace=True)
        df.loc['sum'] = df.apply(lambda x: x.sum())
        for index, row in df.iterrows():
            line=list(row)
            line=[str(x) for x in line]
            f.write(str(index)+'\t'+'\t'.join(line)+'\n')

def main(parser):
    args = parser.parse_args()
    fasta_file_list = args.fasta_list.split(',')
    prefix = args.prefix
    searchlen=args.l
    seedlen=args.seedlength
    final_list=[]
    list1,list2={},{}
    for i in fasta_file_list:
        tmp=RepFind_chl_mito.main(i,i,searchlen,seedlen)
        for j in tmp.keys():
            numb=count1(tmp[j])
            numb.insert(0, j)
            final_list.append(numb)
            list1.update(count2(j,tmp[j]))
            list2.update(count3(j,tmp[j]))
    output_tab(final_list,prefix)
    plot1(list1,prefix)
    plot2(list2,prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for get chloroplast statistics")
    parser.add_argument("fasta_list", type=str, help="input fasta file list,separate with ,")
    parser.add_argument("prefix", type=str, help="output prefix")
    parser.add_argument("-l", type=str,default="30",help="search minimal length. [default: 30]")
    parser.add_argument("-seedlength", type=str,default="7",help="Specify the seed length. [default: 7]")
    main(parser)
