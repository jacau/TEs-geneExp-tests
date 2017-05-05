import gffutils
import pandas as pd
import numpy as np
import argparse

def arguments():
    parser  = argparse.ArgumentParser(description="Loads gene location, gene expression and TE population database files and outputs matrix for smile plot.")
    parser.add_argument("-g","--gene_expression", help="Path for gene expression file", required=True)
    parser.add_argument("-d","--mut_db", help="Path for TE or SNP database file", required=True)
    parser.add_argument("-o","--output", help="Path for output", required=True)
    parser.add_argument("-m","--mutation_type", help="Are the mutations TEs (te) or SNPs (snp)?", choices=['te','snp'], required=True)
    parser.add_argument("-n","--pop_num", help="Number of individuals in the population, default 124 for my current population", type=int, required=True)
    parser.add_argument("-r","--rank_num", help="Number of individuals per rank bin, set by giving an int which is a power of 2", type=int, default=1)
    parser.add_argument("-H","--het_type", help="How to count TEs: 'a' for allele, 'i' for individual (ignoring heterozygosity)", choices=['a','i'], default='a')
    parser.add_argument("-t","--te_type", help="Subset by TE type, default to no subset", choices=['NA','RNA','DNA','RNA-TRIM','DNA-MITE','LINE','SINE','LTR','DIRS','TIR','Helitron','Maverick'], default='NA')
    parser.add_argument("-i","--interval_len", help="Length of interval around genes, default 2000", type=int, default=2000)
    parser.add_argument("-s","--mutation_location", help="Consider only upstream or downstream or within", choices=['n','u','d','i'], default='n')

    args = parser.parse_args()
    return(args)

def loop_gene(gene_expr,mut_db,smile,n_rank,het,interval,up,m_type,n):

    for i in range(gene_expr.shape[0]):
        
        #only columns 4-127 (0-based) contain expression values
        #in subset below used 5:129 but after pac becomes index and no longer column
        gene = gene_expr.iloc[i,5:(5+n)].sort_values()

        pac = gene.name
        scaff_gene = gene_expr.loc[pac,'scaff']
        dir = gene_expr.loc[pac,'dir']
        
        start_interval = gene_expr.loc[pac,'start'] - interval
        end_interval = gene_expr.loc[pac,'end'] + interval
        
        if (m_loc == 'u' and dir == '+') or (m_loc == 'd' and dir == '-'):
            end_interval = start_interval + interval
        elif (m_loc == 'u' and dir == '-') or (m_loc == 'd' and dir == '+'):
            start_interval = end_interval - interval
        elif m_loc == 'i':
            start_interval = gene_expr.loc[pac,'start']
            end_interval = gene_expr.loc[pac,'end']
        
        subset = mut_db[(mut_db.scaff == scaff_gene) & (mut_db.start > start_interval) & (mut_db.start < end_interval)]
        
        if subset.shape[0] != 0:
   
            #build up a list of mutation (TE or SNP) numbers to append to smile df
            mut_list = []
        
            for j in range(0,n,n_rank):
                    
                counter = 0
                
                for x in range(n_rank):
                
                    #this is the name of the jth individual for the ith gene sorted by expression 
                    ind = gene.axes[0][j+x]
                    #number of TEs for individual j at gene i

                    for x in range(subset.shape[0]):
                        genotype = subset.iloc[x].loc[ind]
                        if het == 'a':
                            counter += genotype
                        elif (het == 'i') & ((genotype == 1) | (genotype == 2)):
                            counter += 1
                
                mut_list.append(counter)
        else:
            mut_list = [0]*(n/n_rank)

        temp_series = pd.Series(mut_list,index=smile.columns,name=pac)
        smile = smile.append(temp_series)

    return smile

#command line arguments: input unreduced (raw/manually edited) te_db file, output file name
args = arguments()
mutpath = args.mut_db
exprpath = args.gene_expression
outpath = args.output

n_rank = args.rank_num
het = args.het_type
n = args.pop_num
type = args.te_type
interval = args.interval_len
m_loc = args.mutation_location
mut_type = args.mutation_type

#load te_db as dataframe
#requires sorted te_db dataframe, remember to check whether file is comma or tab separated
mut_external = pd.read_csv(mutpath, sep='\t')

if mut_type == 'te':
    #change columns to remove the X infront of the number (output from R)
    new_columns = mut_external.columns.values
    
    for i in range(4,(4+n)):
        new_columns[i] = new_columns[i][1:]
    
    mut_external.columns = new_columns

    #subset only TEs with minor allele frequency 3% or less
    if het == 'a':
        mut_external = mut_external[(mut_external.allele_count < 8)]
    elif het == 'i':
        mut_external = mut_external[(mut_external.ind_count < 4)]

    #subset TEs by type (if type != 'NA')
    if type != 'NA':
        if type == 'RNA':
            mut_external = mut_external[(mut_external.TE.str[0] == 'R')]
        elif type == 'DNA':
            mut_external = mut_external[(mut_external.TE.str[0] == 'D')]
        elif type == 'RNA-TRIM':
            mut_external = mut_external[(mut_external.TE == 'RXX-TRIM')]
        elif type == 'DNA-MITE':
            mut_external = mut_external[(mut_external.TE == 'DXX-MITE')]
        elif type == 'LINE':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'RI')]
        elif type == 'SINE':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'RS')]
        elif type == 'LTR':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'RL')]
        elif type == 'DIRS':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'RY')]
        elif type == 'TIR':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'DT')]
        elif type == 'Helitron':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'DH')]
        elif type == 'Maverick':
            mut_external = mut_external[(mut_external.TE.str[0:2] == 'DM')]

elif mut_type == 'snp':
    #fix column headers
    new_columns = mut_external.columns.values
    new_columns[0] = 'scaff'
    
    #column name should be pos, but use start to match up with TE database
    new_columns[1] = 'start'
    new_columns[2] = 'allele_count'
    new_columns[3] = 'total_allele_num'

    for i in range(4,(4+n)):
        temp = [int(s) for s in new_columns[i].replace(']',' ').replace(':',' ').replace('_',' ').replace('s',' ').split() if s.isdigit()]
        new_columns[i] = str(temp[0])

    mut_external.columns = new_columns

#load gene expression data
gene_expr = pd.read_csv(exprpath, sep='\t')

#change columns to remove the letter at the end of the number
new_columns = gene_expr.columns.values
for i in range(6,(6+n)):
    new_columns[i] = new_columns[i][0:-1]
    
gene_expr.columns = new_columns
gene_expr = gene_expr.set_index('pac')

#initialize new data frame for the smile plot table
smile_external = pd.DataFrame(columns=range(1,(n/n_rank) + 1))

smile_external = loop_gene(gene_expr,mut_external,smile_external,n_rank,het,interval,m_loc,mut_type,n)
smile_external.to_csv(outpath, index=False, sep='\t')