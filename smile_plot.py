import gffutils
import pandas as pd
import numpy as np
import argparse

def arguments():
    parser  = argparse.ArgumentParser(description="Loads gene location, gene expression and TE population database files and outputs matrix for smile plot.")
    parser.add_argument("-g","--gene_expression", help="Path for gene expression file", required=True)
    parser.add_argument("-te","--tedb", help="Path for TE database file", required=True)
    parser.add_argument("-o","--output", help="Path for output", required=True)
    #other potential arguments: number of rank bins, allele count vs individual count, number of individuals in the population, distance to gene interval
    #... TE type subset
    parser.add_argument("-r","--rank_num", help="Number of individuals per rank bin, set by giving an int which is a power of 2", type=int, default=1)
    parser.add_argument("-H","--het_type", help="How to count TEs: 'a' for allele, 'i' for individual (ignoring heterozygosity)", choices=['a','i'], default='a')
    #parser.add_argument("-n","--pop_num", help="Number of individuals in the population, default 124 for my current population", type=int, default=124)
    #this still needs work
    parser.add_argument("-t","--te_type", help="Subset by TE type, default to no subset", choices=['NA','RNA','DNA','RNA-TRIM','DNA-MITE','LINE','SINE','LTR','DIRS','TIR','Helitron','Maverick'], default='NA')
    parser.add_argument("-i","--interval_len", help="Length of interval around genes, default 2000", type=int, default=2000)
    parser.add_argument("-s","--mutation_location", help="Consider only upstream or downstream or within", choices=['n','u','d','i'], default='n')
    
    args = parser.parse_args()
    return(args)

def loop_gene(gene_expr,te_db,smile,n_rank,het,interval,up):

    for i in range(gene_expr.shape[0]):
        
        #only columns 4-127 (0-based) contain expression values
        #in subset below used 5:129 but after pac becomes index and no longer column
        gene = gene_expr.iloc[i,5:129].sort_values()
        
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
            
        
        subset = te_db[(te_db.scaff == scaff_gene) & (te_db.start > start_interval) & (te_db.start < end_interval)]
        
        if subset.shape[0] != 0:
   
            #build up a list of TE numbers to append to smile df
            te_list = []
        
            for j in range(0,len(gene.axes[0]),n_rank):
                    
                counter = 0
                
                for x in range(n_rank):
                
                    #this is the name of the jth individual for the ith gene sorted by expression 
                    ind = gene.axes[0][j+x]
                    #number of TEs for individual j at gene i

                    for x in range(subset.shape[0]):
                        if het == 'a':
                            counter += subset.iloc[x].loc[ind]
                        elif (het == 'i') & (subset.iloc[x].loc[ind] != 0):
                            counter += 1
                
                te_list.append(counter)
        else:
            te_list = [0]*(124/n_rank)

        temp_series = pd.Series(te_list,index=smile.columns,name=pac)
        smile = smile.append(temp_series)

    return smile

#command line arguments: input unreduced (raw/manually edited) te_db file, output file name
args = arguments()
tepath = args.tedb
exprpath = args.gene_expression
outpath = args.output

n_rank = args.rank_num
het = args.het_type
#n = args.pop_num
type = args.te_type
interval = args.interval_len
m_loc = args.mutation_location

#load te_db as dataframe
#requires sorted te_db dataframe, remember to check whether file is comma or tab separated
tedf_external = pd.read_csv(tepath, sep='\t')

#change columns to remove the X infront of the number (output from R)
new_columns = tedf_external.columns.values
for i in range(4,128):
    new_columns[i] = new_columns[i][1:]
tedf_external.columns = new_columns

#subset only TEs with minor allele frequency 3% or less
if het == 'a':
    tedf_external = tedf_external[(tedf_external.allele_count < 8)]
elif het == 'i':
    tedf_external = tedf_external[(tedf_external.ind_count < 4)]

#subset TEs by type (if type != 'NA')
if type != 'NA':
    if type == 'RNA':
        tedf_external = tedf_external[(tedf_external.TE.str[0] == 'R')]
    elif type == 'DNA':
        tedf_external = tedf_external[(tedf_external.TE.str[0] == 'D')]
    elif type == 'RNA-TRIM':
        tedf_external = tedf_external[(tedf_external.TE == 'RXX-TRIM')]
    elif type == 'DNA-MITE':
        tedf_external = tedf_external[(tedf_external.TE == 'DXX-MITE')]
    elif type == 'LINE':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'RI')]
    elif type == 'SINE':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'RS')]
    elif type == 'LTR':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'RL')]
    elif type == 'DIRS':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'RY')]
    elif type == 'TIR':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'DT')]
    elif type == 'Helitron':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'DH')]
    elif type == 'Maverick':
        tedf_external = tedf_external[(tedf_external.TE.str[0:2] == 'DM')]

#load gene expression data
gene_expr = pd.read_csv(exprpath, sep='\t')

#change columns to remove the letter at the end of the number
new_columns = gene_expr.columns.values
for i in range(6,130):
    new_columns[i] = new_columns[i][0:-1]
gene_expr.columns = new_columns
gene_expr = gene_expr.set_index('pac')

#initialize new data frame for the smile plot table
smile_external = pd.DataFrame(columns=range(1,(124/n_rank) + 1))

smile_external = loop_gene(gene_expr,tedf_external,smile_external,n_rank,het,interval,m_loc)
smile_external.to_csv(outpath, index=False, sep='\t')