#modified from smile_plot.py
#outputs a tab-delim file containing the number and proportion of variants (TEs
#or SNPs) per outlier and non-outlier for each gene

import gffutils
import pandas as pd
import numpy as np
import argparse

def arguments():
    parser  = argparse.ArgumentParser(description="Loads expression outliers "
            "and TE population database files and outputs matrix for smile "
            "plot.")
    parser.add_argument("-g","--expression_outliers", help="Path for "
            "expression outlier file created using find_expression_outliers.py",
            required=True)
    parser.add_argument("-d","--mut_db", help="Path for TE or SNP database "
            "file in tab-delim format", required=True)
    parser.add_argument("-o","--output", help="Path for output", required=True)
    parser.add_argument("-m","--mutation_type", help="Are the mutations TEs "
            "(te) or SNPs (snp)?", choices=['te','snp'], required=True)
    parser.add_argument("-n","--pop_num", help="Number of individuals in the "
            "population, default 124 for my current population", type=int,
            required=True)
    parser.add_argument("-i","--interval_len", help="Length of interval "
            "around genes", type=int, required=True)
    parser.add_argument("-E","--interval_end", help="End of interval, used "
            "for upstream tighter intervals (e.g. 100-200 bp upstream)",
            type=int, default=0)
    parser.add_argument("-H","--het_type", help="How to count TEs: 'a' for "
            "allele, 'i' for individual (ignoring heterozygosity)",
            choices=['a','i'], default='a')
    parser.add_argument("-t","--te_type", help="Subset by TE type, default to "
            "no subset", choices=['NA', 'RNA', 'DNA', 'RNA-TRIM', 'DNA-MITE',
                'LINE', 'SINE', 'LTR', 'DIRS', 'TIR', 'Helitron', 'Maverick'],
            default='NA')
    parser.add_argument("-s","--mutation_location", help="Consider only "
            "upstream or downstream or within or in gene categories (exon, "
            "intron, 5'UTR, 3'UTR)", choices=['n','u','d','i','ie','ii',
                'i5','i3'], default='n')
    parser.add_argument("-e","--gene_breakdown", help="gff file containing "
            "location of exons and UTRs NECESSARY WHEN 'ie' 'ii' 'i5 or 'i3' "
            "SELECTED FOR MUTATION LOCATION", default='NA')
    args = parser.parse_args()
    return(args)

def sub_gft(sub,db,m_loc,pac,m_type):
    """check TE subset of within gene TEs and return further subset with
    those that are in the functional type chosen"""

    if m_loc != 'ii':
        if m_loc == 'ie':
            ft = 'CDS'
        elif m_loc == 'i3':
            ft = 'three_prime_UTR'
        elif m_loc == 'i5':
            ft = 'five_prime_UTR'
        iToDrop = []
        for j in range(sub.shape[0]):
            if m_type == 'te':
                temp = list(db.children('PAC:'+str(pac),
                    limit=(sub.iloc[j].loc['scaff'], sub.iloc[j].loc['start'],
                        sub.iloc[j].loc['end']), featuretype=[ft]))
            else:
                temp = list(db.children('PAC:'+str(pac),
                    limit=(sub.iloc[j].loc['scaff'],
                        sub.iloc[j].loc['start']-1,
                        sub.iloc[j].loc['start']+1), featuretype=[ft]))

            if len(temp) == 0:
                iToDrop.append(j)
    else:

        iToDrop = []
        for j in range(sub.shape[0]):
            if m_type == 'te':
                temp = list(db.children('PAC:'+str(pac),
                    limit=(sub.iloc[j].loc['scaff'], sub.iloc[j].loc['start'],
                        sub.iloc[j].loc['end'])))
            else:
                temp = list(db.children('PAC:'+str(pac),
                    limit=(sub.iloc[j].loc['scaff'],
                        sub.iloc[j].loc['start']-1,
                        sub.iloc[j].loc['start']+1)))

            if len(temp) > 0:
                iToDrop.append(j)

    sub = sub.drop(sub.index[iToDrop])
    return sub

def find_rares(geneSub,mutSub):

    mutCounter = 0

    for j in range(geneSub.size):

        #this is the name of the jth individual for the ith gene
        ind = geneSub.axes[0][j]

        #number of TEs for individual j at gene i
        for x in range(mutSub.shape[0]):

            genotype = mutSub.iloc[x].loc[ind]
            if het == 'a':
                mutCounter += genotype
            elif (het == 'i') & ((genotype == 1) | (genotype == 2)):
                mutCounter += 1

    return mutCounter

#CHANGES: changed interval_len arg to no default
def loop_gene(exprOutlier, mut_db, outdf, het, interval, intervalEnd, m_loc,
        m_type, n, gene_breakdown):

    for i in range(exprOutlier.shape[0]):

        #only columns 5:n (0-based) contain expression values
        #in subset below (where data frame set up) used 6:n 
        #but after pac becomes index and no longer column
        gene = exprOutlier.iloc[i,5:(5+n)]
        gene_out = gene[gene == 1]
        if (gene_out.size == 0):
            continue
        gene_non = gene[gene == 0]

        #solved error with in gene types, pac had .0 at the end which isn't
        #included in gene_breakdown file, so temp was always coming up empty
        pac = int(gene.name)
        scaff_gene = exprOutlier.loc[pac,'scaff']
        dir = exprOutlier.loc[pac,'dir']

        start = exprOutlier.loc[pac,'start']
        end = exprOutlier.loc[pac,'end']

        if (m_loc == 'n'):
            start_interval = start - interval
            end_interval = end + interval
        elif (m_loc[0] == 'i'):
            start_interval = start
            end_interval = end
        elif (m_loc == 'u' and dir == '+') or (m_loc == 'd' and dir == '-'):
            start_interval = start - interval
            end_interval = start - intervalEnd #if intervalEnd not set, will be 0 or equivalent to normal (gene start)
        elif (m_loc == 'u' and dir == '-') or (m_loc == 'd' and dir == '+'):
            start_interval = end + intervalEnd
            end_interval = end + interval

        #check up/downstream variants not in annotated region of another gene
        if i!=0:
            if (exprOutlier.iloc[i-1,3] > start_interval):
                start_interval = exprOutlier.iloc[i-1,3]

        if i!=(exprOutlier.shape[0]-1):
            if (exprOutlier.iloc[i+1,2] < end_interval):
                end_interval = exprOutlier.iloc[i+1,2]

        #when looking up or downstream but genes overlap, end_interval will be
        #smaller value than start_interval which is a nonsensical interval and 
        #whole gene should be disregarded (no variants looked for)
        if end_interval - start_interval <= 0:
            subset = pd.DataFrame()
        else:
            subset = mut_db[(mut_db.scaff == scaff_gene) &
                    (mut_db.start > start_interval) &
                    (mut_db.start < end_interval)]

        if (len(m_loc) == 2) & (subset.shape[0] != 0):
            subset = sub_gft(subset, gene_breakdown, m_loc, pac, m_type)

        if subset.shape[0] != 0:

            #number of rare variants near expression outliers
            outCount = find_rares(gene_out,subset)
            #number of rare variants near non-outliers
            nonCount = find_rares(gene_non,subset)
        else:
            outCount = 0
            nonCount = 0
        outProp = float(outCount) / gene_out.size
        nonProp = float(nonCount) / gene_non.size

        temp_list = exprOutlier.iloc[i,0:4].tolist()
        temp_list.extend((outCount,gene_out.size,outProp,nonCount,gene_non.size,nonProp))
        temp_series = pd.Series(temp_list,index=outdf.columns,name=pac)
        outdf = outdf.append(temp_series)

    return outdf

#command line arguments
args = arguments()
mutpath = args.mut_db
exprpath = args.expression_outliers
outpath = args.output

het = args.het_type
n = args.pop_num
type = args.te_type
interval = args.interval_len
intervalEnd = args.interval_end
m_loc = args.mutation_location
mut_type = args.mutation_type
gbpath = args.gene_breakdown #only when looking into gene sites (e.g. exon)

#load te_db as dataframe
mut_external = pd.read_csv(mutpath, sep='\t')

if mut_type == 'te':
    #change columns to remove the X infront of the number (output from R)
    new_columns = mut_external.columns.values

    for i in range(4,(4+n)):
        new_columns[i] = new_columns[i][1:]

    mut_external.columns = new_columns

    #subset only TEs with minor allele frequency 3% or less
    if het == 'a':
        if n == 124:
            mut_external = mut_external[(mut_external.allele_count < 8)]
        elif n == 109:
            mut_external = mut_external[(mut_external.allele_count < 7)]
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
        temp = [int(s) for s in new_columns[i].replace(']',' '). \
                replace(':',' ').replace('_',' ').replace('s',' ').split()
                if s.isdigit()]
        new_columns[i] = str(temp[0])

    mut_external.columns = new_columns

#load gene expression data
exprOutlier = pd.read_csv(exprpath, sep='\t')

#change columns to remove the letter at the end of the number
new_columns = exprOutlier.columns.values
for i in range(6,(6+n)):
    new_columns[i] = new_columns[i][0:-1]

exprOutlier.columns = new_columns
exprOutlier = exprOutlier.set_index('pac')
exprOutlier = exprOutlier.sort_values(['scaff','start'])

#if we want to find exon, UTR, intron, import the gff file
if len(m_loc) == 2:
    gene_breakdown = gffutils.FeatureDB(gbpath)
else:
    gene_breakdown = 'empty_placeholder'

#initialize new data frame for the smile plot table
outcols = exprOutlier.columns.values[0:4].tolist()
outcols.extend(("outlier_TEnum","o_ind_num","oTEprop","nonOutlier_TEnum","nO_ind_num","nTEprop"))
output_external = pd.DataFrame(columns=outcols)

output_external = loop_gene(exprOutlier, mut_external, output_external, het,
        interval, intervalEnd, m_loc, mut_type, n, gene_breakdown)

output_external.to_csv(outpath, sep='\t')
