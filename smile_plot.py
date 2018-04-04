# Outputs a tab-delim file where columns are expression ranks, rows are genes,
# and entries are the number of variants (TEs or SNPs).

import gffutils
import pandas as pd
import numpy as np
import argparse

def arguments():
    parser  = argparse.ArgumentParser(description="Loads gene location, gene "
            "expression and TE population database files and outputs matrix "
            "for smile plot.")
    parser.add_argument("-g","--gene_expression", help="Path for gene "
            "expression file", required=True)
    parser.add_argument("-d","--mut_db", help="Path for TE or SNP database "
            "file in tab-delim format", required=True)
    parser.add_argument("-o","--output", help="Path for output", required=True)
    parser.add_argument("-m","--mutation_type", help="Are the mutations TEs "
            "(te) or SNPs (snp)?", choices=['te','snp'], required=True)
    parser.add_argument("-n","--pop_num", help="Number of individuals in the "
            "population, default 124 for my current population", type=int,
            required=True)
    parser.add_argument("-i","--interval_len", help="Length of interval "
            "around genes, default 500", type=int, default=500)
    parser.add_argument("-r","--rank_num", help="Number of individuals per "
            "rank bin, set by giving an int which is a power of 2", type=int,
            default=1)
    parser.add_argument("-H","--het_type", help="How to count TEs: 'a' for "
            "allele, 'i' for individual (ignoring heterozygosity)",
            choices=['a','i'], default='a')
    parser.add_argument("-t","--te_type", help="Subset by TE type, default no "
            "subset", choices=['NA','RNA','DNA','RNA-TRIM','DNA-MITE','LINE',
                'SINE','LTR','DIRS','TIR','Helitron','Maverick'], default='NA')
    parser.add_argument("-s","--mutation_location", help="Consider only "
            "upstream or downstream or within or in gene categories (exon, "
            "intron, 5'UTR, 3'UTR)", choices=['n','u','d','i','ie','ii',
                'i5','i3'], default='n')
    parser.add_argument("-e","--gene_breakdown", help="gff file containing "
            "location of exons and UTRs NECESSARY WHEN 'ie' 'ii' 'i5 or 'i3' "
            "SELECTED FOR MUTATION LOCATION", default='NA')
    parser.add_argument("-c","--interCNS", help="Is interCNS the "
            "gene_expression file? 'y' for yes, default is 'n' for no",
            choices=['y','n'], default='n')
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
                        sub.iloc[j].loc['start']-1, sub.iloc[j].loc['start']+1),
                    featuretype=[ft]))

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
                    limit=(sub.iloc[j].loc['scaff'], sub.iloc[j].loc['start']-1,
                        sub.iloc[j].loc['start']+1)))

            if len(temp) > 0:
                iToDrop.append(j)

    sub = sub.drop(sub.index[iToDrop])
    return sub

def loop_gene(gene_expr, mut_db, smile, n_rank, het, interval, up, m_type, n,
        st_i, gene_breakdown):

    for i in range(gene_expr.shape[0]):

        # Only columns 4-127 (0-based) contain expression values.
        # In subset below (where data frame set up) used 5:129 but after pac
        # becomes index and no longer column.
        gene = gene_expr.iloc[i,st_i:(st_i+n)].sort_values()

        pac = gene.name
        scaff_gene = gene_expr.loc[pac,'scaff']

        # If CNS is expr there is no 'dir' column
        if m_loc[0] != 'i':
            dir = gene_expr.loc[pac,'dir']

        start = gene_expr.loc[pac,'start']
        end = gene_expr.loc[pac,'end']

        if (m_loc == 'n'):
            start_interval = start - interval
            end_interval = end + interval
        elif (m_loc[0] == 'i'):
            start_interval = start
            end_interval = end
        elif (m_loc == 'u' and dir == '+') or (m_loc == 'd' and dir == '-'):
            start_interval = start - interval
            end_interval = start
        elif (m_loc == 'u' and dir == '-') or (m_loc == 'd' and dir == '+'):
            start_interval = end
            end_interval = end + interval

        # Check up/downstream variants not in annotated region of another gene
        if i!=0:
            if (gene_expr.iloc[i-1,3] > start_interval):
                start_interval = gene_expr.iloc[i-1,3]

        if i!=(gene_expr.shape[0]-1):
            if (gene_expr.iloc[i+1,2] < end_interval):
                end_interval = gene_expr.iloc[i+1,2]

        # When looking up or downstream but genes overlap, end_interval will be
        # smaller value than start_interval which is a nonsensical interval and 
        # whole gene should be disregarded (no variants looked for).
        if end_interval - start_interval <= 0:
            subset = pd.DataFrame()
        else:
            subset = mut_db[(mut_db.scaff == scaff_gene) &
                    (mut_db.start > start_interval) &
                    (mut_db.start < end_interval)]

        if (len(m_loc) == 2) & (subset.shape[0] != 0):
            subset = sub_gft(subset,gene_breakdown,m_loc,pac,m_type)

        if subset.shape[0] != 0:

            # Build up list of variant (TE or SNP) numbers to append to smile df
            mut_list = []

            for j in range(0,n,n_rank):

                counter = 0

                for x in range(n_rank):

                    # This is the name of the jth individual for the ith gene
                    # sorted by expression.
                    ind = gene.axes[0][j+x]

                    for x in range(subset.shape[0]):
                        genotype = subset.iloc[x].loc[ind]
                        if het == 'a':
                            counter += genotype
                        elif (het == 'i') & ((genotype == 1) |
                                (genotype == 2)):
                            counter += 1

                mut_list.append(counter)
        else:
            mut_list = [0]*(n/n_rank)

        temp_series = pd.Series(mut_list,index=smile.columns,name=pac)
        smile = smile.append(temp_series)

    return smile

# Command line arguments
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
gbpath = args.gene_breakdown # only when looking into gene sites (e.g. exon)
CNS = args.interCNS

# Load te_db as dataframe
# Requires sorted te_db dataframe
mut_external = pd.read_csv(mutpath, sep='\t')

if mut_type == 'te':
    # Change columns to remove the X infront of the number (output from R)
    new_columns = mut_external.columns.values

    for i in range(4,(4+n)):
        new_columns[i] = new_columns[i][1:]

    mut_external.columns = new_columns

    # Subset only TEs with minor allele frequency 3% or less
    if het == 'a':
        if n == 124:
            mut_external = mut_external[(mut_external.allele_count < 8)]
        elif n == 109:
            mut_external = mut_external[(mut_external.allele_count < 7)]
    elif het == 'i':
        mut_external = mut_external[(mut_external.ind_count < 4)]

    # Subset TEs by type (if type != 'NA')
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
    # Fix column headers
    new_columns = mut_external.columns.values
    new_columns[0] = 'scaff'

    # Column name should be pos, but use start to match up with TE database
    new_columns[1] = 'start'
    new_columns[2] = 'allele_count'
    new_columns[3] = 'total_allele_num'

    for i in range(4,(4+n)):
        temp = [int(s) for s in new_columns[i].replace(']',' '). \
                replace(':',' ').replace('_',' ').replace('s',' ').split()
                if s.isdigit()]
        new_columns[i] = str(temp[0])

    mut_external.columns = new_columns

# Load gene expression data
gene_expr = pd.read_csv(exprpath, sep='\t')
if CNS == 'y':
    st_i = 5
    m_loc = 'i' # override m_loc args to ensure only 'in' setting is used
else:
    st_i = 6

# Change columns to remove the letter at the end of the number
new_columns = gene_expr.columns.values
for i in range(st_i,(st_i+n)):
    new_columns[i] = new_columns[i][0:-1]

gene_expr.columns = new_columns
gene_expr = gene_expr.set_index('pac')
gene_expr = gene_expr.sort_values(['scaff','start'])

# If we want to find exon, UTR, intron, import the gff file
if len(m_loc) == 2:
    gene_breakdown = gffutils.FeatureDB(gbpath)
else:
    gene_breakdown = 'empty_placeholder'

# Initialize new data frame for the smile plot table
smile_external = pd.DataFrame(columns=range(1,(n/n_rank) + 1))

# Pass st_i-1 because now pac is set to index so 1 less column before expression
smile_external = loop_gene(gene_expr, mut_external, smile_external, n_rank,
        het, interval, m_loc, mut_type, n, (st_i-1), gene_breakdown)

smile_external.to_csv(outpath, index=False, sep='\t')
