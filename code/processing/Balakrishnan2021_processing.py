#%%
import numpy as np 
import pandas as pd 
data = pd.read_csv('../../data/literature/Balakrishnan2021_compiled_reads_raw.csv')

# Load and tidy the raw data format
data.rename(columns = {'gene name': 'gene_name',
                       'Pre-shift':'preshift',
                       'Post-shift':'postshift',
                       '5 minutes': 5,
                       '15 minutes': 15,
                       '60 minutes': 60,
                       '120 minutes': 120}, 
            inplace=True)
data = data.melt(id_vars=['gene_name'], var_name='sample', value_name='counts')


# Compute the counts of each gene in each sample
frac_df = pd.DataFrame([])
for g, d in data.groupby(['sample']):
    d['fraction'] = d['counts'] / d['counts'].sum()
    frac_df = pd.concat([frac_df, d], sort=False)

# Save the tidied form to disk
frac_df.to_csv('../../data/literature/Balakrishnan2021_compiled_reads_tidy.csv', 
               index=False)


#%%
# Classify each gene as a ribosomal protein or not
ribosomal_genes = ['rpsU', 'rpsT', 'rpsS', 'rpsR', 'rpsQ', 'rpsP' 'rpsO',
'rpsN', 'rpsM', 'rpsL', 'rpsK', 'rpsJ', 'rpsI', 'rpsH', 'rpsG', 'rpsF', 'rpsG',
'rpsF', 'rpsE', 'rpsD', 'rpsC', 'rpsB', 'rpsA', 'rrsA', 'rplL', 'rplJ' 'rpmJ',
'rpmI', 'rpmH', 'rpmG', 'rpmF', 'rpmE' 'rpmD', 'rpmC', 'rpmB', 'rpmA', 'rplY',
'rplX', 'rplW', 'rplV', 'rplU', 'rplT', 'rplS', 'rplR', 'rplQ', 'rplP', 'rplO',
'rplN', 'rplM', 'rplK', 'rplI', 'rplF', 'rplE', 'rplD', 'rplC', 'rplB', 'rplA',
'rrfA', 'rrlA']

coarse_grained_df = pd.DataFrame([])
for g, d in frac_df.groupby('sample'):
    _d = d[d['gene_name'].isin(ribosomal_genes)]
    phiRb = _d['fraction'].sum()
    phiO = 0.55
    phiMb = 1 - phiRb - phiO
    _df = pd.DataFrame({'sample': g, 'phiRb': phiRb, 'phiO': phiO, 'phiMb': phiMb}, index=[0])
    coarse_grained_df = pd.concat([coarse_grained_df, _df], sort=False)

coarse_grained_df.to_csv('../../data/literature/Balakrishnan2021_proteome_fractions.csv',)