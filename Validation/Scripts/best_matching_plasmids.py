#Hungarian for max coverage_percentage and q_perc_cov only
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

file_path = r'Documents\plsMD new val\gplas2_rec_prec_sep.txt'
qseqid_subset_file = r'Documents\plsMD new val\gplas_qseqid.txt'
sseqid_subset_file = r'Documents\plsMD new val\sseqid_new_oct.txt'
output_file_all_samples = r'Documents\plsMD new val\October\gplas_hung_both.txt'

df = pd.read_csv(file_path, delim_whitespace=True)

df_qseqid_subset = pd.read_csv(qseqid_subset_file, delim_whitespace=True, encoding='ISO-8859-1')
df_sseqid_subset = pd.read_csv(sseqid_subset_file, delim_whitespace=True, encoding='ISO-8859-1')

df['coverage_percentage'] = df['coverage_percentage'].str.rstrip('%').astype(float)

unique_samples = df_sseqid_subset['Sample'].unique()

all_best_matches = []  
all_coverage_percentages = []  
all_q_prec_covs = [] 

coverage_weight = 0.5 
q_prec_cov_weight = 0.5  

for sample in unique_samples:
    sample_data = df[df['Sample'] == sample]

    unique_qseqids = df_qseqid_subset[df_qseqid_subset['Sample'] == sample]['qseqid'].unique()
    unique_sseqids = df_sseqid_subset[df_sseqid_subset['Sample'] == sample]['sseqid'].unique()

    num_qseqids = len(unique_qseqids)
    num_sseqids = len(unique_sseqids)

    if num_qseqids > num_sseqids:
        weight_matrix = np.zeros((num_qseqids, num_sseqids + 1))
    elif num_sseqids > num_qseqids:
        weight_matrix = np.zeros((num_qseqids + 1, num_sseqids))
    else:
        weight_matrix = np.zeros((num_qseqids, num_sseqids))

    for i, qseqid in enumerate(unique_qseqids):
        for j, sseqid in enumerate(unique_sseqids):
            match = sample_data[(sample_data['qseqid'] == qseqid) & (sample_data['sseqid'] == sseqid)]
            if not match.empty:
                coverage = match['coverage_percentage'].values[0]
                q_prec_cov = match['q_prec_cov'].values[0]
                weight_matrix[i, j] = coverage_weight * coverage + q_prec_cov_weight * q_prec_cov
            else:
                weight_matrix[i, j] = 0

    row_indices, col_indices = linear_sum_assignment(-weight_matrix)

    matches = []
    for row, col in zip(row_indices, col_indices):
        if row < len(unique_qseqids) and col < len(unique_sseqids): 
            matches.append((unique_qseqids[row], unique_sseqids[col]))

    best_matches_df = sample_data[sample_data.apply(lambda x: (x['qseqid'], x['sseqid']) in matches, axis=1)]


    if not best_matches_df.empty:
        matched_qseqids = best_matches_df['qseqid'].dropna().unique()
        matched_sseqids = best_matches_df['sseqid'].dropna().unique()
    else:
        matched_qseqids = []
        matched_sseqids = []

    unmatched_qseqids = [qseqid for qseqid in unique_qseqids if qseqid not in matched_qseqids]
    unmatched_sseqids = [sseqid for sseqid in unique_sseqids if sseqid not in matched_sseqids]

    unmatched_qseqid_rows = pd.DataFrame({
        'Sample': [sample] * len(unmatched_qseqids),
        'qseqid': unmatched_qseqids,
        'sseqid': [''] * len(unmatched_qseqids),
        'coverage_percentage': [np.nan] * len(unmatched_qseqids),
        'q_prec_cov': [0] * len(unmatched_qseqids),
    })
    unmatched_sseqid_rows = pd.DataFrame({
        'Sample': [sample] * len(unmatched_sseqids),
        'qseqid': [''] * len(unmatched_sseqids),
        'sseqid': unmatched_sseqids,
        'coverage_percentage': [0] * len(unmatched_sseqids),
        'q_prec_cov': [np.nan] * len(unmatched_sseqids),
    })

    best_matches_df = sample_data[sample_data.apply(lambda x: (x['qseqid'], x['sseqid']) in matches, axis=1)]
    best_matches_df = pd.concat([best_matches_df, unmatched_qseqid_rows, unmatched_sseqid_rows])

    all_best_matches.append(best_matches_df)
    all_coverage_percentages.extend(best_matches_df['coverage_percentage'].values)
    all_q_prec_covs.extend(best_matches_df['q_prec_cov'].values)

all_best_matches_df = pd.concat(all_best_matches)

all_best_matches_df.to_csv(output_file_all_samples, sep='\t', index=False)

average_recall = np.nanmean(all_coverage_percentages)
average_precision = np.nanmean(all_q_prec_covs)

print(f"\nAverage Recall (coverage_percentage): {average_recall:.2f}%")
print(f"Average Precision (q_prec_cov): {average_precision:.2f}%")
print(f"Results saved to: {output_file_all_samples}")