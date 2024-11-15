#!/usr/local/bin/python3.9
import argparse
import pandas as pd
import numpy as np
import pybedtools
import sys

from scipy import stats

import warnings
warnings.filterwarnings("ignore")

dict_for_metrics = {
    '0.0 / 0.0' : 'TN',
    '1.0 / 1.0' : 'TP',
    '2.0 / 2.0' : 'TP',
    '-1.0 / -1.0' : 'TP',
    '-2.0 / -2.0' : 'TP',
    '2.0 / 1.0' : 'TP',
    '1.0 / 2.0' : 'TP',
    '-2.0 / -1.0' : 'TP',
    '-1.0 / -2.0' : 'TP',
    '- / 0.0' : 'TN',
    '- / 1.0': 'FN',
    '- / -1.0': 'FN',
    '- / -2.0': 'FN',
    '- / 2.0': 'FN',
    '0.0 / 1.0': 'FN',
    '0.0 / -1.0': 'FN',
    '0.0 / -2.0': 'FN',
    '0.0 / 2.0': 'FN',
    '1.0 / 0.0': 'FP',
    '-1.0 / 0.0': 'FP',
    '-2.0 / 0.0': 'FP',
    '2.0 / 0.0': 'FP',
    '-1.0 / 1.0': 'FP',
    '-1.0 / 2.0': 'FP',
    '1.0 / -1.0': 'FP',
    '1.0 / -2.0': 'FP',
    '-2.0 / 2.0': 'FP',
    '-2.0 / 1.0': 'FP',
    '2.0 / -2.0': 'FP',
    '2.0 / -1.0': 'FP'
}

def exit_with_code(msg: str, error_code: int) -> None:
    print(f"{msg} Failed with code {error_code}")
    sys.exit(error_code)


def intersect_test_ref(df: pd.DataFrame, ref: pd.DataFrame, cell_line: str) -> pd.DataFrame:
    df['norm_total'] = df['norm_total'].astype(float)
    ref_c = ref[['Hugo_symbol', cell_line]].rename(columns={cell_line:'norm_total_ref'})
    df = df.merge(ref_c, on=['Hugo_symbol'], how='outer')

    return df


def extract_counts_cna(df_ref_for_metrics: pd.DataFrame):
    dict_count = dict(df_ref_for_metrics['Status'].value_counts())
    try:
        TP = dict_count['TP']
    except:
        TP = 0
    try:
        FP = dict_count['FP']
    except:
        FP = 0
    try:
        TN = dict_count['TN']
    except:
        TN = 0
    try:
        FN = dict_count['FN']
    except:
        FN = 0
    return TP, FP, FN, TN


def sensitivity_specificity(df: pd.DataFrame, ref: pd.DataFrame, genes: list, cell_line: str):
    df_ref = intersect_test_ref(df, ref, cell_line)
    df_ref['norm_total'] = df_ref['norm_total'].fillna('-')
    df_ref['norm_total_ref'] = df_ref['norm_total_ref'].fillna('-')
    df_ref['concordance'] = df_ref['norm_total'].astype(str) + ' / ' + df_ref['norm_total_ref'].astype(str)
    df_ref_for_metrics = df_ref[(df_ref['Hugo_symbol'].isin(genes)) & (df_ref['norm_total_ref']!= '-')]
    df_ref_for_metrics['Status'] = df_ref_for_metrics['concordance'].map(dict_for_metrics)
    tp, fp, fn, tn = extract_counts_cna(df_ref_for_metrics)

    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)

    return sensitivity, specificity


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input',
                        help='.tsv file with normalized variants',
                        required=True)
    parser.add_argument('--reference', help='Reference file for cnv and neutral genes',
                         required=True)
    parser.add_argument('--genes', help='Target genes for testing',
                       required=True)
    parser.add_argument('--cell_line', help='Cell line for testing',
                        choices=[
                            'COLO829', 'HCC1143',
                            'HCC1937', 'HCC1395',
                            'NCI-H1770'
                        ],
                       required=True)
    parser.add_argument('--output', help='Output .txt file with performance characteristics',
                        required=True)
    return parser.parse_args()

# arguments
args = parse_args()
df_test = pd.read_csv(args.input, sep='\t', comment='#')
reference = pd.read_csv(args.reference, sep='\t', comment='#')
genes = list(pd.read_csv(args.genes, sep='\t', header=None)[0])

# assessment
sensitivity, specificity = sensitivity_specificity(df_test, reference, genes, args.cell_line)

# output
with open(args.output, 'w') as output:
    output.write(f'Sensitivity: {sensitivity * 100}' + '\n')
    output.write(f'Specificity: {specificity * 100}')
