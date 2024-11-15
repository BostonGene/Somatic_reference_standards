#!/usr/local/bin/python3.9
import argparse
import pandas as pd
import numpy as np
import pybedtools
import sys

from scipy import stats

import warnings
warnings.filterwarnings("ignore")

COLUMNS = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']

def exit_with_code(msg: str, error_code: int) -> None:
    print(f"{msg} Failed with code {error_code}")
    sys.exit(error_code)


def intersect_maf_with_target(maf_df: pd.DataFrame, target_bed: str,
                              column_name: str) -> pd.DataFrame:
    """
    Intersects input_maf with target_bed using bedtools intersect
    :param maf_df: dataframe from input .maf file
    :param target_bed: .bed file with target (kit) regions
    :return: dataframe with new column with markup
    """

    target_bed_intervals = pybedtools.BedTool(target_bed)
    maf_intervals = pybedtools.BedTool(maf_df[
        ['Chromosome', 'Start_Position', 'End_Position']
    ].values.tolist())
    try:
        maf_target_intervals = maf_intervals.intersect(target_bed_intervals,
                                                       wa=True,
                                                       u=True).to_dataframe() \
            .rename(columns={
                'chrom': 'Chromosome',
                'start': 'Start_Position',
                'end': 'End_Position'
            })
        maf_target_intervals[column_name] = 'PASS'
        maf_df = maf_df.merge(maf_target_intervals, on=[
            'Chromosome',
            'Start_Position',
            'End_Position'
        ], how='left')
    except KeyError:
        print('There are no variants inside target region')
        maf_df[column_name] = ''
        return maf_df
    except Exception as e:
        exit_with_code(e.__str__() + '. intersect_maf_with_target()', 104)
    return maf_df


def make_ref_sens(ref_sens: pd.DataFrame, target: str, purity: int) -> pd.DataFrame:
    ref_sens = intersect_maf_with_target(ref_sens, target, '_target')
    ref_sens = ref_sens[ref_sens['_target']=='PASS'].drop(columns=['_target'])
    for i in [100, 75, 50, 30, 20, 10]:
        ref_sens['Tumor_VAF_'+str(i)] = ref_sens['Tumor_VAF_median'] * i/100
        print(ref_sens[ref_sens['Tumor_VAF_'+str(i)]>0.05].shape[0])
    ref_sens = ref_sens[ref_sens['Tumor_VAF_'+str(purity)]>0.05]
    return ref_sens


def df_short_prep(df: pd.DataFrame) -> pd.DataFrame:
    df_short = df[COLUMNS]
    return df_short

def tp_fn_assessment(df: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    df_short = df_short_prep(df)
    ref_short = df_short_prep(ref)
    ref_short['IN_REF'] = 'PASS'
    df_short['IN_PIPELINE'] = 'PASS'
    df_tot = df_short.merge(ref_short, on=COLUMNS, how='outer')
    df_tot['TP'] = np.where((df_tot['IN_REF'] == 'PASS') & (df_tot['IN_PIPELINE'] == 'PASS'), 'PASS', '')
    df_tot['FN'] = np.where((df_tot['IN_REF'] == 'PASS') & (df_tot['IN_PIPELINE'] != 'PASS'), 'PASS', '')
    return df_tot

def fp_assessment(df: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    ref_short = df_short_prep(ref)
    ref_short['NON_FP'] = 'PASS'
    df_tot = df.merge(ref_short, on=COLUMNS, how='left')
    df_tot['FP'] = np.where((df_tot['IN_PIPELINE'] == 'PASS') & (df_tot['NON_FP'] != 'PASS'), 'PASS', '')
    return df_tot


def sensitivity_precision(df: pd.DataFrame, ref_s: pd.DataFrame, ref_p: pd.DataFrame):
    df = tp_fn_assessment(df, ref_s)
    df = fp_assessment(df, ref_p)
    try:
        TP = df['TP'].value_counts()['PASS']
    except:
        TP = 0
    try:
        FN = df['FN'].value_counts()['PASS']
    except:
        FN = 0 
    try:
        FP = df['FP'].value_counts()['PASS']
    except:
        FP = 0
    sens = TP / (TP + FN)
    prec = TP / (TP + FP)
    return sens, prec

def prepare_df(df: pd.DataFrame, target: str, type: str) -> pd.DataFrame:
    df = intersect_maf_with_target(df, target, '_target')
    if type == 'SNP':
        f_list = ['SNP']
    else:
        f_list = ['INS', 'DEL']
    df = df[(df['_target']=='PASS') & (df['Variant_Type'].isin(f_list))]
    return df

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input',
                        help='.maf file with annotated variants',
                        required=True)
    parser.add_argument('--target_bed', help='Input target.bed file', required=True)
    parser.add_argument('--reference_sensitivity', help='Reference file for sensitivity',
                         required=True)
    parser.add_argument('--reference_precision', help='Reference file for precision',
                         required=True)
    parser.add_argument('--purity',
                        help='purity of the sample for analysis, \
                            must be eiter 10, 20, 30, 50, 75, 100 %',
                        required=True,
                        type=int)
    parser.add_argument('--type',
                        help='type of the events for analysis',
                        choices=['SNP', 'INDEL'],
                        required=True,
                        type=str.upper)
    parser.add_argument('--output', help='Output .txt file with performance characteristics',
                        required=True)
    return parser.parse_args()

# arguments
args = parse_args()
df_test_full = pd.read_csv(args.input, sep='\t', comment='#')
df_test = prepare_df(df_test_full, args.target_bed, args.type)
ref_sensitivity_full = pd.read_csv(args.reference_sensitivity, sep='\t', comment='#')
ref_sensitivity = make_ref_sens(ref_sensitivity_full, args.target_bed, args.purity)
ref_precision = pd.read_csv(args.reference_precision, sep='\t', comment='#')

# assessment
sensitivity, precision = sensitivity_precision(df_test, ref_sensitivity, ref_precision)

# output
with open(args.output, 'w') as output:
    output.write(args.type + '\n')
    output.write(f'Sensitivity: {sensitivity * 100}' + '\n')
    output.write(f'Precision: {precision * 100}')
