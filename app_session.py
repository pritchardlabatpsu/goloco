# sessions
import os
from os.path import join
import numpy as np
import pandas as pd
import pickle
import logging
import boto3
import sys
import csv
import io
from tqdm import tqdm

#define directories for project:
gene_stats_dir = './data/19q4_sum_stats.csv'
gene_cats_dir = './data/19q4_gene_cats.csv'
model_dir = './ceres-infer/' #model store location for local testing

class infer:
    def __init__(self, df_lossy_200):
        # Settings
        self.ref_data = df_lossy_200
        self.experiments = self.ref_data.columns.tolist()
        self.experiments.remove('feature')
        self.stats_data = pd.read_csv(gene_stats_dir,index_col=0)
        self.gene_categories = pd.read_csv(gene_cats_dir)
                 
    def infer_gene(self, gene2analyz):
        model_key = 'L200_models/model_rd10_' + gene2analyz + '.pkl'
        feats_key = 'L200_models/feats_' + gene2analyz + '.csv'
        
        model_file = join(model_dir, model_key)
        feats_file = join(model_dir, feats_key)
        model_pkl = pickle.load(open(model_file,'rb'))
        feats = pd.read_csv(feats_file, index_col=0)
        top_feats = feats[0:10]['feature'].tolist()
        ceres_top_feats = self.ref_data[self.ref_data.feature.isin(top_feats)]
        pred = model_pkl.predict(ceres_top_feats[self.experiments].to_numpy().T)
        return pred
        
    def calc_z_score(self, gene2analyz, pred_score):
        x_avg = self.stats_data.loc[gene2analyz][0]
        x_std = self.stats_data.loc[gene2analyz][1]
        z_score = (pred_score - x_avg) / x_std
        return x_avg, x_std, z_score
    
    def infer_gene_list(self, genes2analyz):
        genes = []
        ceres_pred = np.zeros(shape=(len(genes2analyz), len(self.experiments)))
        z_scores = np.zeros(shape=(len(genes2analyz), len(self.experiments)))
        x_avgs = []
        x_stds = []
        for i, gene in enumerate(genes2analyz):
            logging.info(gene)
            pred = self.infer_gene(gene)
            x_avg, x_std, z_score = self.calc_z_score(gene, pred)
            genes.append(gene)
            ceres_pred[i] = pred
            z_scores[i] = z_score
            x_avgs.append(x_avg)
            x_stds.append(x_std)
        df_pred = pd.DataFrame()
        df_pred['gene'] = genes
        df_pred['avg'] = x_avgs
        df_pred['std'] = x_stds
        for i, exp in enumerate(self.experiments):
            df_pred[exp] = ceres_pred[:, i]
            df_pred['z_score'+exp] = z_scores[:, i]
        df_pred = pd.merge(df_pred, self.gene_categories, on='gene', how='left')
        df_pred['gene_category'] = df_pred['gene_category'].replace(np.nan, 'conditional essential')
        df_pred.set_index('gene')
        return df_pred