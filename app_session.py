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
model_dir = '/mnt/c/Users/shasa/Desktop/models/' #model store location for local testing

class infer:
    def __init__(self, df_lossy_200):
        # Settings
        self.ref_data = df_lossy_200
        self.experiments = self.ref_data.columns.tolist()
        self.experiments.remove('feature')
        self.stats_data = pd.read_csv(gene_stats_dir,index_col=0)
        self.gene_categories = pd.read_csv(gene_cats_dir)
        self.s3_client = boto3.client('s3')
        self.s3_bucket_name = 'ceres-infer'
        self.s3 = boto3.resource('s3', aws_access_key_id='AKIAWOWWNKWBQ3GACTFT',
                                  aws_secret_access_key='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',
                                  region_name='us-east-1')
        self.my_bucket = self.s3.Bucket(self.s3_bucket_name)
                 
    def infer_gene(self, gene2analyz, aws_s3):
        model_key = 'L200_models/model_rd10_' + gene2analyz + '.pkl'
        feats_key = 'L200_models/feats_' + gene2analyz + '.csv'
        
        if not aws_s3:
            model_file = join(model_dir, model_key)
            feats_file = join(model_dir, feats_key)
            model_pkl = pickle.load(open(model_file,'rb'))
            feats = pd.read_csv(feats_file, index_col=0)
            top_feats = feats[0:10]['feature'].tolist()
            ceres_top_feats = self.ref_data[self.ref_data.feature.isin(top_feats)]
            pred = model_pkl.predict(ceres_top_feats[self.experiments].to_numpy().T)
            return pred

        elif aws_s3:
            obj = self.s3.Object(self.s3_bucket_name, model_key)
            data = obj.get()['Body'].read()
            model_pkl = pickle.load(io.BytesIO(data))
            obj = self.s3.Object(self.s3_bucket_name, feats_key)
            data = obj.get()['Body'].read()
            feats = pd.read_csv(io.BytesIO(data), index_col=0)
            top_feats = feats[0:10]['feature'].tolist();
            ceres_top_feats = self.ref_data[self.ref_data.feature.isin(top_feats)]
            pred = model_pkl.predict(ceres_top_feats[self.experiments].to_numpy().T)
            return pred
        
    def calc_z_score(self, gene2analyz, pred_score):
        x_avg = self.stats_data.loc[gene2analyz][0]
        x_std = self.stats_data.loc[gene2analyz][1]
        z_score = (pred_score - x_avg) / x_std
        return x_avg, x_std, z_score
    
    def infer_gene_list(self, genes2analyz, aws_s3):
        genes = []
        ceres_pred = np.zeros(shape=(len(genes2analyz), len(self.experiments)))
        z_scores = np.zeros(shape=(len(genes2analyz), len(self.experiments)))
        x_avgs = []
        x_stds = []
        for i, gene in enumerate(genes2analyz):
            logging.info(gene)
            pred = self.infer_gene(gene, aws_s3=aws_s3)
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