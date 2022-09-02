#!/usr/bin/env python
import pandas as pd
import os
from glob import glob


class LoadData:
    
    metdata_cons_cols = ["sequence", "charge", "RT", "mz", "quality", 'feature_ids', 'pred_formulas', 'pred_compounds']
    files_ext = "pkl"
    
    sep_fm_name_files = "fm"
    cons_file_name = "FeatureMatrix.tsv"
    feat_maps_folder_name = "feature_map_df"
    ms_df_folder_name = "ms_df"
    sep_sirius_formulas = "_formulas"
    sep_sirius_structure = "_structueçres"
    sep_sirius_sirius = "formulas"
    sirius_folder_name = "Sirius"

    ms2_spectra_col = "MS2_spectra"

    def __init__(self, results_path):
        self.results_path = results_path


    def loadConsensus(self):
        cons_file_path = os.path.join(self.results_path, \
                                    LoadData.cons_file_name)
        ## load consensus data
        consensus_df = pd.read_csv(cons_file_path, sep="\t").set_index(keys='id')

        ## rename columns
        consensus_df.columns = [col.split(".mzML")[0] if ".mzML" in col else col \
                                        for col in consensus_df.columns]
        ## name samples in consenus
        self.intensity_cons_cols = sorted(list(set(consensus_df.columns) - \
                                    set(LoadData.metdata_cons_cols)))
        # number of samples
        self.n_samples_cons = len(self.intensity_cons_cols)
        self.consensus_df = consensus_df

        ## consensus intensity data frame
        self.intensity_df = consensus_df[self.intensity_cons_cols]

    def loadFeatMapsData(self):
        
        str_to_list = lambda x: [int(id_) for id_ in x.split(',') if id_.isdigit()]
        feat_maps_folder_path = os.path.join(self.results_path, LoadData.feat_maps_folder_name)
        feat_maps_files = glob(os.path.join(feat_maps_folder_path, f"*.{LoadData.files_ext}"))

        dct_feat_maps_df = {}
        for file in feat_maps_files:
            file_name = os.path.basename(file)
            sample_name = file_name.split(LoadData.sep_fm_name_files)[0]
            
            feat_maps_df  = pd.read_pickle(file).sort_values(by='mz')
            feat_maps_df[LoadData.ms2_spectra_col + "_array"] = feat_maps_df[LoadData.ms2_spectra_col].\
                                                                    apply(lambda x: str_to_list(x))
            dct_feat_maps_df[sample_name] = feat_maps_df

        # number of samples
        self.n_samples_feat_map = len(dct_feat_maps_df)
        self.dct_feat_maps_df = dct_feat_maps_df


    # def str_list_to_list(str_array):
    #     str_array = str_array.replace("[","").replace("]","").strip()
    #     str_array = re.sub("\s\s+" , " ", str_array).split(" ")
    #     return str_array
    

    def loadMSDF(self):

        ms_df_folder_path = os.path.join(self.results_path, LoadData.ms_df_folder_name)
        ms_df_files = glob(os.path.join(ms_df_folder_path, f"*.{LoadData.files_ext}"))

        rt_col, mzarray_col, intarray_col = "RT", "mzarray", "intarray"

        dct_ms_df = {}
        for file in ms_df_files:

            file_name = os.path.basename(file)
            sample_name = file_name.split("."+LoadData.files_ext)[0]

            ms_df = pd.read_pickle(file)
            dct_au = {}
            for i in ms_df.index:
                ## Convert string arrays to arrays
                rt = ms_df.loc[i, rt_col]
                mzarray = (ms_df.loc[i, mzarray_col])
                intarray = (ms_df.loc[i, intarray_col])
                ## Generate DataFrame with mz and intensity
                mz_int_df = pd.DataFrame.from_dict({"mz": mzarray, "intensity": intarray}).sort_values(by='mz')
                dct_au[i] = {'RT':rt, 'mz_int_data': mz_int_df}
            dct_ms_df[sample_name] = dct_au
        
        # number of samples
        self.n_samples_ms_df = len(dct_ms_df)
        self.dct_ms_df = dct_ms_df
