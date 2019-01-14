import os
import glob


def conf_choose(conf, output_dir):
    """
    Choose specific confounds from csv file. Also computes the
    first derivative of the motion parameters
    """

    import numpy as np
    import pandas as pd

    conf_list = []
    for iconf in conf:
        df = pd.read_csv(iconf, sep='\t', header=0)
        temp_list = np.array([df['CSF'], df['WhiteMatter'],
                              df['X'], df['Y'], df['Z'],
                              df['RotX'], df['RotY'], df['RotZ'],
                              np.insert(np.diff(df['X']), 0, 0),
                              np.insert(np.diff(df['Y']), 0, 0),
                              np.insert(np.diff(df['Z']), 0, 0),
                              np.insert(np.diff(df['RotX']), 0, 0),
                              np.insert(np.diff(df['RotY']), 0, 0),
                              np.insert(np.diff(df['RotZ']), 0, 0)])
        df_out = pd.DataFrame(np.transpose(temp_list))
        df_out.to_csv(output_dir+os.path.basename(iconf),
                      header=None, index=None, sep='\t')


data_path = os.path.abspath(
    '/data/pt_nmc002/other/narps/derivatives/fmriprep/')
out_dir = '/data/pt_nmc002/other/narps/derivatives/confounds/'

dir_list = next(os.walk(data_path, '.'))[1]
subs = [x for x in dir_list if 'sub' in x]

for sub in subs:
    subdir = os.path.join(data_path, sub, 'func/')
    conf = glob.glob(subdir+'*confounds.tsv')
    conf_choose(conf, out_dir)
