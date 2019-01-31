# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
import time
from os.path import join as opj
import os
from nipype.interfaces.io import SelectFiles, DataSink, DataGrabber
from nipype.interfaces.spm import (TwoSampleTTestDesign, EstimateModel,
                                   EstimateContrast, Threshold)
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.pipeline.engine import Workflow, Node
import nipype.interfaces.utility as niu
from nipype.interfaces.fsl import Overlay, Slicer
from processing_functions import choose_subs

# Start timing
start = time.time()

# Specify location of data.
data_dir = os.path.abspath(
    '/data/pt_nmc002/other/narps/derivatives/')
out_dir = os.path.abspath(
    '/data/pt_nmc002/other/narps/derivatives/second_lev_group_compar/')
deriv_dir = os.path.join(data_dir, 'first_lev/equalRange')
fwhm = [5]
mask = '/data/pt_nmc002/other/narps/derivatives/fmriprep/gr_mask_tmax.nii'
#template_image = '/data/pt_neunmc024/SpeechInNoise_Dyslexia_Study/MRI_Data_BIDS/derivatives/mni_template/mni1mm.nii'

# list of contrast identifiers
contrasts = ['con_0003']

# collect all the con images for each contrast.
contrast_ids = list(range(1, len(contrasts) + 1))

# Infosource - a function free node to iterate over the list of subject names
infosource = Node(IdentityInterface(fields=['contrast_id']),
                  name="infosource")
infosource.iterables = [('contrast_id', contrasts)]

# Select files from derivatives.

templates = {'cons_eqind': '/data/pt_nmc002/other/narps/derivatives/first_lev/tmp/equalIndifference/*/contraste_estimate/{contrast_id}.nii',
             'cons_eqrng': '/data/pt_nmc002/other/narps/derivatives/first_lev/tmp/equalRange/*/contraste_estimate/{contrast_id}.nii'}
selectderivs = Node(SelectFiles(templates,
                                sort_filelist=True),
                    name='selectderivs')

# Two Sample T-Test Design
twosampttest = Node(TwoSampleTTestDesign(),
                    # overwrite=True,
                    name="twosampttest")
twosampttest.inputs.explicit_mask_file = mask

# EstimateModel - estimate the parameters of the model
# Even for second level it should be 'Classical': 1.
level2estimate = Node(EstimateModel(estimation_method={'Classical': 1}),
                      name="level2estimate")
# EstimateContrast - estimates simple group contrast
level2conestimate = Node(EstimateContrast(group_contrast=True),
                         name="level2conestimate")
cont1 = ['Eq range vs Eq indiff in loss', 'T', ['mean'], [1, -1]]
level2conestimate.inputs.contrasts = [cont1]

# Create the 2nd level pipeline
secondlev = Workflow(name='secondlev_group_compar', base_dir=out_dir+'/tmp')
secondlev.config['execution']['crashdump_dir'] = base_dir = out_dir + \
    '/tmp/crash_files'

secondlev.connect([
    (infosource, selectderivs, [('contrast_id', 'contrast_id')]),
    (selectderivs, twosampttest, [('cons_eqrng', 'group1_files'),
                                  ('cons_eqind', 'group2_files')]),
    (twosampttest, level2estimate, [('spm_mat_file', 'spm_mat_file')]),
    (level2estimate, level2conestimate, [('spm_mat_file', 'spm_mat_file'),
                                         ('beta_images', 'beta_images'),
                                         ('residual_image', 'residual_image')])
])


'''-----------------
|   Run the jobs   |
-----------------'''
secondlev.run(plugin='MultiProc', plugin_args={'n_procs': 16})

# Time again and spit out difference.
stop = time.time()
if (stop-start) < 3600:
    print('Elapsed time: %.2f mins.' % ((stop-start)/60))
else:
    print('Elapsed time: %.2f hours.' % ((stop-start)/60/60))
