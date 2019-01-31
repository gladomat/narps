# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
import time
from os.path import join as opj
import os
# the workflow and node wrappers
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.utility import Function, IdentityInterface
from nipype.algorithms.modelgen import SpecifySPMModel
from nipype.interfaces.spm import (
    Level1Design, EstimateModel, EstimateContrast, Smooth)
from nipype.interfaces.fsl import Overlay, Slicer
from processing_functions import (
    conf_choose, con_setup, pick_onsets, choose_subs)
from nipype.interfaces.ants import ApplyTransforms
from nipype.algorithms.misc import Gunzip
import nipype.interfaces.utility as niu

# Start timing
start = time.time()

# Specify location of data.
data_path = os.path.abspath('/data/pt_nmc002/other/narps/derivatives/')
event_path = os.path.abspath('/data/pt_nmc002/other/narps/event_tsvs/')
data_dir = os.path.join(data_path, 'fmriprep/')
out_dir = os.path.join(data_path, 'first_lev/')
# Group is chosen to be either equal range or equal indifference
# group = 'eq_indiff'  # eq_indiff
group = 'eq_range'

#subs = choose_subs(group)

subs = [
    "sub-002"
]

task_list = ['']
TR = 1.0  # TR in secs
fwhmlist = [5]
template_image = os.path.join(data_path,
                              'fmriprep/sub-001/anat/sub-001_T1w_space-MNI152NLin2009cAsym_preproc.nii.gz')

'''--------------
|   Grab data   |
--------------'''

# Need an identitiy interface to iterate over subject_id and run
infosource = Node(interface=IdentityInterface(fields=['subject_id']),
                  name="infosource")
infosource.iterables = [('subject_id', subs)]

# Select files from derivatives.
templates = {'func': opj(data_dir,
                         '{subject_id}', 'func',
                         '{subject_id}_task-MGT_run-*_bold_space-MNI152NLin2009cAsym_preproc.nii.gz'),
             'conf': opj(data_path, 'confounds',
                         '{subject_id}_task-MGT_run-*_bold_confounds.tsv'),
             'events': opj(event_path,
                           '{subject_id}_task-MGT_run-*_events.tsv')
             }
selectderivs = Node(SelectFiles(templates,
                                sort_filelist=True),
                    name='selectderivs')

'''-----------------------------
|   The workflow starts here   |
-----------------------------'''

firstlev = Workflow(name='firstlev',
                    base_dir=out_dir+'/tmp')
firstlev.config['execution']['crashdump_dir'] = base_dir = out_dir + \
    '/tmp/crash_files'
firstlev.connect(infosource, 'subject_id',
                 selectderivs, 'subject_id')

'''---------------------
|   First level setup  |
---------------------'''

gunzip = MapNode(Gunzip(), name='gunzip', iterfield=['in_file'])
firstlev.connect(selectderivs, 'func', gunzip, 'in_file')

# Smooth warped functionals. Watch out it smoothes again if you stop here!
smooth = Node(Smooth(),
              overwrite=False,
              name="smooth")
smooth.iterables = ("fwhm", fwhmlist)
firstlev.connect(gunzip, 'out_file', smooth, 'in_files')

getsubinforuns = Node(Function(input_names=["subject_id"],
                               output_names=["subject_info"],
                               function=pick_onsets),
                      name='getsubinforuns')

modelspec = Node(SpecifySPMModel(),
                 overwrite=False,
                 name='modelspec')
modelspec.inputs.concatenate_runs = False
modelspec.inputs.input_units = 'secs'
modelspec.inputs.output_units = 'secs'
modelspec.inputs.time_repetition = TR
modelspec.inputs.high_pass_filter_cutoff = 128

firstlev.connect([
    (infosource, getsubinforuns, [('subject_id', 'subject_id')]),
    (getsubinforuns, modelspec, [('subject_info', 'subject_info')]),
    (selectderivs, modelspec, [('conf',
                                'realignment_parameters')]),
    (smooth, modelspec, [('smoothed_files', 'functional_runs')])
])

level1design = Node(Level1Design(),
                    overwrite=False,
                    name='level1design')
level1design.inputs.bases = {'hrf': {'derivs': [0, 0]}}
level1design.inputs.timing_units = 'secs'
level1design.inputs.interscan_interval = TR
firstlev.connect([
    (modelspec, level1design, [('session_info', 'session_info')])
])

level1estimate = Node(EstimateModel(),
                      overwrite=False,
                      name='level1estimate')
level1estimate.inputs.estimation_method = {'Classical': 1}
firstlev.connect(level1design, 'spm_mat_file',
                 level1estimate, 'spm_mat_file')

contrast_estimate = Node(EstimateContrast(),
                         overwrite=False,
                         name='contraste_estimate')
contrast_estimate.config = {'execution': {'remove_unnecessary_outputs': False}}
firstlev.connect([
    (level1estimate, contrast_estimate,
     [('spm_mat_file', 'spm_mat_file'),
      ('beta_images', 'beta_images'),
      ('residual_image', 'residual_image')])
])

contrasts = Node(Function(function=con_setup,
                          input_names=['subject_id'],
                          output_names=['contrasts']),
                 name='contrasts')
firstlev.connect([
    (infosource, contrasts, [('subject_id', 'subject_id')]),
    (contrasts, contrast_estimate, [('contrasts', 'contrasts')])
])


'''--------------------------------------
|   Create overlays for quick checking  |
--------------------------------------'''
# Use :class: nipype.interfaces.utility.Select to select each contrast
# for reporting.
selectcontrast = Node(niu.Select(), name="selectcontrast")
# Iterate over each contrast and create report images.
selectcontrast.iterables = ('index', [[x] for x in range(3)])

# Use nipype.interfaces.fsl.Overlay to combine the statistical output of the
# contrast estimate and a background image into one volume.
overlaystats = Node(Overlay(),
                    # iterfield=['background_image'],
                    name="overlaystats")
overlaystats.inputs.stat_thresh = (1.5, 10)
overlaystats.inputs.show_negative_stats = True
overlaystats.inputs.auto_thresh_bg = True
overlaystats.inputs.background_image = template_image

# Use nipype.interfaces.fsl.Slicer to create images of the overlaid
# statistical volumes for a report of the first-level results.
slicestats = Node(Slicer(),
                  # iterfield='in_file',
                  name="slicestats")
slicestats.inputs.sample_axial = 4
slicestats.inputs.image_width = 1500

firstlev.connect([
    (contrast_estimate, selectcontrast, [('con_images', 'inlist')]),
    (selectcontrast, overlaystats, [('out', 'stat_image')]),
    (overlaystats, slicestats, [('out_file', 'in_file')]),
])

'''-----------------
|   Run the jobs   |
-----------------'''
firstlev.write_graph(dotfilename='func_preproc_spespk.dot',
                     graph2use='orig', format='pdf')  # , simple_form=True
# plugin='MultiProc', plugin_args={'n_procs' : 8}
firstlev.run(plugin='MultiProc', plugin_args={'n_procs': 32})

# Time again and spit out difference.
stop = time.time()
if (stop-start) < 3600:
    print('Elapsed time: %.2f mins.' % ((stop-start)/60))
else:
    print('Elapsed time: %.2f hours.' % ((stop-start)/60/60))
