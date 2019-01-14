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
    Level1Design, EstimateModel, EstimateContrast)
from nipype.interfaces.fsl import Overlay, Slicer
from processing_functions import (
    conf_choose, con_setup, pick_onsets_spespk)
from nipype.interfaces.ants import ApplyTransforms
import nipype.interfaces.utility as niu

# Start timing
start = time.time()

# Specify location of data.
data_path = os.path.abspath('/data/pt_nmc002/other/narps/derivatives/')
event_path = os.path.abspath('/data/pt_nmc002/other/narps/event_tsvs/')
data_dir = os.path.join(data_path, 'fmriprep/')
out_dir = os.path.join(data_path, 'first_lev/')

subs = [
    "sub-001"
]

task_list = ['']
TR = 1.0  # TR in secs
fwhm = [4]

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
                           '{subject_id}_task-MGT_run-01_events.tsv')
}
selectderivs = Node(SelectFiles(templates,
                                sort_filelist=True),
                    name='selectderivs')

'''-----------------------------
|   The workflow starts here   |
-----------------------------'''

             firstlev = Workflow(name='firstlev',
                               base_dir=out_dir+'/tmp')
             firstlev.config['execution']['crashdump_dir'] = base_dir = out_dir +
             '/tmp/crash_files'
             firstlev.connect(infosource, 'subject_id',
                              selectderivs, 'subject_id')

             '''---------------------
|   First level setup  |
---------------------'''

             modelspec = Node(SpecifySPMModel(),
                            overwrite=False,
                            name='modelspec')
             modelspec.inputs.concatenate_runs = False
             modelspec.inputs.input_units = 'secs'
             modelspec.inputs.output_units = 'secs'
             modelspec.inputs.time_repetition = TR
             modelspec.inputs.high_pass_filter_cutoff = 128

             firstlev.connect([
                 getsubinforuns, modelspec, [('subject_info', 'subject_info')]),
    (selectderivs, modelspec, [('conf',
                                'realignment_parameters')]),
    (selectderivs, modelspec, [('func', 'functional_runs')])
])


'''-----------------
|   Run the jobs   |
-----------------'''
firstlev.write_graph(dotfilename = 'func_preproc_spespk.dot',
                     graph2use = 'orig', format = 'pdf')  # , simple_form=True
# plugin='MultiProc', plugin_args={'n_procs' : 8}
firstlev.run(plugin = 'MultiProc', plugin_args = {'n_procs': 40})

# Time again and spit out difference.
stop=time.time()
if (stop-start) < 3600:
    print('Elapsed time: %.2f mins.' % ((stop-start)/60))
else:
    print('Elapsed time: %.2f hours.' % ((stop-start)/60/60))
