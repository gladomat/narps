# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-


def get_niftis(datapath, subject_id=[]):
    """
    Loads BIDS formatted files.
    """

    from bids.layout import BIDSLayout
    # Get rid of "sub" from "sub-xx"
    subject_id = subject_id[-2:]
    layout = BIDSLayout(datapath)

    from bids.layout import BIDSLayout
    # Get rid of "sub" from "sub-xx"
    subject_id = subject_id[-2:]
    layout = BIDSLayout(datapath)

    func = layout.get(subject=subject_id,
                      modality='func',
                      type='bold',
                      task='SPINN',
                      return_type='file',
                      extensions=['nii', 'nii.gz'])
    # Whole brain EPI
    wb = layout.get(subject=subject_id,
                    modality='func',
                    type='bold',
                    task='wholehead',
                    return_type='file',
                    extensions=['nii', 'nii.gz'])[0]

    func_des = layout.get(subject=subject_id,
                          modality='func',
                          type='bold',
                          return_type='file',
                          extensions=['json'])

    UNI = layout.get(subject=subject_id,
                     modality='anat',
                     type='UNI',
                     return_type='file',
                     extensions=['nii', 'nii.gz'])

    events = layout.get(subject=subject_id,
                        modality='func',
                        type='events',
                        return_type='file',
                        extensions=['txt', 'tsv'])

    physio = layout.get(subject=subject_id,
                        modality='func',
                        type='physio',
                        return_type='file',
                        extensions=['txt', 'tsv'])

    fmaps = layout.get(subject=subject_id,
                       modality='fmap',
                       return_type='file',
                       extensions=['nii', 'nii.gz'])

    subject_id = 'sub-' + subject_id

    fmapMag = [f for f in fmaps if layout.get_metadata(f)[
        'ImageType'][2] == 'M']
    fmapPh = [f for f in fmaps if layout.get_metadata(
        f)['ImageType'][2] == 'P'][0]

    if not fmapMag == []:
        fmapM1, fmapM2 = fmapMag
    else:
        fmapM1, fmapM2 = [], []

    return (func, wb, events, physio, UNI, fmapM1, fmapPh, subject_id)


def merge_regress(realignment_parameters, physio):
    """
    Merge realignment_parameters and physiological regressors.
    """

    import tempfile
    import csv

    merged_files = []
    for files in zip(realignment_parameters, physio):
        dd = []
        for file in files:
            with open(file, 'rt') as handle:
                filerows = list(handle)
            data = []
            for row in filerows:
                data.append(row.split())
            dd.append(data)
        ddmerged = [d0+d1 for (d0, d1) in zip(dd[0], dd[1])]
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as handle:
            writer = csv.writer(handle, delimiter='\t')
            for row in ddmerged:
                writer.writerow(row)
            merged_files.append(handle.name)

    return merged_files


def choose_subs(group):
    """
    Choose subjects based on group characteristics:
    eq_range: equal range group
    eq_indiff: equal indifference group
    """

    import pandas as pd

    subs_df = pd.read_csv(
        '/data/pt_nmc002/other/narps/event_tsvs/participants.tsv',
        sep='\t', header=0)

    if group is 'eq_range':
        subs = subs_df.participant_id[
            subs_df['group'] == 'equalRange'
        ]
    else:
        subs = subs_df.participant_id[
            subs_df['group'] == 'equalIndifference'
        ]

    return list(subs)


def conf_choose(conf):
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
        conf_list.append(temp_list)
    return conf_list


def con_setup(subject_id):
    '''
    ======================================================================
    Contrast Setup
    ======================================================================
    This is a list of lists. The inner list specifies the contrasts and
    has the following format -
    [Name,Stat,[list of condition names],[weights on those conditions].
    The condition names must match the names listed in the subjectinfo
    function described above.
    '''

    runs = 4

    gamble = []
    gain = []
    loss = []
    for ir in range(runs):
        ir += 1
        gamble.append('gamble_run%i' % ir)
        gain.append('gamble_run%ixgain_run%i^1' % (ir, ir))
        loss.append('gamble_run%ixloss_run%i^1' % (ir, ir))

    pos_1 = [1] * runs
    neg_1 = [-1] * runs

    gamble = (
        'gamble', 'T',
        gamble, pos_1)

    gain = (
        'gain', 'T',
        gain, pos_1)

    loss = (
        'loss', 'T',
        loss, pos_1)

    contrasts = [
        gamble, gain, loss
    ]

    return contrasts


def get_mask_size(mask):
    # Generate output regarding size of ROI from fslstats, and get rid
    # of the last two numbers which represent t_min and t_size.

    import subprocess
    command = "fslstats " + mask + " -w"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    # get rid of last 2 numbers and extra chars
    output = output[:-6].decode("utf-8")
    return output


def pick_onsets(subject_id):
    '''Picks onsets and durations per condition and adds them to lists.
    This function specifically picks onsets for the speech vs speaker
    where the presentation is clear or in noise.
    The function accepts event files.

    'subject_id' is a string, i.e., sub-001
    '''

    cond_names = ['gamble']
    onset = {}
    duration = {}
    weights_gain = {}
    weights_loss = {}
    runs = ['01', '02', '03', '04']

    for r in range(len(runs)):  # Loop over number of runs.
        onset.update({s + '_run' + str(r+1): [] for s in cond_names})
        duration.update({s + '_run' + str(r+1): [] for s in cond_names})
        weights_gain.update({'gain_run' + str(r+1): []})
        weights_loss.update({'loss_run' + str(r+1): []})

    base_name = '/data/pt_nmc002/other/narps/event_tsvs/'
    # subject_id = 'sub-001'
    for ir, run in enumerate(runs):
        f_events = base_name + subject_id + \
            '_task-MGT_run-' + runs[ir] + '_events.tsv'
        with open(f_events, 'rt') as f:
            next(f)  # skip the header
            for line in f:
                info = line.strip().split()
                for cond in cond_names:
                    val = cond + '_run' + str(ir+1)
                    val_gain = 'gain_run' + str(ir+1)
                    val_loss = 'loss_run' + str(ir+1)
                    onset[val].append(float(info[0]))
                    duration[val].append(float(info[1]))
                    weights_gain[val_gain].append(float(info[2]))
                    weights_loss[val_loss].append(float(info[3]))
    #                if cond == 'gain':
    #                    weights[val].append(float(info[2]))
    #                elif cond == 'loss':
    #                    weights[val].append(float(info[3]))
    #                elif cond == 'task-activ':
    #                    weights[val].append(float(1))
    from nipype.interfaces.base import Bunch

    # Bunching is done per run, i.e. cond1_run1, cond2_run1, etc.
    subjectinfo = []
    for r in range(len(runs)):

        cond = [c + '_run' + str(r+1) for c in cond_names]
        gain = 'gain_run' + str(r+1)
        loss = 'loss_run' + str(r+1)

        subjectinfo.insert(r,
                           Bunch(conditions=cond,
                                 onsets=[onset[k] for k in cond],
                                 durations=[duration[k] for k in cond],
                                 amplitudes=None,
                                 tmod=None,
                                 pmod=[Bunch(name=[gain, loss],
                                             poly=[1, 1],
                                             param=[weights_gain[gain],
                                                    weights_loss[loss]])],
                                 regressor_names=None,
                                 regressors=None))

    return subjectinfo
