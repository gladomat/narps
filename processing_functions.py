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
    if subject_id in ['sub-01', 'sub-02']:
        runs = 3
    else:
        runs = 4

    spe_cl = []
    spe_no = []
    spk_cl = []
    spk_no = []
    instr = []
    for ir in range(runs):
        ir += 1
        spe_cl.append('speech_clear_run%i' % ir)
        spe_no.append('speech_noise_run%i' % ir)
        spk_cl.append('speaker_clear_run%i' % ir)
        spk_no.append('speaker_noise_run%i' % ir)
        instr.append('instruction_run%i' % ir)

    pos_1 = [1] * runs
    neg_1 = [-1] * runs

    spe_clear = (
        'Speech clear', 'T',
        spe_cl, pos_1)

    spe_noise = (
        'Speech noise', 'T',
        spe_no, pos_1)

    spk_clear = (
        'Speaker clear', 'T',
        spk_cl, pos_1)

    spk_noise = (
        'Speaker noise', 'T',
        spk_no, pos_1)

    spe_list = [spe_cl, spe_no]
    spe_1s = [pos_1, pos_1]
    spe = (
        'Speech', 'T',
        [y for x in spe_list for y in x],
        [y for x in spe_1s for y in x])

    spk_list = [spk_cl, spk_no]
    spk_1s = [pos_1, pos_1]
    spk = (
        'Speaker', 'T',
        [y for x in spk_list for y in x],
        [y for x in spk_1s for y in x])

    clear_list = [spe_cl, spk_cl]
    clear_1s = [pos_1, pos_1]
    clear = (
        'Clear', 'T',
        [y for x in clear_list for y in x],
        [y for x in clear_1s for y in x])

    noise_list = [spe_no, spk_no]
    noise_1s = [pos_1, pos_1]
    noise = (
        'Noise', 'T',
        [y for x in noise_list for y in x],
        [y for x in noise_1s for y in x])

    spe_vs_spk_list = [spe_cl, spe_no, spk_cl, spk_no]
    spe_vs_spk_1s = [pos_1, pos_1, neg_1, neg_1]
    spe_vs_spk = (
        'Speech vs Speaker', 'T',
        [y for x in spe_vs_spk_list for y in x],
        [y for x in spe_vs_spk_1s for y in x],
    )

    spk_vs_spe_list = [spe_cl, spe_no, spk_cl, spk_no]
    spk_vs_spe_1s = [neg_1, neg_1, pos_1, pos_1]
    spk_vs_spe = (
        'Speaker vs Speech', 'T',
        [y for x in spk_vs_spe_list for y in x],
        [y for x in spk_vs_spe_1s for y in x]
    )

    clear_vs_noise_list = [spe_cl, spe_no, spk_cl, spk_no]
    clear_vs_noise_1s = [pos_1, neg_1, pos_1, neg_1]
    clear_vs_noise = (
        'Clear vs Noise', 'T',
        [y for x in clear_vs_noise_list for y in x],
        [y for x in clear_vs_noise_1s for y in x]
    )

    noise_vs_clear_list = [spe_cl, spe_no, spk_cl, spk_no]
    noise_vs_clear_1s = [neg_1, pos_1, neg_1, pos_1]
    noise_vs_clear = (
        'Noise vs Clear', 'T',
        [y for x in noise_vs_clear_list for y in x],
        [y for x in noise_vs_clear_1s for y in x]
    )

    interaction_1_list = [spe_cl, spe_no, spk_cl, spk_no]
    interaction_1_1s = [pos_1, neg_1, neg_1, pos_1]
    interaction_1 = (
        'Interaction spe vs spk (clear vs noise)', 'T',
        [y for x in interaction_1_list for y in x],
        [y for x in interaction_1_1s for y in x]
    )

    interaction_2_list = [spe_cl, spe_no, spk_cl, spk_no]
    interaction_2_1s = [neg_1, pos_1, pos_1, neg_1]
    interaction_2 = (
        'Interaction spe vs spk (noise vs clear)', 'T',
        [y for x in interaction_2_list for y in x],
        [y for x in interaction_2_1s for y in x]
    )

    spe_cl_vs_no_list = [spe_cl, spe_no]
    spe_cl_vs_no_1s = [pos_1, neg_1]
    spe_cl_vs_no = (
        'Speech: clear vs noise', 'T',
        [y for x in spe_cl_vs_no_list for y in x],
        [y for x in spe_cl_vs_no_1s for y in x]
    )

    spk_cl_vs_no_list = [spk_cl, spk_no]
    spk_cl_vs_no_1s = [pos_1, neg_1]
    spk_cl_vs_no = (
        'Speaker: clear vs noise', 'T',
        [y for x in spk_cl_vs_no_list for y in x],
        [y for x in spk_cl_vs_no_1s for y in x],
    )

    spe_no_vs_cl_list = [spe_no, spe_cl]
    spe_no_vs_cl_1s = [pos_1, neg_1]
    spe_no_vs_cl = (
        'Speech: noise vs clear', 'T',
        [y for x in spe_no_vs_cl_list for y in x],
        [y for x in spe_no_vs_cl_1s for y in x]
    )

    spk_no_vs_cl_list = [spk_no, spk_cl]
    spk_no_vs_cl_1s = [pos_1, neg_1]
    spk_no_vs_cl = (
        'Speaker: noise vs clear', 'T',
        [y for x in spk_no_vs_cl_list for y in x],
        [y for x in spk_no_vs_cl_1s for y in x],
    )

    instr_vs_task_list = [instr, spe_cl, spe_no, spk_cl, spk_no]
    instr_vs_task_1s = [pos_1, neg_1*4]
    instr_vs_task = (
        'Instruction vs Task', 'T',
        [y for x in instr_vs_task_list for y in x],
        [y for x in instr_vs_task_1s for y in x]
    )

    task_vs_instr_list = [instr, spe_cl, spe_no, spk_cl, spk_no]
    task_vs_instr_1s = [neg_1, pos_1*4]
    task_vs_instr = (
        'Task vs Instruction', 'T',
        [y for x in task_vs_instr_list for y in x],
        [y for x in task_vs_instr_1s for y in x]
    )

    contrasts = [
        spe_clear, spe_noise, spk_clear, spk_noise,
        spe, spk, clear, noise,
        spe_vs_spk, spk_vs_spe, clear_vs_noise, noise_vs_clear,
        interaction_1, interaction_2,
        spe_cl_vs_no, spe_no_vs_cl, spk_cl_vs_no, spk_no_vs_cl,
        instr_vs_task, task_vs_instr,
    ]

    return contrasts


def pick_onsets_spespk(events, subject_id):
    '''Picks onsets and durations per condition and adds them to lists.
    This function specifically picks onsets for the speech vs speaker
    where the presentation is clear or in noise.

    The function accepts event files.
    '''

    cond_names = ['speech_clear', 'speech_noise',
                  'speaker_clear', 'speaker_noise',
                  'instruction']
    onset = {}
    duration = {}

    if subject_id in ['sub-01', 'sub-02']:
        runs = ['01', '02', '03']
    else:
        runs = ['01', '02', '03', '04']

    for r in range(len(runs)):  # Loop over number of runs.
        onset.update({s + '_run' + str(r+1): [] for s in cond_names})
        duration.update({s + '_run' + str(r+1): [] for s in cond_names})

    print(runs)

    for ir, run in enumerate(runs):
        with open(events[ir], 'rt') as f:
            if run == '01':
                for line in f:
                    info = line.strip().split()
                    if info[2] == 'speech':
                        onset['speech_clear_run1'].append(float(info[0]))
                        duration['speech_clear_run1'].append(float(info[1]))
                    elif info[2] == 'speech_noise':
                        onset['speech_noise_run1'].append(float(info[0]))
                        duration['speech_noise_run1'].append(float(info[1]))
                    elif info[2] == 'speaker':
                        onset['speaker_clear_run1'].append(float(info[0]))
                        duration['speaker_clear_run1'].append(float(info[1]))
                    elif info[2] == 'speaker_noise':
                        onset['speaker_noise_run1'].append(float(info[0]))
                        duration['speaker_noise_run1'].append(float(info[1]))
                    elif info[2] == 'instruction':
                        onset['instruction_run1'].append(float(info[0]))
                        duration['instruction_run1'].append(float(info[1]))
            elif run == '02':
                for line in f:
                    info = line.strip().split()
                    if info[2] == 'speech':
                        onset['speech_clear_run2'].append(float(info[0]))
                        duration['speech_clear_run2'].append(float(info[1]))
                    elif info[2] == 'speech_noise':
                        onset['speech_noise_run2'].append(float(info[0]))
                        duration['speech_noise_run2'].append(float(info[1]))
                    elif info[2] == 'speaker':
                        onset['speaker_clear_run2'].append(float(info[0]))
                        duration['speaker_clear_run2'].append(float(info[1]))
                    elif info[2] == 'speaker_noise':
                        onset['speaker_noise_run2'].append(float(info[0]))
                        duration['speaker_noise_run2'].append(float(info[1]))
                    elif info[2] == 'instruction':
                        onset['instruction_run2'].append(float(info[0]))
                        duration['instruction_run2'].append(float(info[1]))
            elif run == '03':
                for line in f:
                    info = line.strip().split()
                    if info[2] == 'speech':
                        onset['speech_clear_run3'].append(float(info[0]))
                        duration['speech_clear_run3'].append(float(info[1]))
                    elif info[2] == 'speech_noise':
                        onset['speech_noise_run3'].append(float(info[0]))
                        duration['speech_noise_run3'].append(float(info[1]))
                    elif info[2] == 'speaker':
                        onset['speaker_clear_run3'].append(float(info[0]))
                        duration['speaker_clear_run3'].append(float(info[1]))
                    elif info[2] == 'speaker_noise':
                        onset['speaker_noise_run3'].append(float(info[0]))
                        duration['speaker_noise_run3'].append(float(info[1]))
                    elif info[2] == 'instruction':
                        onset['instruction_run3'].append(float(info[0]))
                        duration['instruction_run3'].append(float(info[1]))
            elif run == '04' and subject_id not in ['sub-01', 'sub-02']:
                for line in f:
                    info = line.strip().split()
                    if info[2] == 'speech':
                        onset['speech_clear_run4'].append(float(info[0]))
                        duration['speech_clear_run4'].append(float(info[1]))
                    elif info[2] == 'speech_noise':
                        onset['speech_noise_run4'].append(float(info[0]))
                        duration['speech_noise_run4'].append(float(info[1]))
                    elif info[2] == 'speaker':
                        onset['speaker_clear_run4'].append(float(info[0]))
                        duration['speaker_clear_run4'].append(float(info[1]))
                    elif info[2] == 'speaker_noise':
                        onset['speaker_noise_run4'].append(float(info[0]))
                        duration['speaker_noise_run4'].append(float(info[1]))
                    elif info[2] == 'instruction':
                        onset['instruction_run4'].append(float(info[0]))
                        duration['instruction_run4'].append(float(info[1]))

    from nipype.interfaces.base import Bunch

    # Bunching is done per run, i.e. cond1_run1, cond2_run1, etc.
    subjectinfo = []
    for r in range(len(runs)):
        cond = [c + '_run' + str(r+1) for c in cond_names]
        subjectinfo.insert(r,
                           Bunch(conditions=cond,
                                 onsets=[onset[k] for k in cond],
                                 durations=[duration[k] for k in cond],
                                 amplitudes=None,
                                 tmod=None,
                                 pmod=None,
                                 regressor_names=None,
                                 regressors=None))
    return subjectinfo


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
