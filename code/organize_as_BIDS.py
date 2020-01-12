#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# Author: Meng Du
# November 2017
# python organize_as_BIDS.py ###


"""
This script renames and reorganizes the fMRI data into the Brain Imaging Data Structure (BIDS).
It's very rudimentary so you may need to change it (a lot) according to your need.

Usage:
        python organize_as_BIDS.py <subject_id_1> <subject_id_2> ...
    OR
        python organize_as_BIDS.py --all

It's not thoroughly tested for all cases -- use it with caution, and run a test before you
use it on your actual files. You can test it by changing the generate_test_files() function
to generate a file structure similar to yours, then running generate_test_files() followed
by main() (see the last two lines of this file), and see if there's anything wrong in the
generated directories/files.

This script also assumes that one task has a latin square design made by subject_id % num_runs.
Change it based on your design.
"""

from __future__ import print_function
import os
import sys
import re
from funcs import *

# DIRECTORY NAMES
BIDS_DIR = '../bids/' #os.environ['BIDS_DIR'] + '/'
SUBID_PREFIX = 'sub-' #os.environ['SUBID_PREFIX']
PATH_BETWEEN_SUBJECT_AND_TASK_DIR = '/raw'

# FILE NAMES
# {current name: (new name, number of runs)}
FUNC_NAME_DICT = {
    'friend': ('encoding_mb_tr750_', 4),
    'number': ('encoding_mb_tr750_', 4)
}
# {scanner label: BIDS label}
ANAT_NAME_DICT = {'MPRAGE': 'T1w'}
FMAP_NAME_DICT = {'SpinEchoFieldMap_': 'epi'}
SBREF_NAME_DICT = {'_SBRef': '_sbref'}

# SCANNER SETTINGS
TOTAL_READOUT_TIME = '0.000580009'  # TODO read EffectiveEchoSpacing from file




def rename(old_item, new_item):
    os.rename(old_item, new_item)
    print('Renamed "%s" to "%s".' % (old_item, new_item))


def get_task_name(path, run_num):
    """
    Finds TSV file from path with specific run number and returns the
    relevant task name. Checks to make sure there is exactly 1 file.

    :param path: string path to participants' functional directory
    :param run_num: string run number to find matching behavioral file. Assumes
                    run_number is full number following 'run-' (i.e., must
                    include leading 0s)
    :param return: the task name for this run as a string
    """
    # get filname for this run (includes run-## and .tsv)
    behav_fname = [f for f in os.listdir(path)
                    if "_run-" + str(run_num) in f and ".tsv" in f]

    # Check for errors
    # check if more than one behavioral file exists for this run. If so, throw error
    if len(behav_fname) > 1:
        raise RuntimeError(
            'More than one matching behavioral file found: %s.' % (behav_fname))
    # check if fewer than one behavioral file exists for this run. If so, throw error
    if len(behav_fname) < 1:
        raise RuntimeError(
            'No behavioral file found: %s.' % (behav_fname))

    # save as string rather than list
    behav_fname = behav_fname[0]
    # split into discrete filename parts
    fname_parts = behav_fname.split('_')
    # find the task section and get relevant task for this run
    for i in fname_parts:
        if "task-" in i:
            run_task = i.split('-')[1]
    # return task name
    return run_task



def rename_func_dirs(path, sid, task_prefix, task_name, run_num_dict=None, multi_run=True):
    """
    Rename the folder names in path that start with task_prefix, based on FUNC_NAME_DICT
    Assuming 1) the last number in the folder name indicates the position of this
                run in the scanning sequence
             2) if multi_run is True, the first number in the folder name indicates
                the number of run, starting from 1
    :param path: string path to data directory for one subject
    :param sid: string subject id
    :param task_prefix: string prefix of the task name
    :param run_num_dict: an optional dictionary of strings {old_run_#: new_run_#},
                         if run # needs to be changed
    :param multi_run: whether the task contains multiple runs (run number 'runX' has to
                      be present in the file name)
    :return: a dictionary {new_folder_name: old_folder_name}
    """

    # get list of all files that start with functional image prefix string
    dir_list = [dir_name for dir_name in sorted(os.listdir(path))
                if dir_name.startswith(task_prefix)]
    # subselect only those that contain r0 (i.e. runs) and are not the SBRef files
    #dir_list = [d for d in dir_list if 'r0' in d and 'SBRef' not in d]  # TODO this is just for the encoding study
    dir_list = [d for d in dir_list if 'r0' in d]

    # remove runs that are not in this condition
    for d in reversed(dir_list):
        run_num = re.search(r'r0\d+', d).group()
        # get condition from tsv file
        run_task = get_task_name(path+'/../func/', run_num[1:3])
        # if this run is not in this task, delete from dir_list
        if run_task != task_name:
            dir_list.remove(d)


    # error checking
    if run_num_dict is not None and len(dir_list) != len(run_num_dict):
        raise RuntimeError(
            'Number of %s folders (%d) does not equal to number of runs (%d).'
            % (task_prefix, len(dir_list), len(run_num_dict)))
    run_order = [int(re.search(r'\d+', d[::-1]).group()[::-1]) for d in dir_list]  # last number

    for i in range(len(run_order) - 2):
        if run_order[i + 2] <= run_order[i]:
            print(run_order)
            raise RuntimeError('Folders are not ordered by time.')
    if multi_run:
        run_ids = [re.search(r'r0\d+', d).group()[1:3] for d in dir_list]  # first number that matches r0#
        # for i, run in enumerate(run_ids):
        #     if i != int(run) and i + 1 != int(run) and i + 4 != int(run) and i + 5 != int(run):
        #         print("WARNING: Something's wrong with the run numbers. \
        #         Check the behavioral tsv files.")
        #         print(run_ids)
        #         #raise RuntimeError('Run numbers in the folder names are not ordered:')

    # renaming
    folder_dict = {}

    # iterate through each file
    for i, folder in enumerate(dir_list):

        # start new filename string
        new_name = 'sub-{0}_task-{1}'.format(sid, task_name)

        # get run number
        if multi_run:
            run_num = run_num_dict[run_ids[i]] if run_num_dict is not None else run_ids[i]
            new_name += '_run-' + run_num.zfill(2)

        # add suffix (e.g., _bold, _sbref)
        func_run=True
        # check if sbref file
        for suffix in SBREF_NAME_DICT:
            if suffix in folder:
                # not the functional run
                func_run=False
                # add approapriate suffix
                new_name += SBREF_NAME_DICT[suffix]
        if func_run:
            # functional run
            new_name += '_bold'

        # rename
        rename(path + folder, path + new_name)
        # add to dictionary
        folder_dict[new_name] = folder

    return folder_dict


def rename_anat_dirs(path, sid):
    """
    Rename the anatomical scan directory based on ANAT_NAME_DICT.
    See rename_func_dirs() for info on parameters and return value.
    """

    folder_dict = {}
    for prefix in ANAT_NAME_DICT:
        new_prefix = ANAT_NAME_DICT[prefix]
        anat_name = None
        # iterate through all files in directory
        for dir_name in os.listdir(path):
            # check if this has prefix
            if dir_name.startswith(prefix):
                anat_name = dir_name
                break
            # check if already renamed
            if dir_name.startswith(new_prefix):
                # if so, simply return the already converted filename
                print('Anatomical already renamed: ' + dir_name)
                folder_dict[dir_name] = dir_name
                return folder_dict
        # if old or new filename was not found, throw error
        if anat_name is None:
            raise RuntimeError('Anatomical scan not found.')

        # rename anatomical
        new_name = 'sub-{0}_'.format(sid) + new_prefix
        rename(path + anat_name, path + new_name)
        # add to dictionary
        folder_dict[new_name] = anat_name
    return folder_dict


def rename_fmap_dirs(path, sid):
    """
    Rename the fieldmap scan directory based on FMAP_NAME_DICT.
    See rename_func_dirs() for info on parameters and return value.
    """
    folder_dict = {}
    # get list of directories
    dir_list = [dir_name for dir_name in os.listdir(path)]
    # iterate through directories
    for prefix in FMAP_NAME_DICT:
        fmap_name = None
        for folder in dir_list:
            # check if starts with prefix
            if folder.startswith(prefix):
                # find last underscore to cutoff numbers the scanner adds
                iend = folder.rfind('_')
                # extract the direction (AP or PA)
                direction = folder[len(prefix):iend]
                # create new name string
                fmap_name = 'sub-{0}_dir-{1}_'.format(sid, direction) + FMAP_NAME_DICT[prefix]
                # rename file
                rename(path + folder, path + fmap_name)
                # add to dictionary
                folder_dict[fmap_name] = folder

            # check if already renamed
            elif folder.startswith(FMAP_NAME_DICT[prefix]):
                fmap_name = folder
                # if so, simply return the already converted filename
                print('Fieldmap already renamed: ' + folder)
                folder_dict[folder] = folder
    # check if folder_dict is empty
    if fmap_name is None:
        # if old or new name was not found, throw error
        raise RuntimeError('Fieldmap scan not found: ' + folder)
    return folder_dict


def rename_files(path, folder_name_dict):
    """
    Rename all files according to their parent folder name,
    based on folder_name_dict
    :param path: string path to the folders where files need renaming
    :param folder_name_dict: {current_folder_name: old_name}
    """
    for folder_name in folder_name_dict:
        old_name = folder_name_dict[folder_name]
        for filename in os.listdir(path + folder_name):
            if filename.startswith(old_name):
                new_name = filename.replace(old_name, folder_name, 1)
                rename(path + folder_name + '/' + filename, path + folder_name + '/' + new_name)


def reorganize_files(subj_dir, sid, dir_list, file_extensions=('.json', '.nii.gz')):
    """
    Reorganize files into BIDS (Brain Imaging Data Structure), i.e. move data from functional
    scans, anatomical scans and fieldmaps to sub-<id>/func, sub-<id>/anat, sub-<id>/fmap,
    respectively. Only files that match <parent_folder_name>.<extension> (where <extension>
    has to be one of the strings specified in file_extensions parameter) are moved. The other
    files stay at where they are.
    :param subj_dir: string path to the directory of a subject that needs to be reorganized
    :param sid: string subject id
    :param dir_list: a list of directory names where the files are (i.e. directories in
                     subj_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR)
    :param file_extensions: a list of file extensions that need to be moved
    """
    dir_lists = {'func': [], 'anat': [], 'fmap': []}
    for folder in dir_list:
        if 'bold' in folder:
            dir_lists['func'].append(folder)
        elif any(postfix in folder for postfix in ANAT_NAME_DICT.values()):
            dir_lists['anat'].append(folder)
        elif any(postfix in folder for postfix in FMAP_NAME_DICT.values()):
            dir_lists['fmap'].append(folder)
        elif 'sbref' in folder:          
            dir_lists['func'].append(folder)

    data_path = subj_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR + '/'
    for dir_type in dir_lists:
        if not os.path.isdir(subj_dir + dir_type):
            os.makedirs(subj_dir + dir_type)
        for folder in dir_lists[dir_type]:
            for f in os.listdir(data_path + folder):
                if any(f == folder + ext for ext in file_extensions):
                    rename(data_path + folder + '/' + f, subj_dir + dir_type + '/' + f)

    rename(subj_dir, BIDS_DIR + 'sub-' + sid)


def append_to_json(filename, contents):
    """
    Append contents to a json file, i.e. insert the content string(s) to the second last line
    (before the "}") of the file.
    :parameter filename: string path and file name of the json
    :parameter contents: (string, or a list of strings) content to insert
    """
    # read
    with open(filename, 'r') as json_file:
        json_content = json_file.readlines()
    # change
    if len(json_content[-2]) > 2 and json_content[-2][-1] != ',':
        json_content[-2] = json_content[-2][:-1] + ',\n'
    if type(contents) is str:
        contents = [contents]
    for line in contents:
        json_content.insert(-1, line)
    # write
    with open(filename, 'w') as json_file:
        json_file.write(''.join(json_content))


def fix_fmap_json(sid, total_readout_time=None):
    """
    Add the "IntendedFor" and optionally "TotalReadoutTime" parameters to the json files for
    fieldmap data so the fieldmaps are intended for all functional scans.
    Assuming the files are already organized as BIDS.
    :parameter sid: string subject id
    :parameter total_readout_time: string or float number, or None if unnecessary
    """
    # get functional scan names
    subject_path = BIDS_DIR + 'sub-%s/' % sid
    func_names = ['func/' + f for f in os.listdir(subject_path + 'func/') if f.endswith('bold.nii.gz')]
    intended_for = '\t"IntendedFor": ["' + '",\n\t\t"'.join(func_names) + '"\n\t],\n'

    # change json
    json_filenames = [f for f in os.listdir(subject_path + 'fmap/') if f.endswith('json')]
    for json_name in json_filenames:
        json_name = subject_path + 'fmap/' + json_name
        contents = intended_for if total_readout_time is None \
                   else [intended_for, '\t"TotalReadoutTime": %s\n' % str(total_readout_time)]
        append_to_json(json_name, contents)


def fix_func_json(sid):
    """
    Add the "TaskName" parameter to json files for functional scans.
    Assuming the files are already organized as BIDS.
    :parameter sid: string subject id
    """
    subject_path = BIDS_DIR + 'sub-%s/' % sid
    json_filenames = [f for f in os.listdir(subject_path + 'func/') if f.endswith('json')]
    for json_name in json_filenames:
        task_name = re.search(r'task-\w+_', json_name).group()[5:-1]
        append_to_json(subject_path + 'func/' + json_name, '\t"TaskName": "%s"\n' % task_name)


def generate_test_files(subject_ids):
    """
    Generate a bunch of directories and files to test main().
    """
    try:
        for sid in subject_ids:
            subject_dir = BIDS_DIR + SUBID_PREFIX + str(sid)
            os.makedirs(subject_dir)
            os.makedirs(subject_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR)
            os.makedirs(subject_dir + '/irrelevant_folder')

            # a one-run task
            run_name = FUNC_NAME_DICT.keys()[2] + '_20'
            run_dir = subject_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR + '/' + run_name
            os.makedirs(run_dir)
            for file_postfix in ['.nii.gz', '.json', '_yo.ica', '_sth_else.pdf']:  # arbitrary stuff
                open(run_dir + '/' + run_name + file_postfix, 'a').close()
            with open(run_dir + '/' + run_name + '.json', 'w') as f:
                f.write('{\n\t"ABC": 123\n}\n')
            # anatomical
            anat_name = ANAT_NAME_DICT.keys()[0] + '_17'
            anat_dir = subject_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR + '/' + anat_name
            os.makedirs(anat_dir)
            for file_postfix in ['.nii.gz', '.json', '_yo.nii.gz', '_sth_else.pdf']:  # arbitrary stuff
                open(anat_dir + '/' + anat_name + file_postfix, 'a').close()
            with open(anat_dir + '/' + anat_name + '.json', 'w') as f:
                f.write('{\n\t"ABC": 123\n}\n')
            # fieldmaps
            for d in ('PA', 'AP'):
                fmap_name = FMAP_NAME_DICT.keys()[0] + d + '_5'
                fmap_dir = subject_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR + '/' + fmap_name
                os.makedirs(fmap_dir)
                for file_postfix in ['.nii.gz', '.json', '_yo.nii.gz', '_sth_else.pdf']:  # arbitrary stuff
                    open(fmap_dir + '/' + fmap_name + file_postfix, 'a').close()
                with open(fmap_dir + '/' + fmap_name + '.json', 'w') as f:
                    f.write('{\n\t"ABC": 123\n}\n')

    except OSError as err:
        print('Error when generating test files: {0}'.format(err))
        return


def main():
    try:
        if len(sys.argv) < 2:
            raise RuntimeError()
        subject_ids = sys.argv[1:]
        if len(subject_ids) == 1 and subject_ids[0] == '--all':
            subject_ids = None
    except RuntimeError:
        print('Usage:\n\t\tpython organize_as_BIDS.py <subject_id_1> <subject_id_2> ...'
              '\n\tOR\n\t\tpython organize_as_BIDS.py --all')
        quit()

    print('Running subject(s): ',subject_ids)

    for subj_dir in os.listdir(BIDS_DIR):
        #print(subj_dir)
        # TODO iterate through subject_ids instead (and print errors)
        if not subj_dir.startswith(SUBID_PREFIX):
            #print(subj_dir +' does not start with '+ SUBID_PREFIX)
            continue
        sid = subj_dir[len(SUBID_PREFIX):]
        if subject_ids is not None and sid not in subject_ids:
            #print(sid + ' is not in ',subject_ids)
            continue

        path = BIDS_DIR + subj_dir + PATH_BETWEEN_SUBJECT_AND_TASK_DIR + '/'
        sid = str(sid)

        print(path)
        print(sid)

        folder_dict = {}
        try:
            # rename anatomical files
            folder_dict = rename_anat_dirs(path, sid)
            # rename functional files
            for task_name in FUNC_NAME_DICT:
                num_runs = FUNC_NAME_DICT[task_name][1]
                task_prefix = FUNC_NAME_DICT[task_name][0]
                run_dict = None
                folder_dict.update(rename_func_dirs(path, sid, task_prefix, task_name, run_dict, multi_run=(num_runs > 1)))
            # rename fieldmap files
            folder_dict.update(rename_fmap_dirs(path, sid))
        except RuntimeError as err:
            print('Error in %s:' % subj_dir, err, 'Skipping %s.' % subj_dir)
            # reverse renamed folders
            for item in folder_dict:
                rename(path + item, path + folder_dict[item])
        else:  # no error
            print('Renaming files')
            rename_files(path, folder_dict)
            print('Reorganizing files')
            reorganize_files(BIDS_DIR + subj_dir + '/', sid, folder_dict.keys())
            fix_fmap_json(sid, total_readout_time=TOTAL_READOUT_TIME)
            fix_func_json(sid)
    print('Done')


if __name__ == '__main__':
    #generate_test_files(list(range(101, 103)))
    main()
