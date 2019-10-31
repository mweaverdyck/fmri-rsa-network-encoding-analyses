#!/bin/python
# funcs.py
import os
from os.path import expanduser
HOME = expanduser("~")

# keywords
ALL = os.environ['ALL']
NEW = os.environ['NEW']
PARC = os.environ['PARC']
SL = os.environ['SL']

# arrays that must be manually updated
EXCLUDE_SUBS=["204"]
TASKS=['friend', 'number']
PROCEDURES_ALL=[PARC, SL]

# get project directories
PROJECT_DIR = os.environ['PROJECT_DIR'] + '/'
BIDS_DIR = os.environ['BIDS_DIR'] + '/'
CODE_DIR = os.environ['CODE_DIR'] + '/'
QA_DIR = os.environ['QA_DIR'] + '/'
RECON_DIR = os.environ['RECON_DIR'] + '/'
FMRIPREP_DIR = os.environ['FMRIPREP_DIR'] + '/'
FIRST_LEVEL_DIR = os.environ['FIRST_LEVEL_DIR'] + '/'
SECOND_LEVEL_DIR = os.environ['SECOND_LEVEL_DIR'] + '/'
DERIVATIVES_DIR = os.environ['DERIVATIVES_DIR'] + '/'
GLM_DIR = os.environ['GLM_DIR'] + '/'
RSA_DIR = os.environ['RSA_DIR'] + '/'
MNI_DIR = os.environ['MNI_DIR'] + '/'

# get project variables
SUBID_PREFIX = os.environ['SUBID_PREFIX']
N_NODES = int(os.environ['N_NODES'])
N_TRS = int(os.environ['N_TRS'])
TR = float(os.environ['TR'])
N_TRS_DEL = int(os.environ['N_TRS_DEL'])
N_PARCELS = int(os.environ['N_PARCELS'])
MNI_PARCELLATION = os.environ['MNI_PARCELLATION']
SPACE = os.environ['SPACE']
STAT = os.environ['STAT']
PROCEDURE = os.environ['PROCEDURE']
if PROCEDURE == ALL:
    PROCEDURES = PROCEDURES_ALL
else:
    PROCEDURES = [PROCEDURE]

# filenames
CFD_FNAME = os.environ['CFD_FNAME']
ATTS_FNAME = os.environ['ATTS_FNAME']

# are you on a local computer or Hoffman2
if HOME[:3] == '/u/':
    isLocal=0
else:
    isLocal=1
