#!/bin/python
import os
import nibabel as nib
from templateflow.api import get
from funcs import *

parcellations = [f for f in get('MNI152NLin2009cAsym', resolution=2, atlas="Schaefer2018")]

mni_dir = MNI_DIR + 'Schaefer/'

for parc in parcellations:
    img=nib.load(str(parc))
    img.to_filename(mni_dir+parc.name)
    print("Saved " + mni_dir + parc.name)
