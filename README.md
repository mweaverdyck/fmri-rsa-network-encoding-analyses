Neural Encoding of Social Network Position

Started: Fall 2018

Ended:

Behavioral pilot completed November 2018

2 fMRI subjects November 2018 (sub-185, sub-186)

This study builds off of Carolyn Parkinson's Spontaneous Neural Encoding paper (2017).
In the current design, participants undergo 2 sessions.


Session 1:

The first is a behavioral session in which they learn a novel, artificial network of friends.
There are 10 people in the network.
The edges are: [1,2], [1,3], [1,4], [1,9], [2,4], [3,4], [5,6], [5,7], [5,8], [5,9], [6,7], [7,8], [9,0].
Faces from the Chicago Face Database (background altered by Meng Du to look more realistic)
paired with common college-aged names were randomly assigned to each node across participants.
Participants first completed the Pairwise Learning Task, reported each node's ego network,
and drew the full network. If participants passed the ego network (100% accuracy) and the
full network (70% accuracy), then they practiced the fMRI task. If they passed the fMRI task,
they were scheduled for an fmri session 1-7 days from that time.


Session 2:

Participants completed all 4 of the above tasks again (shortened versions) and then were scanned.
Scanning protocol:
8 runs: 4 runs of each question/task (number and friend). Runs were blocked so first 4 was one question,
the last 4 were the other. The order was counterbalanced across subjects (odd vs even subject numbers).

mprage anatomical

HCP functional sequences

TR=750ms

6.5 min runs

Spin Echo EPIs

Auto Align used for all sequences (HCP version)


Analysis:

Convert neuro data to BIDS using custom script (see code/ folder)
Use fmriprep to preprocess data (see settings in run_fmriprep.sh)
Project variables and functions are sourced from funcs
For outline of analysis steps, type:  

    source funcs

    get_outline
