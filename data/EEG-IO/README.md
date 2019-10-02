# [Dataset] EEG-IO 
EEG dataset recording for 20 subjects for involuntary eye blinks (i.e., forced, when an external stimulation was given) collected on OpenBCI headset

## Citing the Dataset
----------------------------

The dataset is freely available for research use. Please cite the following publication if using
<blockquote>
  <p>Mohit Agarwal, Raghupathy Sivakumar<br />
BLINK: A Fully Automated Unsupervised Algorithm for Eye-Blink Detection in EEG Signals<br />
2019 56th Annual Allerton Conference on Communication, Control, and Computing (Allerton). IEEE, 2019.</p>
</blockquote>

## Dataset Specifications
----------------------------
*  **Sampling Frequency:** 250.0 samples per second
*  **Electrodes:** Fp1 and Fp2
*  **Total Subjects:** 20
*  **Experiment:** Subjects were asked to blink when an external stimulation was provided

## Dataset Files
----------------
Raw EEG data is stored as `S<Sub_ID>_data.csv` and corresponding labels are stored in `S<Sub_ID>_labels.csv`

## Reading Data
----------------

A script `read_data.py` is provided to read the dataset 

#### Raw EEG data
Data is stored in a .csv format where column 0 represents time, and column 1 and 2 represents raw EEG potentials (in uV) for channel 1 (Fp1) and channel 2 (Fp2) respectively.

#### Labels
* The first line is `corrupt, <n>` where n represents the total number of corrupt intervals
* n following lines represent the start and end time of corrupt interval in seconds. A value of -1 means until the end.
* The next line is 'blinks' which marks the starting of blinks
* Blinks are arranged as `<blink_time>, <code>` where
  *  `<blink_time>`: middle point of blink in seconds (where the minima appears in EEG)
  * `<code>`: `0` is normal blink, `1` is blink when stimulation was given, `2` is soft blink


## Contact
----------------

For any queries, contact Mohit Agarwal
Email: me.agmohit@gmail.com
