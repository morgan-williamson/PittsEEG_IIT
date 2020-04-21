# PittsEEG_IIT
code for IIT analysis of Shafto &amp; Pitts (2015) IB paradigm under TLab during Summer 2019/20
there are 2 main sections: analysis on the basis of SCALP EEG recordings, and analysis on the basis of SOURCE reconstructed EEG waveforms
Under SOURCE, the MATLAB pipeline is as follows:   (currently only 4ch IIT analysis)
1. concatenate & binarise the waveforms in <concatenate.m>
2. estimate a per-phase TPM of binarised data in <source_TPM_generator.m>
3. export TPMs to PyPhi on a linux or mac system and run <4ch_IIS_generator.py> to generate sias for all states of each estimated TPM 
4. extract the cause-effect-structures and average over 100/200/400/600ms epochs for SVM decoding in <generating_ces_4ch_source.m>
5. run LIBSVM linear SVM decoding on these ces epochs in <svm_decoding_source_4ch_IIS.m>
6. plot all of these decoding timeseries in <plot_ces_decoding.m>
7. run a control on the frequencies of binarised system states by runnning <svm_source_4ch_states.m> and plotting with <plot_state_decoding.m>

SCALP is pending refinement
