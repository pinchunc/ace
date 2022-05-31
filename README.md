# ace (autonomic-central events)

# Please cite 
Naji, M., Krishnan, G. P., McDevitt, E. A., Bazhenov, M., & Mednick, S. C. (2019). Coupling of autonomic and central events during sleep benefits declarative memory consolidation. Neurobiology of learning and memory, 157, 139-150.

Chen, P. C., Whitehurst, L. N., Naji, M., & Mednick, S. C. (2020). Autonomic/central coupling benefits working memory in healthy young adults. Neurobiology of Learning and Memory, 173, 107267.

# To run this code, the following files are required: 
sleep edf file
single column text file for sleep stages (one for every 30 sec)
.mat file of Kubios output of ECG R-R intervals

# Output
Output will be named as .mat
This version of code gives your output per quartile (approximately each overnight sleep cycle).
EEG frequency bands calculated: swa=[0.1 4]; alpha=[8 13]; theta=[4 7]; sigma=[12 16];
EEG events: spindles
