# CAR T-Cell MFA
 
MFA from isotope tracing in CAR-T Cells

Code and data used in:

Extracellular Domains of CAR Reprogram T-Cell Metabolism Without Antigen Stimulation

Aliya Lakhani, Ximin Chen, Laurence C. Chen, Mobina Khericha, Yvonne Y. Chen, and Junyoung O. Park

University of California, Los Angeles

In the CAR_T_Cell folder,

1. EGFRt.xlsx and Rituximab.xlsx contain measured labeling fractions and net/exchange flux constraints that went into the model as input
2. .xml shows the metabolites, reactions, and carbon mapping
3. .m functions contain the stoichiometric matrix, its kernel, and the EMU (elementary metabolite unite) model
4. .mat file contains stoichiometric matrix, its kernel, and measured metabolite labeling


To run MFA,

1. Open Matlab 2019b or newer
2. Set Path -> Add Folder -> Choose ./src -> Save
3. Open CARTscript.m and execute
4. To run it parallel, enter 'parpool('Processes')' on the command line. MATLAB initiates a pool with one worker per physical core on the local machine
5. Running on clusters requires an additional setup. Please contact system administrator