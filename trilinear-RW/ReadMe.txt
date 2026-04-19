#================================================================================================
# *********************************
# INSTRUCTIONS FOR REWEIGHTING CODE
# *********************************
#
# References:
#1. Giuseppe Degrassi, Pier Paolo Giardino, Fabio Maltoni, Davide Pagani (JHEP 1612, 080 (2016))
#2. Fabio Maltoni, Davide Pagani, Ambresh Shivaji and Xiaoran Zhao (arXiv:1707.XXXXX)
#
# Contact: xiaoran.zhao@uclouvain.be,ambresh.shivaji@uclouvain.be
#================================================================================================

With this code following files are provided:

1. 'hhh-model': the UFO model file to be used 
2. "gevirt.sh" : auxiliary script to generate virtual EW supbrocesses
3. "vvh-loop_diagram_generation.py" : to select right set of diagrams in VH, VBF and tHj
4  "tth-loop_diagram_generation.py" : to select right set of diagrams in ttH
5  "check_OLP.f": reweighting code
6. "makefile"   : makefile for reweighting code
7. "example_hz" : contains script ("run_l3_hz.sh") to build ZH reweighting code. Read 
                  comments in the script for the correct usage. 

*******************************************
# Needed : gcc version 4.6.0 or higher is 
           required for loop-calculation 
           in MG5_aMC_v2_5_5
*******************************************

=========================================
Steps to follow: ZH is used as an example  
=========================================

You are inside the 'trilinear-RW' folder.

1. Copy hhh-model in 'MG5_aMC_v2_5_5/models/'. 
   launch './bin/mg5_aMC' in 'MG5_aMC_v2_5_5' and generate 
   ZH process with following syntax, 
----------------------------------
>import model hhh-model
>generate p p > h z [LOonly= QCD]
>output hz_MC
>quit
----------------------------------

2. copy "gevirt.sh" from 'trilinear-RW' in 'MG5_aMC_v2_5_5' and run './gevirt.sh hz_MC'
   (Gives two outputs: "check_olp.inc" & "proc_ml")

3. copy "vvh-loop_diagram_generation.py"  from 'trilinear-RW' in 'madgraph/loop/' and 
   rename it as "loop_diagram_generation.py" 

4. generate EW virtual subprocesses collected in "proc_ml" using 'hhh-model' with output 'hz_ME'
   (DO NOT INSTALL 'collier'. In case you end up installing it, disable it in 
    "MG5_aMC_v2_5_5/input/mg5_configuration.txt" by setting 'collier = None' and remove the # in 
    front of it.)

5. copy following files in 'hz_ME/SubProcesses/'
 "makefile", "check_OLP.f", "check_olp.inc" (provided+generated),
 "pmass.inc", "nsqso_born.inc", "nsquaredSO.inc" (from one of the subprocess folders in 'hz_ME'),
 "c_weight.inc" (from 'hz_MC/SubProcesses') and "nexternal.inc" (from one of the subprocess folders in 'hz_MC')

6. copy "libpdf.a", "libLHAPDF.a", 'Pdfdata', 'PDFsets' from the 'lib' folder of any process already generated 
   in Madgraph to 'hz_ME/lib/'

7. Go to 'hz_ME/SubProcesses/' folder and,  
----------------
make OLP_static
make check_OLP
----------------
   (the output is an executable file 'check_OLP')

8. set 'True = store rwgt info' in "hz_MC/Cards/run_card.dat"

9. generate LO events in 'hz_MC' with following options 
-----------------
fixedorder = OFF
shower     = OFF
reweight   = OFF
order      = LO
madspin    = OFF
-----------------

9. move the LO lhe event file (don't forget to unzip it!) to 'hz_ME/SubProcesses/' and execute './check_OLP'
   (note that the input file name should be "events.lhe"(unweighted) and you get an output file named 
   "events_rwgt.lhe" (weighted).)

10. The steps can be repeated for WH, VBF, tHj and ttH processes.

=====================
# BENCHMARKING
=====================

Using following inputs you should get for ZH 
-----------------------------------------------------------------------------------
default "param_card.dat", 13 TeV, 500K events, lhapdf:90500, Fixed scale: (mh+mz)/2
-----------------------------------------------------------------------------------
 SUM OF ORIGINAL WEIGHTS:   315494.71080292721  (leading order, XLO)
 SUM OF NEW WEIGHTS:   3750.0139794108868       (Only O(\lambda correction, XK3)

 (The weights are cross section(pb) times number of events.)

=> C1 =  XK3/XLO = 0.0119 

========================
# IMPORTANT 
========================

1. For VBF (p p > H j j), DO NOT remove the VH configuration at the diagram level, rather 
   apply suitable kinematic cuts.

2. For ttH, the diagram filter MAY NOT work properly in "gg" channel. Please check that 
   all the loop-diagrams involve only Higgs trilinear self coupling. In case, you find
   diagrams with configuration like g g > H > Htt~, please contact us.

3. For tHj, set top width to zero in "tHj_ME/cards/param_card.dat".
