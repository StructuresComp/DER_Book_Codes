In this example, we are comparing the shapes of two helices (one from DER and the other one from an "exact" solution). However, comparing two helices aren't trivial- we have to properly align them. This alignment is done in Step 1. Follow the following steps:

(1) Step 1: Run Step1_GetHelixParams.m. It will generate the helixParams.mat file.
(2) Step 2: Run Step2_GenerateFig_ch1_helixTerminalMoment.m. It will generate the figure.
You can generate the datafile "ch1_helixTerminalMomentData.txt" by running problem 2 of the DER code (located in DER_C++, also read DER_Book_2018/README_MAIN.txt for instructions).
The exact solution is contained in "kirchoffrodexact.mat". It can be generated using RodUnderMoment.m.
