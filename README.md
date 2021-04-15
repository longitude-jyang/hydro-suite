# hydro-suite
A suite of codes for hydrodynamic analysis of offshore slender structures. You are advised to start with [Introduction](/docs/ASuiteOfCodeForDynamicModellingOfSlenderOffshoreStructures.pdf) to get the big picture of what is inside the code 

This suite consists of eight nested codes of different fidelities (see the possible combination of arrow routes in the chart below). For example, the simplest one is for regular harmonic wave ➡️ rigid body stcuture ➡️ linear frequency response analysis; while the most complicated one is random wave ➡️ flexible structure ➡️ nonlinear time domain analysis. 

![suite of codes flowchart](/docs/codesuite.png) 

 
There are total 35 matlab codes in the folder, see the dependency graph below. Everything starts with [callcode](/Codes/callcode.m). Note that in addition to basic Matlab, there are two files, [cal_BM_AC](/Codes/cal_BM_AC.m) and [solve4RandomTimeResponse](/Codes/solve4RandomTimeResponse.m), that will need Singal Processing Toolbox becaue the function pwelch is used to convert time data to power spectra density. 

![code dependency](/docs/HydroSuiteCodeDependencyGraph.png)
