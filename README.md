# Constraint steady-state BN inference

This repository stores source codes and data of the Bachelor thesis of Andrej Šimurka, a student of the Faculty of Informatics, Masaryk University, Brno. The thesis was submitted on 14.12.2022. The submitted version is stored in the branch **bachelor_thesis**. This branch will no longer be updated after the submission of the thesis. The content of the repo:

* data
  * target_networks - 6 Boolean models used for the evaluation in .aeon format
  * input_networks - input data of given 6 models used for the evaluation
    * sinks - steady-states of given 6 models
    * constraints - files for specification of constraints
*src - contains source code of implemented inference algorithm
  * classes - classes that store the majority of the functionality of the algorithm
  * main.py - the entry point of the algorithm
  * ...
  
