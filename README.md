Quantum communication simulation code
This project is a modified version of the original simulator developed by Jianyao Jin during an internship at LQUOM, Inc.

Original Repository:https://gitlab.com/jianyaojin/lquom-simulations

--analystical_solution.py:Author Jianyao Jin---(this is original code)
This code is mathematical simulation by Jianyao Jin.And this code's simulation is based on the thesis Investigating Entanglement Distribution Rates for Multi Platform First Generation Quantum Repeater Networks.
We updated this code to consider Fidelity! So codes described below is inspired by Jianyao Jin and thesis.


from now our contribution


---Montecallo_LQUOM_1usper_sameconditionas_takahasi.py  Author joe---
This code simulate LQUOM simulation as 1 step is 1us.Former simulation is simulated 1 step as 100us. so this simulation is more exact one than other simulation.
and this protocol is simulated samecondition as takahasi protol for the comparing the architecture

---Montecallo_distillitin_value.py :Launched Feb 2026 21:03 Author joe---
This code add distillation efect to Montecallo_simplest_ver.py.We get database which can make simulation time shorter in the case of the long term simulation! and this define 1 step as 100us

---Montecallo_distillition_value_1us_onlymetodA.py Author joe---
This code add distillation efect to Montecallo_simplest_ver.py.We get database which can make simulation time shorter in the case of the long term simulation! and this define 1 step as 1us

---Montecallo_simplest_ver.py:Launched Feb 2026 21:03 Author joe---
This code is simulating quantum communication more strictly than matemathical approach by reproduceing realistic architechre on Python.
we call this simulation "Montecallo" from below comment! But now Method B is inperfect... Please stay until update!!

---Montecallo_single_photon.py Author joe---
this code is almost same condition as Montecallo_simplest_ver.py! but the number of photon when we implement bell measeurement is one.

---Montecallo_single_takahasi.py Author joe---
this code is  simulating quantum communication more strictly than matemathical approach in the ondition of takahasi method by reproduceing realistic architechre on Python.

---Matoryshika.py :Launched May 2 2026 14:59 Author joe---
This code ennables more faster simulation in the situation of caluculating ARC >= 2 than normal simulation by assembling 100000 data from normal simulation in the case of ARC=1.

---Numerical_solution_Qutip_distilliation_value.py  :Launched Feb 2026 21:03 Author joe---
The simulation of distillation circulit.This code is incorporatede to Montecallo_distillitin_value.py.So please use just reference

---lenth_optimazation_analyticaly.py Author joe---
this code simulate lenth optimazation analyticaly.but this method is imperfect now.

---takahasi_analystical_solution_LQUOMver Author joe---
this code is analytical simulation code which smulate takahasi protocol by cordinating LQUOM's simulation

---tau_caluculation.py AUthor joe---
this is code simulate the parameter of tau (https://gitlab.com/jianyaojin/lquom-simulations) for the finding memory time of Matoryosica method.

---Montecalo_single_photon_1us_ion_caluculation.py Author joe---
this code simulates Montecallo and count ion num sametime.

---Montecalo_single_photon_1us_ion_caluculation_ionlimit.py Author joe---
this code simulates Montecallo and count ion num sametime,and consider ion limit

---Matuzaki_ion.py Author joe---
this code simulates Montecallo and count ion num sametime whith easier protocol than Montecalo_single_photon_1us_ion_caluculation.py

---Matuzaki_ion.py Author joe---
this code simulates Montecallo and count ion num sametime whith easier protocol than Montecalo_single_photon_1us_ion_caluculation.py but more stlict condition than Matuzaki_ion.py.








---hybrid(lquom edr saigenn).py:Launched Mar 2026 20:00 Author arata---
This program calculates the EDR by applying LQUOM parameters to the equations of the hybrid repeater model.

