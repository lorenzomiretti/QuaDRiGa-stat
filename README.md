# QuaDRiGa-stat

This is a fork of the original QuaDRiGa (short for QUAsi Deterministic RadIo channel GenerAtor) repository (v2.8.1) hosted by the Fraunhofer HHI (https://github.com/fraunhoferhhi/QuaDRiGa). 

I created this fork to develop new functions that allow to use QuaDRiGa also for performing conventional statistical simulations of MIMO fading channels, which typically involve the evaluation of 
statistical performance metrics (such as the ergodic capacity) by averaging over multiple i.i.d. small-scale fading realizations. As the name suggests, QuaDRiGa models the impact of small-scale 
fading in a quasi-deterministic manner, by updating the channel parameters with sofisticated techniques as the terminals move along the service area. While this is indeed the killer feature of 
QuaDRiGa, it is intrinsically incompatible with the abovementioned statistical approaches. Fortunately, many of the inner building blocks of QuaDRiGa can be readily reused to perform statistical simulations 
that closely follow 3GPP guidelines while maintaning the convenient properties of simplified models such as the omnipresent spatially correlated Rayleigh fading model.

I am developing this experimental version mostly for my own research. It is currently in a very early stage of development, and it has been tested only under the restricted setup documented in `QuaDriGa_stat_tutorial.m`. 

QuaDRiGa v2.8.1 is released under the
Software License for The QuaDRiGa Channel Model  
© Copyright 2011 - 2021 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V., All rights reserved.

