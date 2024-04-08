---
title: "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator"
tags:
  - C++
  - Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS)
  - Message Passing Interface (MPI)
  - Computational Materials Science
  - Off-lattice kinetic Monte Carlo
  - Molecular Dynamics
  
authors:
  - name: Stavros Ntioudis
    orcid: 0009-0000-8095-1727
    affiliation: 1
  - name: James P. Ewen
    orcid: 0000-0001-5110-6970
    affiliation: 1
  - name: Daniele Dini
    orcid: 0000-0002-5518-499X
    affiliation: 1
  - name: C. Heath Turner
    orcid: 0000-0002-5707-9480
    affiliation: 2
affiliations:
 - name: Department of Mechanical Engineering, Imperial College London, London, SW7 2BX, United Kingdom
   index: 1
 - name: Department of Chemical and Biological Engineering, University of Alabama, Tuscaloosa, Alabama 35487, United States of America
   index: 2
date: 25 March 2024
bibliography: paper.bib
---

# Summary

Kinetic Monte Carlo (kMC) is an atomistic and stochastic simulation technique that captures the temporal evolution of various systems in Materials Science, Physics, Biology, and Engineering. An array of kMC open-source packages is currently distributed online. Nevertheless, such implementations are typically lattice-based and are mostly suitable for the sake of investigating ordered materials. In this work, we present an easy-to-use and completely lattice-free open-source kMC software called ```PAPRECA``` that enables kMC simulations on amorphous materials or systems characterized by a low degree of crystallinity. ```PAPRECA``` is a massively parallel C++ software using the Message Passing Interface (MPI) protocol and coupled with the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) [@Thompson2022Feb] to enable not only pure kMC runs but additionally hybrid kMC/MD simulations.

# Statement of need

kMC models have been deployed to investigate non-equilibrium dynamics and properties of thin films [@Ntioudis2023Oct], nanoparticles [@Turner2016Apr], quantum dots [@Zhu2007May], semiconductors [@vanderKaap2016Feb], catalysts [@Stamatakis2012Dec], energy-storage devices [@Abbott2022Feb], interstellar grain chemistry [@Cuppen2013Dec], protein folding [@Makarov2001Jun], and enzyme reactions [@Slepoy2008May]. KMC algorithms offer a tradeoff between model accuracy and simulation performance. This is justified by the fact that kMC does not describe atomic vibrations explicitly but evolves the system through discrete elementary processes (e.g., diffusions, depositions, reactions etc.) [@Gillespie1976Dec ; @Fichthorn1991Jul]. Overall, kMC techniques are on one hand less accurate but on the other hand more computationally efficient than other atomistic techniques such as the Density Functional Theory (DFT) or Molecular dynamics (MD). The efficiency of kMC models unlocks the possibility for long-timescale simulations with molecular-level accuracy beyond the limits of DFT and MD (typically a few microseconds [@Trochet2020Mar,@Vakis2018Sep]).

Typically, on-lattice kMC models select atomistic events from a predefined table and execute them on fixed lattice sites [@Andersen2019Apr]. The use of fixed lattice sites contributes to the computational efficiency of on-lattice algorithms but introduces obstacles associated with the study of unordered materials.
Several lattice-based open-source kMC software are identified in the bibliography. Examples include the KMCLib [@Leetmaa2014Sep], lattice_mc [@Morgan2017May], KMC_Lattice v2.0 [@Heiber2019Jan], Excimontec v1.0 [@Heiber2020Sep], MonteCoffee [@Jorgensen2018Sep], and KIMERA [@Martin2020Sep]. Additionally, the Stochastic Parallel PARticle Kinetic Simulator (SPPARKS) [@Mitchell2023May] offers solely on-lattice kMC modeling capabilities, since the only currently available off-lattice solver is a Metropolis Monte Carlo relaxation scheme. Furthermore, a wide range of on-lattice kMC packages such as kmcos [@kmcos], PyCD [@PyCD], VIS-A-VIS [VISAVIS], MulSKIPS [@MulSKIPS], Kimocs [@Kimocs], KSOME [@KSOME], kMCpy [@kMCpy], Morphokinetics [@Morphokinetics] are detected in open-source software repositories.

EON [@Chill2014May] is the only identified off-lattice package distributed under an open-source license. Nonetheless, EON is an Adaptive kinetic Monte Carlo (AkMC) software. AkMC implementations discover as well as store atomic-scale processes (e.g., reactions, diffusions hops) throughout the simulation instead of stochastically selecting transition events from a predefined table [@Henkelman2001Dec]. Such an instance elevates the accuracy of AkMC approaches but decreases their computational efficiency and elevates their implementation complexity compared to predefined table schemes.

To the best of our knowledge, a completely lattice-free kMC code with predefined table of events is currently unavailable in open-source repositories. ```PAPRECA``` aims to fill that gap by providing the scientific community with a general and easy-to-use solution for performing long-timescale atomistic simulations on complex Materials Science, Physics, Biology, and Engineering problems involving amorphous materials or non-crystalline systems.  As far as we are aware, ```PAPRECA``` is the only off-lattice kMC code with predefined table of events available under an open-source license. 

```PAPRECA``` (in its initial version) can capture four distinct classes of predefined transition events: 1) reactions (bonding and scission), 2) depositions (of molecules and atoms), 3) diffusion hops, and 4) monoatomic desorptions. Virtually any system whose temporal evolution can be described by these atomic-scale processes can be effectively simulated by ```PAPRECA```. Example applications of ```PAPRECA``` include but are not limited to: adsorption/desorption on catalytic surfaces, amorphous thin films (e.g., phosphate films, solid electrolyte interphases, oxide layers), modeling self-diffusion of gases. Furthermore, ```PAPRECA``` allows for the extension of the source code to include other classes of transition events (e.g., reaction chains). 

Finally, ```PAPRECA``` is a massively parallelized code (uses the MPI protocol) that couples an off-lattice kMC solver with a MD solver (LAMMPS [@Thompson2022Feb]) to enable hybrid off-lattice kMC/MD simulations with a predefined table of events. Deploying a dual-solver approach can be beneficial for capturing both fast atomic-scale processes (whose frequency is comparable to that of atomic vibrations) through the MD stage as well as elementary events of elevated activation energies via the off-lattice kMC stage.

# Scalability of PAPRECA

The scalability of ```PAPRECA``` was investigated by performing hybrid kMC/MD simulations on thin films grown from the decomposition of lubricant additive tricresyl phosphate (TCP) molecules on an Fe 110 substrate. For further information regarding the system setup refer to the ```PAPRECA``` documentation (Example Applications section) and our previous study [@Ntioudis2023Oct]. 

Two independent scalability tests were performed. The first scalability test was conducted locally, on a personal computer (Intel(R) Core(TM) i9-10980XE CPU @ 3.00GHz, 128 Gb RAM DDR4 \verb|@| 3200 MHz). Four runs were performed with 1, 4, 9, and 16 MPI processes, respectively.
The local tests simulated 1000 ```PAPRECA``` steps. The second scalability test was performed on Imperial College's HPC service (CPU type: Rome, RAM usage: ~1GB per MPI process, filesystem type: GPFS, interconnect: Ethernet). Such scalability test was conducted with the same parameters as the local scalability test but with a different number of total ```PAPRECA``` steps (i.e., 9000 instead of 1000). Five runs were performed with 1, 4, 16, 64, and 144 MPI processes, respectively. Figure \ref{scalability} illustrates the results of both tests:

![Hybrid kMC/MD scalability tests of PAPRECA for TCP on Fe110 conducted on a workstation (left) and on Imperial College HPC (right).](./scalability.jpg){#scalability width="100%"}


Where the speed-up value of N MPI processes was calculated as the $t_N = \frac{T_1}{T_N}$, with $t_1$ being the total walltime of 1 MPI process. It can be observed that in both tests the kMC stage does not scale as effectively as the MD stage (performed in LAMMPS). Nonetheless, the total speed-up (i.e., combine kMC and MD) of a ```PAPRECA``` run is comparable to the MD stage speed-up. This can be justified by the fact that the kMC stages last significantly less than the MD stages regardless of the number of MPI processes. For instance, the total walltimes of the kMC and MD stages of the
64 MPI processes example (second scalability test on Imperial's HPC service) were 0.226 and 2.71 hours, respectively. Overall, improving the scalability of the kMC stage will be prioritized in the upcoming versions of ```PAPRECA```.

# Acknowledgements

S.N. thanks Shell and the Engineering and Physical Sciences Research Council, United Kingdom (EPSRC) for PhD funding via the InFUSE Prosperity Partnership (EP/V038044/1).
J.P.E. acknowledges the support of the Royal Academy of Engineering (RAEng) for support through their Research Fellowship scheme. 
D.D. acknowledges a Shell/RAEng Research Chair in Complex Engineering Interfaces.
We acknowledge computational resources and support provided by the [Imperial College Research Computing Service](http://doi.org/10.14469/hpc/2232).

# References
