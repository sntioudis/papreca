# Introduction

\anchor main

\section intro What is PAPRECA?

%PAPRECA off-lattice kMC stands for Parallel Predefined Catalog off-lattice kinetic Monte Carlo.

%PAPRECA is a hybrid rejection-free kinetic Monte Carlo/molecular dynamics (kMC/MD) code. It can be used to capture the long-timescale dynamics
of systems with atomic-level resolution by performing a series of elementary processes (e.g., depositions/desorptions, reactions, diffusion steps).
%PAPRECA can be relevant to scientists conducting research in the broad field of Materials Science and Engineering.
Example applications of %PAPRECA include but are not limited to catalysts, amorphous thin films (e.g., phosphate films, solid electrolyte interphases, oxide layers), modeling self-diffusion of gases.

A typical %PAPRECA step comprises two distinct parts: a kMC stage and an MD stage [1], [2].
During the kMC stage, system sites are scanned to detect atomistic events. Then, a single event is stochastically selected (based on a rejection-free kMC scheme), executed, and finally, the simulation clock is advanced [3].
The kMC stage of the code is completely lattice-free and does not simulate events on "fixed lattice sites". Rather, all particles (e.g., atoms) in the simulation box are treated as "mobile event sites".
Furthermore, a predefined catalog approach is implemented, i.e., kMC events  are selected from event templates defined before the start of the simulation.
The kMC stage is usually followed by the MD stage that allows for the relaxation of interatomic forces or/and the equilibration of the system to the target ambient conditions (i.e., temperature, pressure, etc.).
Of course, by choosing not to invoke the MD stage, %PAPRECA can be used as a pure off-lattice kMC simulator.

%PAPRECA is an extension of the [Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS)](https://docs.lammps.org/Manual.html#) [4]. Native LAMMPS functions as well as %PAPRECA wrappers around LAMMPS functions are utilized for
domain decomposition, treatment of periodic boundaries, update of neighbors lists, addition/deletion of particles/bonds, and modification of particle properties (e.g., positions, velocities, charges).
Also, the MD stage is simulated entirely using the LAMMPS library. Finally, detection, selection, and execution of predefined kMC events, system clock advancement, and communication of event data among the MPI processes and between the kMC and MD stages 
is performed by the %PAPRECA library and the %PAPRECA driver code (papreca.cpp).

%PAPRECA is written entirely in C++ and uses the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface#:~:text=MPI%20%22is%20a%20message%2Dpassing,in%20high%2Dperformance%20computing%20today.) protocol to enable multi-processor hybrid kMC/MD simulations. During the kMC stage, each subdomain (MPI process) discovers predefined kMC events on local (per MPI process) "mobile event sites" (i.e., atoms owned by the respective MPI process).
On the other hand, the efficient parallelization techniques of LAMMPS are utilized during the MD stage.

%PAPRECA is effectively designed for parallel computers but it can also run in serial. %PAPRECA can be built on any LINUX machine that supports an MPI/C++ compiler (e.g., mpicxx or mpiCC) that is at least compatible with the C++-11 standard.

\section authors Authors of PAPRECA and contact details

The following people are the main/original contributors of the %PAPRECA project:

- Stavros Ntioudis ([stavros.ntioudis20@imperial.ac.uk](mailto:stavros.ntioudis20@imperial.ac.uk)): Development and maintenance of the %PAPRECA C++ code, parallel algorithms, conceptualization as well as methodology of the hybrid kMC/MD model, and %PAPRECA website design/maintenance.
- Daniele Dini ([d.dini@imperial.ac.uk](mailto:d.dini@imperial.ac.uk)): Methodology of the hybrid kMC/MD model.
- James P. Ewen ([j.ewen@imperial.ac.uk](mailto:j.ewen@imperial.ac.uk)): Methodology of the hybrid kMC/MD model.
- C. Heath Turner ([hturner@eng.ua.edu](mailto:hturner@eng.ua.edu)): Conceptualization and methodology of the hybrid kMC/MD model.

\section open_source Open-source and community guidelines

%PAPRECA is an open-source software that is distributed under the terms of the [GNU General Public License, version 2](https://en.wikipedia.org/wiki/GNU_General_Public_License).
%PAPRECA is designed as an easy-to-modify and easy-to-extend software. For example, new predefined event templates can be added to the original code.

In the spirit of open source, please let us know if you have any ideas about new features or if you can contribute to the code in any other way (e.g., bug fixes, code optimization, new examples, corrections to documentation).

GitHub users can visit the [PAPRECA repository](https://github.com/sntioudis/papreca) and create:

- GitHub Pull Requests to inform other users about changes to the code.
- GitHub Issues to track ideas, feedback, tasks, or bugs.

Additionally, enquiries can be directed to the main developer of this software, Stavros Ntioudis, via email: [stavros.ntioudis20@imperial.ac.uk](mailto:stavros.ntioudis20@imperial.ac.uk).

\section citations Acknowledgments and citations

If you use %PAPRECA for your research, kindly give credits to the %PAPRECA as well as LAMMPS project authors by citing [1], [2], and [4].

\section bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

[3] Fichthorn, K.A. and Weinberg, W.H. "Theoretical foundations of dynamic Monte Carlo simulations." Journal of Chemical Physics, vol. 95, 1991

[4] Thompson, A.P. et al. "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales." Computer Physics Communications, vol. 272 2022
