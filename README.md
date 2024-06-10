PAPRECA
[![Tests](https://github.com/sntioudis/papreca/actions/workflows/run_ALLtests.yml/badge.svg)](https://github.com/sntioudis/papreca/actions/workflows/run_ALLtests.yml)
[![Compatibility](https://img.shields.io/badge/Earliest%20Compatible%20Lammps%20Version-patch_15Sep2022-BrightGreen)](https://github.com/sntioudis/papreca/actions/runs/9228674083)
[![Docs](https://github.com/sntioudis/papreca/actions/workflows/documentation.yml/badge.svg)](https://github.com/sntioudis/papreca/actions/workflows/documentation.yml)
[![status](https://joss.theoj.org/papers/f00ac3f3856e2f369c96646b66a1581b/status.svg)](https://joss.theoj.org/papers/f00ac3f3856e2f369c96646b66a1581b)
========
PAPRECA stands for PArallel PREdefined CAtalog off-lattice kinetic Monte Carlo.

PAPRECA is a hybrid rejection-free kinetic Monte Carlo/molecular dynamics (kMC/MD) software built around [LAMMPS](https://github.com/lammps/lammps). PAPRECA uses the [Message Passaging Interface (MPI) protocol](https://en.wikipedia.org/wiki/Message_Passing_Interface)
to enable parallel hybrid off-lattice kMC/MD or pure off-lattice kMC studies related to a variety of systems in materials science and engineering (e.g., thin film growth and heterogeneous catalysis).

Open-source and community guidelines
-------------
PAPRECA is an open-source software that is distributed under the terms of the [GNU General Public License, version 2](https://en.wikipedia.org/wiki/GNU_General_Public_License). PAPRECA is designed as an easy-to-modify and easy-to-extend software. For example, new predefined event templates can be added to the original code.

In the spirit of open source, we encourage contributions from other Authors! Kindly use pull requests to contribute to the code after ensuring your local fork passes ALL tests. Furthermore, you can use GitHub issues for other forms of contribution (e.g., bug fixes, code optimization, new examples, corrections to documentation). Note that, development will solely be performed on the main branch. Merges with the release branch will only happen after ensuring that the code on the main branch is ready for production. Appropriate (semantic) versioning tags will be used to link the released version against the relevant commit.

Additionally, enquiries can be directed to the original developer of this software, Stavros Ntioudis, via e-mail: [stavros.ntioudis20@imperial.ac.uk](mailto:stavros.ntioudis20@imperial.ac.uk).

In any case, kindly review our [Code of Conduct](CODE_OF_CONDUCT.md) before participating in our community.

Documentation
-------------
For additional information (installation, dependencies, theory, algorithms, examples, and more) about PAPRECA please visit our [documentation pages](https://sntioudis.github.io/papreca/).

Citing PAPRECA
---------------

If you use PAPRECA for your research, kindly give credit to the PAPRECA [1] and LAMMPS [2] developers by citing the following papers:

[1] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

[2] Thompson, A.P. et al. "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales", Computer Physics Communications, 272, 108171 (2022). https://doi.org/10.1016/j.cpc.2021.108171
