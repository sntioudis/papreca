# Contributing to PAPRECA

Thank you for your interest in contributing to **PAPRECA**! Community contributions are warmly welcomed and greatly appreciated. Please take a moment to read these guidelines before submitting your contribution.

Take a look at this repository's open issues. There is usually something you could help with!

---

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [A Note on AI-Generated Contributions](#a-note-on-ai-generated-contributions)
3. [How Can I Contribute?](#how-can-i-contribute)
4. [Getting Started](#getting-started)
5. [Branch Strategy](#branch-strategy)
6. [Coding Standards](#coding-standards)
7. [Testing Requirements](#testing-requirements)
8. [Submitting a Pull Request](#submitting-a-pull-request)
9. [Reporting Bugs and Requesting Features](#reporting-bugs-and-requesting-features)
10. [Documentation Contributions](#documentation-contributions)
11. [Citing PAPRECA](#citing-papreca)
12. [Contact](#contact)

---

## Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md). Please review it before engaging with the community.

---

## A Note on AI-Generated Contributions

AI tools can be a powerful aid in development, and contributors are welcome to use them. 
However, any code submitted (AI-assisted or not) must be fully understood 
and explainable by the contributor. Contributors are expected to be able 
to walk through their code and justify their implementation decisions. The 
contributor is solely responsible for the quality of their submission. Low effort or 
poorly understood AI-generated contributions that do not meet the standards of this 
project will be rejected, and the accountability for such submissions lies entirely 
with the contributor, not the tool used to generate them.

---

## How Can I Contribute?

- **Bug fixes** - Found a bug? Fix it and submit a pull request.
- **Code optimisation** - Improve performance, memory usage, or MPI communication efficiency.
- **New examples** - Add new simulation examples (e.g., thin film growth, catalysis, diffusion studies).
- **New event templates** - New predefined event types (depositions, desorptions, reactions, diffusion hops, bond formation/breaking) are highly encouraged.
- **Documentation** - Fix typos, clarify explanations, or expand the documentation.
- **New features** - Propose and implement new simulation capabilities.

---

## Getting Started

1. Fork the repository on GitHub.
2. Clone your fork locally:

git clone https://github.com/<your-username>/papreca.git
cd papreca

3. Set up the upstream remote:

git remote add upstream https://github.com/sntioudis/papreca.git

4. Follow the Installation Guide in the Installation/ directory to build PAPRECA with its LAMMPS dependency.
5. Verify your setup by running the existing test suite.

---

## Branch Strategy

- main - Active development branch. All contributions should target this branch.
- release - Stable, production-ready branch. Merges into release are performed only by the core maintainers.
- Semantic versioning tags (e.g., v2.0.0) are used to link released versions to specific commits.

Do not submit pull requests directly to the release branch.

When working on a contribution, create a dedicated feature branch from main:

git checkout main
git pull upstream main
git checkout -b feature/my-new-feature

---

## Coding Standards

PAPRECA is written in C++ (C++11 or later). Please adhere to the following conventions:

- Consistency - Follow the style and conventions present in the existing source files under source/.
- Object-Oriented Programming (OOP) - PAPRECA uses an OOP design. New classes must inherit from the existing base event classes found in source/. Do not duplicate 
functionality that already exists in a parent class, and create a new class only if it is absolutely necessary.
- Naming - Use descriptive, camelCase names for variables and functions. Follow the present naming styles closely.
- Comments - Document non-trivial logic, especially for new event templates or MPI communication routines. Do not add comments to explain trivial logic.
- MPI safety - Ensure that any new code involving parallel execution is MPI-safe. Pay careful attention to process synchronisation, global reductions, and output (handled by Process 0 only).
- LAMMPS integration - Use existing PAPRECA wrappers around LAMMPS functions wherever possible.
- No breaking changes - Contributions should not break backward compatibility with existing PAPRECA input files or APIs without prior discussion.
- PAPRECA input files - Certain changes will require the addition or modification 
of PAPRECA commands. Make sure your features/changes are accessible to the user. 
Input file parsing logic is handled in [input_file.cpp](source/libraries/PAPRECA/input_file.cpp) 
and [input_file.h](source/libraries/PAPRECA/input_file.h).
- Documentation - Document all changes properly in the PAPRECA documentation 
found in the documentation/ directory. This includes new commands, new event 
types, and any changes to existing behaviour.

---

## Testing Requirements

Before submitting a pull request, ensure that your local fork passes all existing tests:

cd tests/
# Follow the test instructions in the tests/ directory README

- If you are adding a new feature or event template, please also add corresponding tests that clearly demonstrate that the code functions as expected.
- Tests should cover both serial and parallel (MPI) execution where applicable.
- Verify that your changes compile cleanly without warnings using a C++11-compatible MPI compiler (e.g., mpicxx or mpiCC).

Pull requests that fail existing tests will not be merged.

---

## Submitting a Pull Request

1. Ensure your branch is up to date with upstream/main:

git fetch upstream
git rebase upstream/main

2. Push your branch to your fork:

git push origin feature/my-new-feature

3. Open a Pull Request against the main branch of sntioudis/papreca.
4. In your PR description, please include:
   - A clear summary of the changes made.
   - The motivation and context for the contribution.
   - Reference to any related GitHub Issues (e.g., Closes #42).
   - Confirmation that all tests pass.
5. Be responsive to review feedback. The maintainers may request changes before merging.

---

## Reporting Bugs and Requesting Features

Please use GitHub Issues to:

- Report bugs - include steps to reproduce, expected vs. actual behaviour, and your build details.
- Suggest new features or improvements.
- Track ongoing tasks or discussions.

When reporting a bug, please provide:

- PAPRECA version (e.g., v2.0.0)
- LAMMPS version used
- Compiler and MPI implementation (e.g., mpicxx, OpenMPI version)
- Operating system
- A minimal reproducible example if possible

---

## Documentation Contributions

Documentation lives in the documentation/ directory and is hosted on the PAPRECA GitHub Pages site at https://sntioudis.github.io/papreca/.

- Corrections to existing documentation can be submitted via a Pull Request.
- For larger documentation additions (e.g., new tutorials or theory sections), please open a GitHub Issue first to discuss the scope.

---

## Citing PAPRECA

If you use PAPRECA in your research, please cite the following:

Ntioudis, S., et al. PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator. Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

Ntioudis, S., et al. A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth. Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

Thompson, A.P. et al. LAMMPS - a flexible simulation tool for particle-based materials modeling. Computer Physics Communications, 272, 108171 (2022). https://doi.org/10.1016/j.cpc.2021.108171

---

## Contact

For questions not suited for a public GitHub Issue, you can reach the main developer directly:

Stavros Ntioudis - stavros.ntioudis20@imperial.ac.uk

---

Thank you for helping make PAPRECA better for the entire computational materials science community!
