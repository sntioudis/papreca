# Contributing to PAPRECA

Thank you for your interest in contributing to **PAPRECA**. Contributions that improve the code, tests, examples, documentation, and scientific usability of the project are welcome.

Before starting substantial work, please review the open issues and consider opening an issue to discuss your proposal. Early discussion is especially helpful for new features, changes to user-facing commands, architectural changes, and anything that may affect backward compatibility.

## Table of contents

1. [Code of conduct](#code-of-conduct)
2. [Ways to contribute](#ways-to-contribute)
3. [AI-assisted contributions](#ai-assisted-contributions)
4. [Development requirements](#development-requirements)
5. [Getting started](#getting-started)
6. [Branch and commit workflow](#branch-and-commit-workflow)
7. [Coding standards](#coding-standards)
8. [Testing](#testing)
9. [Documentation](#documentation)
10. [Submitting a pull request](#submitting-a-pull-request)
11. [Reporting bugs](#reporting-bugs)
12. [Requesting features](#requesting-features)
13. [Contact](#contact)

## Code of conduct

By participating in this project, you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

Please communicate respectfully, assume good intent, and keep technical discussions focused on improving the project.

## Ways to contribute

Useful contributions include:

- Fixing bugs.
- Improving performance, memory use, or MPI communication.
- Adding or improving tests.
- Adding simulation examples.
- Adding predefined event types, such as deposition, desorption, reaction, diffusion, bond-formation, or bond-breaking events.
- Improving error messages and input validation.
- Fixing or expanding documentation.
- Proposing and implementing new simulation capabilities.

Small corrections may be submitted directly. For larger changes, please open an issue before implementation so that the intended design and scope can be discussed.

## AI-assisted contributions

AI tools may be used to assist development, but the contributor remains fully responsible for every submitted change.

By submitting AI-assisted work, you confirm that you:

- Understand the code and documentation you are submitting.
- Can explain the implementation and the decisions behind it.
- Have reviewed the output for correctness, security, licensing, and project compatibility.
- Have removed fabricated APIs, unnecessary abstractions, irrelevant comments, and unrelated changes.
- Have built and tested the contribution yourself.

Low-effort or poorly understood generated contributions may be closed without further review.

## Development requirements

PAPRECA requires:

- A compiler with **C++17** support.
- CMake **3.16 or later** when using the CMake build.
- An MPI implementation and compiler wrapper, such as Open MPI with `mpicxx`.
- A compatible shared LAMMPS library and its source headers.
- Bash for the supplied build and test scripts.
- Python 3 and NumPy for the predefined-event tests.

Follow the instructions under [`Installation/`](Installation/) to build PAPRECA and its LAMMPS dependency.

## Getting started

### 1. Fork and clone the repository

Fork `sntioudis/papreca` on GitHub, then clone **your fork**:

```bash
git clone https://github.com/<YOUR-GITHUB-USERNAME>/papreca.git
cd papreca
```

Replace `<YOUR-GITHUB-USERNAME>` with your GitHub username.

### 2. Add the upstream repository

```bash
git remote add upstream https://github.com/sntioudis/papreca.git
git fetch upstream
```

Verify the remotes:

```bash
git remote -v
```

`origin` should refer to your fork, and `upstream` should refer to `sntioudis/papreca`.

### 3. Create a branch from the latest `main`

```bash
git switch main
git fetch upstream
git merge --ff-only upstream/main
git push origin main
git switch -c feature/short-description
```

Use a descriptive branch name, for example:

- `fix/deposition-boundary-check`
- `feature/new-diffusion-event`
- `docs/installation-clarification`
- `test/mpi-regression`

Do not commit contribution work directly to `main` or `release`.

### 4. Build PAPRECA

Follow one of the supported build procedures under [`Installation/`](Installation/).

Keep note of these absolute paths because they are used by the tests:

```bash
export PAPRECA_BUILD_DIR=/absolute/path/to/papreca/build/PAPRECA
export LAMMPS_SRC_DIR=/absolute/path/to/lammps/src
export LAMMPS_LIB_DIR=/absolute/path/to/lammps/build
```

`PAPRECA_BUILD_DIR` must be the directory containing the `papreca` executable.

If the LAMMPS shared library is not installed in a standard library location, make it discoverable before running PAPRECA or the tests:

```bash
export LD_LIBRARY_PATH="$LAMMPS_LIB_DIR:${LD_LIBRARY_PATH:-}"
```

## 6. Branch and commit workflow

The repository uses:

* `main` for active development and pull requests.
* `release` for stable releases maintained by the core maintainers.
* Semantic-version tags, such as `v2.0.0`, for released versions.

Pull requests must target `main`. Do not open pull requests against `release` unless a maintainer explicitly requests it.

Keep commits focused and use clear commit messages. Avoid mixing formatting changes, generated files, refactoring, and functional changes in the same commit unless they are inseparable.

### 6.1 Developer Certificate of Origin

PAPRECA uses the Developer Certificate of Origin (DCO) 1.1.

By contributing to PAPRECA, you certify that you have the legal right to submit your contribution under the GNU General Public License, version 2 (GPL-2.0).

Every commit submitted to PAPRECA must include a `Signed-off-by` line:

```text
Signed-off-by: Your Name <your.email@example.com>
```

You can add this automatically when committing with:

```bash
git commit -s
```

The sign-off certifies your agreement with the Developer Certificate of Origin 1.1:

https://developercertificate.org/

If your contribution is made as part of your employment, you are responsible for ensuring that you have authorization from your employer to submit the contribution under PAPRECA's GPL-2.0 license.

Before opening or updating a pull request, rebase your branch onto the latest upstream `main`:

```bash
git fetch upstream
git rebase upstream/main
```

Resolve any conflicts locally and ensure that the resulting commit history is clean before pushing your changes.

If you have already pushed the branch before rebasing, update your fork using:

```bash
git push --force-with-lease
```

Use `--force-with-lease` rather than `--force` to reduce the risk of overwriting changes you do not have locally.

## Coding standards

### General style

- Follow the conventions already used in the surrounding files under [`source/`](source/).
- Use descriptive names and keep functions focused.
- Prefer simple changes over unnecessary abstractions.
- Avoid unrelated cleanup in a functional pull request.
- Do not leave debugging output, commented-out code, temporary files, or generated build artefacts in the repository.
- Document non-trivial algorithms, scientific assumptions, event-selection logic, and MPI communication.
- Do not add comments that merely restate obvious code.

### C++ requirements

- Write standard **C++17**.
- Ensure the project compiles without new warnings.
- Use existing PAPRECA types, helpers, and wrappers where appropriate.
- Follow the ownership and lifetime patterns used by nearby code.
- Avoid introducing undefined behaviour, unchecked indexing, or unsafe resource handling.

### Classes and event types

PAPRECA uses an object-oriented design for its event system.

New event implementations should derive from the appropriate existing event base class and reuse inherited behaviour rather than duplicating it. This requirement applies to event implementations; unrelated utility or infrastructure classes should use the design that best fits their purpose.

Before adding a new class, check whether the change can be implemented cleanly by extending an existing component.

### MPI safety

Changes involving parallel execution must be MPI-safe.

Pay particular attention to:

- Collective operations being called by all required ranks.
- Consistent ordering of collective calls.
- Synchronisation and data ownership.
- Global reductions and broadcasts.
- Deterministic handling of shared simulation state.
- Output that should be written only by rank 0.
- Serial and multi-rank behaviour.

Do not assume that code working with one MPI rank is correct for multiple ranks.

### LAMMPS integration

Use existing PAPRECA wrappers around LAMMPS functionality wherever possible. Avoid bypassing established interfaces without a clear reason.

Changes that depend on particular LAMMPS packages, commands, or versions must be documented.

### Architectural changes

Changes that affect PAPRECA's core architecture, existing abstractions, class hierarchy, component boundaries, or overall design should be discussed with the maintainers before implementation.

Architectural refactoring should not be bundled into a feature or bug-fix pull request unless it is necessary for the change and has been discussed beforehand. If a proposed feature would benefit from a broader redesign, please open an issue first so that the design and scope can be agreed upon before substantial implementation work begins.

### Input files and backward compatibility

Changes that add or modify user-facing functionality may require updates to PAPRECA input commands. Input parsing is primarily implemented in:

- [`source/libraries/PAPRECA/input_file.cpp`](source/libraries/PAPRECA/input_file.cpp)
- [`source/libraries/PAPRECA/input_file.h`](source/libraries/PAPRECA/input_file.h)

Do not intentionally break existing PAPRECA input files or public behaviour without discussing the change with the maintainers first.

When a breaking change is approved, clearly document:

- What changed.
- Why it changed.
- Which users are affected.
- How existing input files or code should be migrated.

## Testing

All relevant tests must pass before a pull request is submitted.

PAPRECA currently has three main test families:

1. Predefined-event tests.
2. Hybrid kMC/MD tests.
3. Source-level tests.

Run commands from the repository root unless the instructions explicitly change directory.

### Predefined-event tests

These tests require Python 3 and NumPy. Choose an integer seed between `0` and `900000000`.

```bash
(
  cd "tests/event tests"
  bash test_all_events.sh 49123 "$PAPRECA_BUILD_DIR" python3
)
```

The test script creates working `in_kmc.ppc` files from templates. Review `git status` after running the tests and do not commit generated test files unless they are intentionally part of the contribution.

### Hybrid kMC/MD test

```bash
(
  cd "tests/hybrid kMC MD"
  bash test_hybrid.sh "$PAPRECA_BUILD_DIR"
)
```

### Source-level tests

Configure and build the source-test executable:

```bash
cmake \
  -S "tests/source tests" \
  -B "tests/source tests/build" \
  -DLAMMPS_SRC_DIR="$LAMMPS_SRC_DIR" \
  -DLAMMPS_LIB_DIR="$LAMMPS_LIB_DIR" \
  -DPAPRECA_SRC_DIR="$PWD/source"

cmake --build "tests/source tests/build"
```

Run it from the source-test directory so that its relative input paths resolve correctly:

```bash
(
  cd "tests/source tests"
  bash run_tests.sh "$PWD/build"
)
```

### Test expectations for changes

- Bug fixes should include a regression test whenever practical.
- New features and event types must include tests that demonstrate the new behaviour.
- Changes affecting MPI execution should be tested with one rank and multiple ranks where applicable.
- Changes affecting input parsing should include valid-input and invalid-input coverage where practical.
- Tests should be deterministic, isolated, and reasonably fast.
- Do not weaken, skip, or delete an existing test solely to make a change pass.

In the pull-request description, list the exact commands you ran and their results. If a test cannot be run locally, explain why rather than claiming that all tests passed.

## Documentation

Update documentation whenever a contribution changes user-visible behaviour.

Relevant documentation is stored under [`documentation/`](documentation/) and published on the PAPRECA documentation website.

Documentation updates may include:

- New or changed input commands.
- New event types.
- New dependencies or installation steps.
- New examples.
- Behavioural or compatibility changes.
- Scientific assumptions, limitations, or units.

Small documentation corrections can be submitted directly. For larger tutorials or theory sections, please open an issue first to discuss scope and placement.

## Submitting a pull request

Before opening a pull request:

- Rebase your branch onto the latest `upstream/main`.
- Build PAPRECA successfully.
- Run all relevant tests.
- Review the complete diff.
- Remove accidental, generated, or unrelated changes.
- Update documentation and examples where required.
- Confirm that no secrets, credentials, personal paths, or large binary artefacts are included.

Open the pull request against the `main` branch of `sntioudis/papreca`.

The pull-request description should include:

- A concise summary of the change.
- The problem or motivation.
- The implementation approach.
- Any user-visible or compatibility impact.
- Tests performed, including the commands and results.
- Documentation updated.
- Related issues, using wording such as `Closes #42` where appropriate.
- Any remaining limitations, follow-up work, or areas where reviewer attention is requested.

A useful template is:

````markdown
## Summary

Describe what changed.

## Motivation

Explain the problem or requested capability.

## Implementation

Summarise the technical approach and important design decisions.

## Testing

- [ ] PAPRECA builds successfully
- [ ] Predefined-event tests pass, where applicable
- [ ] Hybrid kMC/MD test passes, where applicable
- [ ] Source-level tests pass, where applicable
- [ ] Serial and MPI behaviour tested, where applicable

Commands run:

```text
Paste the exact commands and results here.
```

## Documentation

Describe the documentation or examples updated, or explain why none were needed.

## Compatibility

Describe any effect on existing input files, APIs, output, or scientific behaviour.

## Related issues

Closes #<issue-number>
````

Pull requests may be asked to change scope, add tests, improve documentation, or revise the implementation. Please respond to review comments and keep the branch up to date while review is ongoing.

Maintainers may close pull requests that are abandoned, unrelated to the project, insufficiently tested, or not understood by their author.

## Reporting bugs

Use GitHub Issues for reproducible bugs.

Before opening an issue, search existing issues to avoid duplicates. Include:

- PAPRECA version, branch, tag, or commit hash.
- LAMMPS version or commit.
- Operating system.
- Compiler and version.
- MPI implementation and version.
- Build method and relevant CMake options.
- Number of MPI ranks used.
- Minimal input files or a minimal reproducible example.
- Exact command executed.
- Expected behaviour.
- Actual behaviour.
- Complete error output or stack trace.

Please paste text logs as text or attach them as files rather than providing only screenshots.

Do not include credentials, private data, or information you are not authorised to share.

## Requesting features

Use GitHub Issues to propose features or significant design changes.

Describe:

- The scientific or technical problem.
- The desired behaviour.
- A representative use case.
- Possible alternatives or workarounds.
- Potential effects on input syntax, output, performance, MPI execution, LAMMPS integration, and backward compatibility.

Please wait for maintainer feedback before investing substantial effort in a large implementation.

## Contact

Public technical questions should normally be raised in a GitHub Issue so that answers can benefit other users and contributors.

For matters that are not suitable for a public issue, contact:

**Stavros Ntioudis**  
stavros.ntioudis20@imperial.ac.uk

Thank you for helping improve PAPRECA.
