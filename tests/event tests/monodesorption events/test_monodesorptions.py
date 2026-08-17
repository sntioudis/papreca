import sys
import subprocess


def validatePaprecaSteps(steps, event_name):
    if not steps:
        raise RuntimeError(
            f"No {event_name} events were found in papreca.log."
        )

    expected = list(range(steps[0], steps[0] + len(steps)))
    if steps != expected:
        raise RuntimeError(
            f"PAPRECA {event_name} steps are not consecutive/in order: "
            f"{steps}"
        )


def readPaprecaDesorptions(file_path):
    """
    Read monoatomic-desorption events from PAPRECA's rank-0 structured log.

    papreca.log rows have the form:
        step  MonoDesorption  time  parent_id  parent_type
    """
    atom_ids = []
    steps = []

    with open(file_path, "r") as file:
        for line in file:
            fields = line.split()

            if len(fields) < 5 or fields[1:2] != ["MonoDesorption"]:
                continue

            try:
                step = int(fields[0])
                atom_id = int(fields[3])
            except ValueError:
                continue

            steps.append(step)
            atom_ids.append(atom_id)

    validatePaprecaSteps(steps, "MonoDesorption")
    return atom_ids


def readLammpsDumpFrames(file_path):
    """Parse LAMMPS atom-ID frames using explicit ITEM boundaries."""
    frames = []

    with open(file_path, "r") as file:
        while True:
            line = file.readline()
            if not line:
                break

            if line.strip() != "ITEM: TIMESTEP":
                continue

            timestep = int(file.readline().strip())

            line = file.readline()
            if line.strip() != "ITEM: NUMBER OF ATOMS":
                raise RuntimeError(
                    f"Expected NUMBER OF ATOMS after timestep {timestep}."
                )

            number_of_atoms = int(file.readline().strip())

            line = file.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(
                    f"Expected BOX BOUNDS at timestep {timestep}."
                )

            for _ in range(3):
                if not file.readline():
                    raise RuntimeError(
                        f"Unexpected EOF while reading bounds at "
                        f"timestep {timestep}."
                    )

            header = file.readline().strip()
            if not header.startswith("ITEM: ATOMS"):
                raise RuntimeError(
                    f"Expected ATOMS header at timestep {timestep}."
                )

            columns = header.split()[2:]
            if "id" not in columns:
                raise RuntimeError(
                    f"LAMMPS dump has no 'id' column at timestep {timestep}."
                )

            id_col = columns.index("id")
            atom_ids = set()

            for _ in range(number_of_atoms):
                fields = file.readline().split()

                if not fields:
                    raise RuntimeError(
                        f"Unexpected EOF in atom rows at timestep {timestep}."
                    )

                atom_ids.add(int(fields[id_col]))

            if len(atom_ids) != number_of_atoms:
                raise RuntimeError(
                    f"Duplicate/missing atom IDs at timestep {timestep}."
                )

            frames.append((timestep, atom_ids))

    if not frames:
        raise RuntimeError(
            f"No LAMMPS atom-dump frames found in {file_path!r}."
        )

    return frames


def getDesorbedIdsFromLammps(frames):
    """Reconstruct removed atom IDs from consecutive LAMMPS frames."""
    desorbed = []

    for (previous_step, previous_ids), (current_step, current_ids) in zip(
        frames,
        frames[1:],
    ):
        removed = previous_ids - current_ids
        added = current_ids - previous_ids

        if added:
            raise RuntimeError(
                "Mono-desorption test detected atom(s) being added between "
                f"timesteps {previous_step} and {current_step}: "
                f"{sorted(added)}"
            )

        if len(removed) > 1:
            raise RuntimeError(
                "Mono-desorption test expected at most one removed atom "
                f"between timesteps {previous_step} and {current_step}, "
                f"but found {sorted(removed)}"
            )

        if len(removed) == 1:
            desorbed.append(next(iter(removed)))

    return desorbed


def compareEvents(papreca_ids, lammps_ids):
    print(" ")
    print("PARSED EVENT COUNTS")
    print("----------------------------------------------------------------")
    print("PAPRECA event count:", len(papreca_ids))
    print("LAMMPS event count: ", len(lammps_ids))
    print("PAPRECA IDs:", papreca_ids)
    print("LAMMPS IDs: ", lammps_ids)
    print("----------------------------------------------------------------")
    print(" ")

    if len(papreca_ids) != len(lammps_ids):
        print("ERROR: PAPRECA and LAMMPS event counts differ.")
        return 0.0

    if not papreca_ids:
        print("ERROR: No mono-desorption events were detected.")
        return 0.0

    successes = 0

    for i, (papreca_id, lammps_id) in enumerate(
        zip(papreca_ids, lammps_ids)
    ):
        print(
            f"Step {i}: PAPRECA desorbed atom {papreca_id}; "
            f"LAMMPS lost atom {lammps_id}"
        )

        if papreca_id == lammps_id:
            successes += 1

    success = 100.0 * successes / len(papreca_ids)

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        "Monoatomic desorption tests after comparing "
        f"{len(papreca_ids)} chronological events..."
    )
    print("The test success rate was:", success)
    print("----------------------------------------------------------------")
    print(" ")

    return success


def main():
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_monodesorptions.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    print("Running PAPRECA...")

    # Keep merged MPI stdout for diagnostics, but do not use its line order as
    # the event chronology.
    command = (
        "mpiexec "
        + papreca_path
        + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    )

    mpi_return = subprocess.run(command, shell=True)

    if mpi_return.returncode != 0:
        print("Error: PAPRECA did not finish successfully!")
        sys.exit(1)

    print(
        "Papreca finished successfully..."
        "initiating monoatomic desorptions test!"
    )
    print(" ")

    try:
        papreca_ids = readPaprecaDesorptions("papreca.log")
        frames = readLammpsDumpFrames("lammps_full.dat")
        lammps_ids = getDesorbedIdsFromLammps(frames)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        print("ERROR while parsing mono-desorption outputs:")
        print(str(exc))
        sys.exit(1)

    success = compareEvents(papreca_ids, lammps_ids)
    sys.exit(0 if success == 100.0 else 1)


if __name__ == "__main__":
    main()
