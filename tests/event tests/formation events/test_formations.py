import sys
import subprocess


def canonicalBond(bond_type, atom1, atom2):
    """Return an order-independent representation of an undirected bond."""
    atom_low = min(int(atom1), int(atom2))
    atom_high = max(int(atom1), int(atom2))
    return (int(bond_type), atom_low, atom_high)


def readPaprecaFormationEvents(file_path):
    """
    Read bond-formation events from PAPRECA's rank-0 structured log.

    papreca.log rows have the form:
        step  Bond-form  time  atom1_id  atom2_id  bond_type
    """
    events = []
    steps = []

    with open(file_path, "r") as file:
        for line in file:
            fields = line.split()

            if len(fields) < 6 or fields[1:2] != ["Bond-form"]:
                continue

            try:
                step = int(fields[0])
                atom1 = int(fields[3])
                atom2 = int(fields[4])
                bond_type = int(fields[5])
            except ValueError:
                continue

            steps.append(step)
            events.append(canonicalBond(bond_type, atom1, atom2))

    validatePaprecaSteps(steps, "Bond-form")
    return events


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


def readLammpsLocalDumpFrames(file_path):
    """
    Parse a LAMMPS dump-local file using explicit ITEM frame boundaries.

    Returns a list of:
        (timestep, set[(bond_type, atom_low, atom_high), ...])

    Row ordering inside a LAMMPS local dump is intentionally ignored.
    """
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
            if line.strip() != "ITEM: NUMBER OF ENTRIES":
                raise RuntimeError(
                    f"Expected NUMBER OF ENTRIES after timestep {timestep}."
                )

            number_of_entries = int(file.readline().strip())

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
            if not header.startswith("ITEM: ENTRIES"):
                raise RuntimeError(
                    f"Expected ENTRIES header at timestep {timestep}."
                )

            columns = header.split()[2:]
            required = (
                "c_bondsinfo[1]",
                "c_bondsinfo[2]",
                "c_bondsinfo[3]",
            )

            missing = [column for column in required if column not in columns]
            if missing:
                raise RuntimeError(
                    f"Missing bond columns at timestep {timestep}: {missing}"
                )

            btype_col = columns.index("c_bondsinfo[1]")
            atom1_col = columns.index("c_bondsinfo[2]")
            atom2_col = columns.index("c_bondsinfo[3]")

            bonds = set()

            for _ in range(number_of_entries):
                fields = file.readline().split()
                if not fields:
                    raise RuntimeError(
                        f"Unexpected EOF in entries at timestep {timestep}."
                    )

                bond_type = int(float(fields[btype_col]))
                atom1 = int(float(fields[atom1_col]))
                atom2 = int(float(fields[atom2_col]))

                bonds.add(canonicalBond(bond_type, atom1, atom2))

            frames.append((timestep, bonds))

    if not frames:
        raise RuntimeError(
            f"No LAMMPS bond-dump frames were found in {file_path!r}."
        )

    return frames


def getFormedBondsFromLammps(frames):
    """Reconstruct one bond formation from each changing frame transition."""
    formed_bonds = []

    previous_step, previous_bonds = frames[0]

    # If the first frame is non-empty, those are pre-existing bonds and are
    # not formation events from this run.
    for current_step, current_bonds in frames[1:]:
        removed_bonds = previous_bonds - current_bonds
        added_bonds = current_bonds - previous_bonds

        if removed_bonds:
            raise RuntimeError(
                "Formation test detected bond(s) disappearing between "
                f"timesteps {previous_step} and {current_step}: "
                f"{sorted(removed_bonds)}"
            )

        if len(added_bonds) > 1:
            raise RuntimeError(
                "Formation test expected at most one newly formed bond "
                f"between timesteps {previous_step} and {current_step}, "
                f"but found {sorted(added_bonds)}"
            )

        if len(added_bonds) == 1:
            formed_bonds.append(next(iter(added_bonds)))

        previous_step = current_step
        previous_bonds = current_bonds

    return formed_bonds


def compareEvents(papreca_bonds, lammps_bonds):
    print(" ")
    print("PARSED EVENT COUNTS")
    print("----------------------------------------------------------------")
    print("PAPRECA event count:", len(papreca_bonds))
    print("LAMMPS event count: ", len(lammps_bonds))
    print("PAPRECA bonds:", papreca_bonds)
    print("LAMMPS bonds: ", lammps_bonds)
    print("----------------------------------------------------------------")
    print(" ")

    if len(papreca_bonds) != len(lammps_bonds):
        print("ERROR: PAPRECA and LAMMPS event counts differ.")
        return 0.0

    if not papreca_bonds:
        print("ERROR: No formation events were detected.")
        return 0.0

    successes = 0

    for i, (papreca_bond, lammps_bond) in enumerate(
        zip(papreca_bonds, lammps_bonds)
    ):
        print(
            f"Step {i}: PAPRECA formed {papreca_bond}; "
            f"LAMMPS newly contains {lammps_bond}"
        )

        if papreca_bond == lammps_bond:
            successes += 1

    success = 100.0 * successes / len(papreca_bonds)

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        "Bond formation tests after comparing "
        f"{len(papreca_bonds)} chronological events..."
    )
    print("The test success rate was:", success)
    print("----------------------------------------------------------------")
    print(" ")

    return success


def main():
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_formations.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    print("Running PAPRECA...")

    # papreca_full.log remains useful diagnostic console output, but event
    # chronology is read from the structured rank-0 papreca.log below.
    command = (
        "mpiexec "
        + papreca_path
        + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    )

    mpi_return = subprocess.run(command, shell=True)

    if mpi_return.returncode != 0:
        print("Error: PAPRECA did not finish successfully!")
        sys.exit(1)

    print("Papreca finished successfully...initiating formations test!")
    print(" ")

    try:
        papreca_bonds = readPaprecaFormationEvents("papreca.log")
        frames = readLammpsLocalDumpFrames("bonds_full.log")
        lammps_bonds = getFormedBondsFromLammps(frames)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        print("ERROR while parsing formation outputs:")
        print(str(exc))
        sys.exit(1)

    success = compareEvents(papreca_bonds, lammps_bonds)
    sys.exit(0 if success == 100.0 else 1)


if __name__ == "__main__":
    main()
