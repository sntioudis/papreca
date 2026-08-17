import numpy as np
import re
import sys
import subprocess


def canonicalBond(bond_type, atom1, atom2):
    """
    Return a canonical representation of an undirected bond.

    LAMMPS may report the two atoms of the same physical bond in either
    orientation.  For comparison purposes, (1, 5) and (5, 1) must therefore
    be treated as the same bond.
    """
    atom_low = min(int(atom1), int(atom2))
    atom_high = max(int(atom1), int(atom2))

    return (int(bond_type), atom_low, atom_high)


def readLammpsLocalDumpFrames(file_path):
    """
    Parse a LAMMPS 'dump local' file using its explicit ITEM headers.

    Returns
    -------
    list[tuple[int, set[tuple[int, int, int]]]]
        A list of (timestep, bonds) frames, where each bond is stored as
        (bond_type, min(atom1, atom2), max(atom1, atom2)).

    This deliberately does NOT rely on the ordering of local entries within
    a frame.  LAMMPS does not guarantee a consistent ordering of local data
    from one timestep to the next.
    """

    frames = []

    with open(file_path, "r") as file:
        while True:
            line = file.readline()

            if not line:
                break

            if line.strip() != "ITEM: TIMESTEP":
                continue

            timestep_line = file.readline()
            if not timestep_line:
                raise RuntimeError(
                    "Unexpected end of bonds dump after 'ITEM: TIMESTEP'."
                )

            try:
                timestep = int(timestep_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid timestep value in bonds dump: "
                    f"{timestep_line.strip()!r}"
                ) from exc

            line = file.readline()
            if not line or line.strip() != "ITEM: NUMBER OF ENTRIES":
                raise RuntimeError(
                    f"Expected 'ITEM: NUMBER OF ENTRIES' after timestep "
                    f"{timestep}."
                )

            entries_line = file.readline()
            if not entries_line:
                raise RuntimeError(
                    f"Unexpected end of bonds dump while reading the number "
                    f"of entries for timestep {timestep}."
                )

            try:
                number_of_entries = int(entries_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid entry count at timestep {timestep}: "
                    f"{entries_line.strip()!r}"
                ) from exc

            line = file.readline()
            if not line or not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(
                    f"Expected 'ITEM: BOX BOUNDS' at timestep {timestep}."
                )

            # Standard 3-D LAMMPS dump: one bounds line for x, y, and z.
            for axis in ("x", "y", "z"):
                bounds_line = file.readline()

                if not bounds_line:
                    raise RuntimeError(
                        f"Unexpected end of bonds dump while reading {axis} "
                        f"bounds at timestep {timestep}."
                    )

            entries_header = file.readline()

            if not entries_header or not entries_header.startswith(
                "ITEM: ENTRIES"
            ):
                raise RuntimeError(
                    f"Expected 'ITEM: ENTRIES ...' at timestep {timestep}."
                )

            columns = entries_header.split()[2:]

            required_columns = (
                "c_bondsinfo[1]",
                "c_bondsinfo[2]",
                "c_bondsinfo[3]",
            )

            missing_columns = [
                column for column in required_columns if column not in columns
            ]

            if missing_columns:
                raise RuntimeError(
                    f"Missing required bond columns at timestep {timestep}: "
                    f"{missing_columns}. Columns were: {columns}"
                )

            btype_column = columns.index("c_bondsinfo[1]")
            atom1_column = columns.index("c_bondsinfo[2]")
            atom2_column = columns.index("c_bondsinfo[3]")

            bonds = set()

            for entry_index in range(number_of_entries):
                entry_line = file.readline()

                if not entry_line:
                    raise RuntimeError(
                        f"Unexpected end of bonds dump at timestep {timestep}; "
                        f"expected {number_of_entries} local entries but "
                        f"reached EOF after {entry_index}."
                    )

                data = entry_line.split()

                max_column = max(
                    btype_column,
                    atom1_column,
                    atom2_column,
                )

                if max_column >= len(data):
                    raise RuntimeError(
                        f"Malformed local entry at timestep {timestep}: "
                        f"{entry_line.strip()!r}"
                    )

                try:
                    # int(float(...)) also tolerates integer-valued data if a
                    # dump format happens to write it as e.g. "5.0".
                    bond_type = int(float(data[btype_column]))
                    atom1 = int(float(data[atom1_column]))
                    atom2 = int(float(data[atom2_column]))
                except ValueError as exc:
                    raise RuntimeError(
                        f"Invalid bond entry at timestep {timestep}: "
                        f"{entry_line.strip()!r}"
                    ) from exc

                bonds.add(canonicalBond(bond_type, atom1, atom2))

            frames.append((timestep, bonds))

    if not frames:
        raise RuntimeError(
            f"No LAMMPS local-dump frames were found in {file_path!r}."
        )

    return frames


def getFormedBondsFromLammps(frames):
    """
    Reconstruct the chronological bond-formation sequence from consecutive
    LAMMPS local-dump frames.

    The physical bond set is compared between frames, so the order of rows
    inside a frame is irrelevant.

    For this formation-only test:
      * bonds must never disappear;
      * a frame transition that changes topology should add exactly one bond.
    """

    formed_bonds = []

    # Preserve the semantics of the original test: before the first observed
    # bonds frame, the test assumes that no formation event has yet been
    # accounted for.
    previous_bonds = set()
    previous_step = None

    for current_step, current_bonds in frames:
        removed_bonds = previous_bonds - current_bonds
        added_bonds = current_bonds - previous_bonds

        if removed_bonds:
            raise RuntimeError(
                "Formation test detected bonds disappearing between "
                f"LAMMPS timesteps {previous_step} and {current_step}: "
                f"{sorted(removed_bonds)}"
            )

        # No formation occurred between these observations.
        if len(added_bonds) == 0:
            previous_bonds = current_bonds
            previous_step = current_step
            continue

        # Each kMC formation event in this test is expected to add one bond.
        if len(added_bonds) != 1:
            raise RuntimeError(
                "Formation test expected exactly one newly formed bond "
                f"between LAMMPS timesteps {previous_step} and "
                f"{current_step}, but found {len(added_bonds)}: "
                f"{sorted(added_bonds)}"
            )

        formed_bonds.append(next(iter(added_bonds)))

        previous_bonds = current_bonds
        previous_step = current_step

    return formed_bonds


def getFormedBondsFromPapreca(file_path):
    """
    Read PAPRECA bond-formation events in chronological order.

    Atom IDs are canonicalized so that atom1/atom2 orientation does not affect
    the physical bond comparison.
    """

    formed_bonds = []

    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith(" Executing bond formation event"):
                continue

            bond_type = re.search(r"BOND_TYPE=(\d+)", line)
            atom1_id = re.search(r"ATOM1_ID=(\d+)", line)
            atom2_id = re.search(r"ATOM2_ID=(\d+)", line)

            if bond_type and atom1_id and atom2_id:
                formed_bonds.append(
                    canonicalBond(
                        int(bond_type.group(1)),
                        int(atom1_id.group(1)),
                        int(atom2_id.group(1)),
                    )
                )

    return formed_bonds


def compareArraysAndPrintStats(papreca_bonds, lammps_bonds):
    """
    Compare the chronological bond-formation sequence exactly.

    The list itself is intentionally NOT sorted.  Only the two atom IDs within
    an individual bond are canonicalized.  Sorting the event list would hide a
    genuine difference in event chronology.
    """

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
        print(
            "ERROR: PAPRECA and LAMMPS contain a different number of "
            "bond-formation events."
        )
        print("The event sequences cannot be compared safely.")
        print(" ")

        print("PRINTING TEST SUMMARY")
        print("----------------------------------------------------------------")
        print(
            "Bond formation test failed because the two outputs contain "
            "different numbers of detected formation events."
        )
        print("The average test success rate was: 0.0")
        print("---------------------------------------------------------------- \n \n \n")

        return 0.0

    if len(papreca_bonds) == 0:
        print("ERROR: No bond-formation events were detected.")
        return 0.0

    btype_success = 0
    atom1_success = 0
    atom2_success = 0

    for i, (papreca_bond, lammps_bond) in enumerate(
        zip(papreca_bonds, lammps_bonds)
    ):
        btype_PAPRECA, atom1_PAPRECA, atom2_PAPRECA = papreca_bond
        btype_LAMMPS, atom1_LAMMPS, atom2_LAMMPS = lammps_bond

        if btype_PAPRECA == btype_LAMMPS:
            btype_success += 1

        if atom1_PAPRECA == atom1_LAMMPS:
            atom1_success += 1

        if atom2_PAPRECA == atom2_LAMMPS:
            atom2_success += 1

        print(
            "PAPRECA formed the following bond: ",
            "btype=",
            btype_PAPRECA,
            " atom1_id=",
            atom1_PAPRECA,
            " atom2_id=",
            atom2_PAPRECA,
            " in step ",
            i,
        )

        print(
            "The following bond: ",
            "btype=",
            btype_LAMMPS,
            " atom1_id=",
            atom1_LAMMPS,
            " atom2_id=",
            atom2_LAMMPS,
            " was newly observed in the LAMMPS bond list in step ",
            i,
        )

    number_of_events = len(papreca_bonds)

    btype_success = 100.0 * float(btype_success) / number_of_events
    atom1_success = 100.0 * float(atom1_success) / number_of_events
    atom2_success = 100.0 * float(atom2_success) / number_of_events

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        "Bond formation test results after comparing the bond types, "
        "atom1_ids, and atom2_ids of ",
        number_of_events,
        " PAPRECA events...",
    )
    print("Successful btypes: ", btype_success)
    print("Successful atom1_ids: ", atom1_success)
    print("Successful atom2_ids: ", atom2_success)

    success_avg = (
        btype_success + atom1_success + atom2_success
    ) / 3.0

    print("The average test success rate was: " + str(success_avg))
    print("---------------------------------------------------------------- \n \n \n")

    return success_avg


def main():
    # Check if the correct number of command-line arguments is provided.
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_formations.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    print("Running PAPRECA...")

    command = (
        "mpiexec "
        + papreca_path
        + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    )

    mpi_return = subprocess.run(command, shell=True)

    if mpi_return.returncode != 0:
        print(
            "Error: PAPRECA did not finish successfully! "
            "Please check your papreca executable path"
        )
        sys.exit(1)

    print(
        "Papreca finished successfully..."
        "initiating formations test!"
    )
    print(" ")

    try:
        papreca_bonds = getFormedBondsFromPapreca(
            "papreca_full.log"
        )

        frames = readLammpsLocalDumpFrames(
            "bonds_full.log"
        )

        print(
            "Parsed ",
            len(frames),
            " LAMMPS bond-dump frames spanning timesteps ",
            frames[0][0],
            " to ",
            frames[-1][0],
            ".",
        )

        lammps_bonds = getFormedBondsFromLammps(frames)

    except (OSError, RuntimeError) as exc:
        print("ERROR while parsing bond-formation outputs:")
        print(str(exc))
        sys.exit(1)

    success = compareArraysAndPrintStats(
        papreca_bonds,
        lammps_bonds,
    )

    if success < 100.0:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
