import numpy as np
import re
import sys
import subprocess


def readLammpsDumpFrames(file_path):
    """
    Parse a standard LAMMPS custom dump using the explicit ITEM headers.

    Returns
    -------
    list[tuple[int, set[int]]]
        A list of (timestep, atom_ids) frames.

    This deliberately does NOT infer frame boundaries from atom-ID ordering.
    The atom rows may be sorted or unsorted without affecting the parser.
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
                    "Unexpected end of LAMMPS dump after 'ITEM: TIMESTEP'."
                )

            try:
                timestep = int(timestep_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid timestep value in LAMMPS dump: "
                    f"{timestep_line.strip()!r}"
                ) from exc

            line = file.readline()
            if not line or line.strip() != "ITEM: NUMBER OF ATOMS":
                raise RuntimeError(
                    f"Expected 'ITEM: NUMBER OF ATOMS' after timestep {timestep}."
                )

            number_line = file.readline()
            if not number_line:
                raise RuntimeError(
                    f"Unexpected end of LAMMPS dump while reading the atom count "
                    f"for timestep {timestep}."
                )

            try:
                number_of_atoms = int(number_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid atom count at timestep {timestep}: "
                    f"{number_line.strip()!r}"
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
                        f"Unexpected end of LAMMPS dump while reading {axis} "
                        f"bounds at timestep {timestep}."
                    )

            atom_header = file.readline()
            if not atom_header or not atom_header.startswith("ITEM: ATOMS"):
                raise RuntimeError(
                    f"Expected 'ITEM: ATOMS ...' header at timestep {timestep}."
                )

            columns = atom_header.split()[2:]
            if "id" not in columns:
                raise RuntimeError(
                    f"LAMMPS dump has no 'id' column at timestep {timestep}. "
                    f"Columns were: {columns}"
                )

            id_column = columns.index("id")
            ids = set()

            for atom_index in range(number_of_atoms):
                atom_line = file.readline()
                if not atom_line:
                    raise RuntimeError(
                        f"Unexpected end of LAMMPS dump at timestep {timestep}; "
                        f"expected {number_of_atoms} atom rows but reached EOF "
                        f"after {atom_index}."
                    )

                data = atom_line.split()

                if id_column >= len(data):
                    raise RuntimeError(
                        f"Malformed atom row at timestep {timestep}: "
                        f"{atom_line.strip()!r}"
                    )

                try:
                    atom_id = int(data[id_column])
                except ValueError as exc:
                    raise RuntimeError(
                        f"Invalid atom ID at timestep {timestep}: "
                        f"{data[id_column]!r}"
                    ) from exc

                if atom_id in ids:
                    raise RuntimeError(
                        f"Duplicate atom ID {atom_id} found in LAMMPS dump "
                        f"at timestep {timestep}."
                    )

                ids.add(atom_id)

            if len(ids) != number_of_atoms:
                raise RuntimeError(
                    f"LAMMPS dump frame at timestep {timestep} reports "
                    f"{number_of_atoms} atoms but contains {len(ids)} unique IDs."
                )

            frames.append((timestep, ids))

    if not frames:
        raise RuntimeError(
            f"No LAMMPS dump frames were found in {file_path!r}."
        )

    return frames


def getDesorbedIdsFromLammps(frames):
    """
    Reconstruct the chronological monoatomic-desorption sequence.

    Consecutive LAMMPS dump frames are compared directly. Frames with no
    topology change are ignored. A frame transition that removes more than
    one atom, or adds any atom, is considered invalid for this test.
    """

    id_LAMMPS = []

    for (previous_step, previous_ids), (current_step, current_ids) in zip(
        frames, frames[1:]
    ):
        removed_ids = previous_ids - current_ids
        added_ids = current_ids - previous_ids

        if added_ids:
            raise RuntimeError(
                "Monoatomic-desorption test detected atoms being added between "
                f"LAMMPS timesteps {previous_step} and {current_step}: "
                f"{sorted(added_ids)}"
            )

        # It is valid to have dump frames between events where the atom set
        # does not change.
        if len(removed_ids) == 0:
            continue

        if len(removed_ids) != 1:
            raise RuntimeError(
                "Monoatomic-desorption test expected exactly one atom to be "
                f"removed between LAMMPS timesteps {previous_step} and "
                f"{current_step}, but found {len(removed_ids)}: "
                f"{sorted(removed_ids)}"
            )

        id_LAMMPS.append(next(iter(removed_ids)))

    return np.array(id_LAMMPS, dtype=int)


def getDesorbedIdsFromPapreca(file_path):
    """Read the chronological desorption IDs printed by PAPRECA."""

    id_PAPRECA = []

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith(" Executing "):
                atom_id = re.search(r"ATOM_ID=(\d+)", line)

                if atom_id:
                    id_PAPRECA.append(int(atom_id.group(1)))

    return np.array(id_PAPRECA, dtype=int)


def compareArraysAndPrintStats(id_PAPRECA, id_LAMMPS):
    """
    Compare event chronology exactly.

    IMPORTANT: the final event arrays are intentionally NOT sorted. Sorting
    them would hide a genuine difference in event order.
    """

    print(" ")
    print("PARSED EVENT COUNTS")
    print("----------------------------------------------------------------")
    print("PAPRECA event count:", len(id_PAPRECA))
    print("LAMMPS event count: ", len(id_LAMMPS))
    print("PAPRECA IDs:", id_PAPRECA.tolist())
    print("LAMMPS IDs: ", id_LAMMPS.tolist())
    print("----------------------------------------------------------------")
    print(" ")

    if len(id_PAPRECA) != len(id_LAMMPS):
        print("ERROR: PAPRECA and LAMMPS contain a different number of events.")
        print("The event sequences cannot be compared safely.")
        print(" ")

        print("PRINTING TEST SUMMARY")
        print("----------------------------------------------------------------")
        print(
            "Monoatomic desorption test failed because the two outputs contain "
            "different numbers of detected events."
        )
        print("Successful id comparisons: 0.0")
        print("The test success rate was: 0.0")
        print("----------------------------------------------------------------  \n \n \n")

        return 0.0

    if len(id_PAPRECA) == 0:
        print("ERROR: No monoatomic desorption events were detected.")
        return 0.0

    success_count = 0

    for i, (papreca_id, lammps_id) in enumerate(
        zip(id_PAPRECA, id_LAMMPS)
    ):
        if papreca_id == lammps_id:
            success_count += 1

        print(
            "PAPRECA desorbed atom ID for step ",
            i,
            ": ",
            papreca_id,
        )
        print(
            "LAMMPS desorbed atom ID for step ",
            i,
            ": ",
            lammps_id,
        )

    success = 100.0 * float(success_count) / len(id_PAPRECA)

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        "Monoatomic desorption tests after comparing the ids of ",
        len(id_PAPRECA),
        " events...",
    )
    print("Successful id comparisons: ", success)
    print("The test success rate was: " + str(success))
    print("----------------------------------------------------------------  \n \n \n")

    return success


def main():
    # Check if the correct number of command-line arguments is provided.
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_monodesorptions.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    # Run PAPRECA with MPI and send screen output to papreca_full.log.
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
        "initiating monoatomic desorptions test!"
    )
    print(" ")

    try:
        frames = readLammpsDumpFrames("lammps_full.dat")

        print(
            "Parsed ",
            len(frames),
            " LAMMPS dump frames spanning timesteps ",
            frames[0][0],
            " to ",
            frames[-1][0],
            ".",
        )

        id_LAMMPS = getDesorbedIdsFromLammps(frames)
        id_PAPRECA = getDesorbedIdsFromPapreca("papreca_full.log")

    except (OSError, RuntimeError) as exc:
        print("ERROR while parsing monoatomic-desorption outputs:")
        print(str(exc))
        sys.exit(1)

    success = compareArraysAndPrintStats(id_PAPRECA, id_LAMMPS)

    if success < 100.0:
        sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
