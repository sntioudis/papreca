import sys
import subprocess
import numpy as np


# papreca.log currently writes event coordinates to four digits after the
# decimal point.  Use an absolute tolerance consistent with that precision.
COORD_ATOL = 1.0e-4
COORD_RTOL = 1.0e-5


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


def readPaprecaDiffusions(file_path):
    """
    Read diffusion vacancy positions from PAPRECA's rank-0 structured log.

    papreca.log rows begin:
        step  Diffusion  time  vac_x  vac_y  vac_z  ...
    """
    positions = []
    steps = []

    with open(file_path, "r") as file:
        for line in file:
            fields = line.split()

            if len(fields) < 6 or fields[1:2] != ["Diffusion"]:
                continue

            try:
                step = int(fields[0])
                position = (
                    float(fields[3]),
                    float(fields[4]),
                    float(fields[5]),
                )
            except ValueError:
                continue

            steps.append(step)
            positions.append(position)

    validatePaprecaSteps(steps, "Diffusion")
    return positions


def readLastLammpsAtomFrame(file_path):
    """Read the final atom-dump frame as dictionaries keyed by column name."""
    last_frame = None

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
            rows = []

            for _ in range(number_of_atoms):
                fields = file.readline().split()

                if len(fields) != len(columns):
                    raise RuntimeError(
                        f"Malformed atom row at timestep {timestep}: "
                        f"{fields}"
                    )

                rows.append(dict(zip(columns, fields)))

            last_frame = (timestep, rows)

    if last_frame is None:
        raise RuntimeError(
            f"No LAMMPS atom-dump frame found in {file_path!r}."
        )

    return last_frame


def getDiffusedPositions(rows):
    """
    Collect final positions of type-2 atoms.

    Their row order is deliberately ignored later because LAMMPS atom-row
    order is not an event chronology.
    """
    required = ("type", "xu", "yu", "zu")

    for column in required:
        if rows and column not in rows[0]:
            raise RuntimeError(
                f"LAMMPS dump is missing required column {column!r}."
            )

    positions = []

    for row in rows:
        if int(row["type"]) == 2:
            positions.append(
                (
                    float(row["xu"]),
                    float(row["yu"]),
                    float(row["zu"]),
                )
            )

    return positions


def matchPositions(papreca_positions, lammps_positions):
    """
    Match PAPRECA diffusion coordinates to unused LAMMPS type-2 coordinates.

    The final LAMMPS atom dump does not encode event chronology, so comparing
    its row order to PAPRECA event order would be invalid.
    """
    if len(papreca_positions) != len(lammps_positions):
        return [], False

    unused = list(enumerate(lammps_positions))
    matches = []

    for papreca_index, papreca_pos in enumerate(papreca_positions):
        found = None

        for list_index, (lammps_index, lammps_pos) in enumerate(unused):
            if np.allclose(
                papreca_pos,
                lammps_pos,
                atol=COORD_ATOL,
                rtol=COORD_RTOL,
            ):
                found = (list_index, lammps_index, lammps_pos)
                break

        if found is None:
            return matches, False

        list_index, lammps_index, lammps_pos = found
        matches.append(
            (papreca_index, papreca_pos, lammps_index, lammps_pos)
        )
        unused.pop(list_index)

    return matches, len(unused) == 0


def comparePositions(papreca_positions, lammps_positions):
    print(" ")
    print("PARSED EVENT COUNTS")
    print("----------------------------------------------------------------")
    print("PAPRECA diffusion count:", len(papreca_positions))
    print("LAMMPS type-2 count:    ", len(lammps_positions))
    print("----------------------------------------------------------------")
    print(" ")

    matches, success_bool = matchPositions(
        papreca_positions,
        lammps_positions,
    )

    for pap_idx, pap_pos, lmp_idx, lmp_pos in matches:
        print(
            f"PAPRECA event {pap_idx} position {pap_pos} matched "
            f"LAMMPS type-2 row {lmp_idx} position {lmp_pos}"
        )

    success = 100.0 if success_bool else 0.0

    if not success_bool:
        print("ERROR: not all diffusion coordinates could be matched.")
        print("PAPRECA positions:", papreca_positions)
        print("LAMMPS positions: ", lammps_positions)

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        f"Diffusion test after comparing {len(papreca_positions)} "
        "diffusion events..."
    )
    print("The test success rate was:", success)
    print("Coordinate atol:", COORD_ATOL)
    print("----------------------------------------------------------------")
    print(" ")

    return success


def main():
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_diffusions.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    print("Running PAPRECA...")

    # Preserve merged MPI stdout for diagnostics only.  The authoritative
    # PAPRECA event order is read from papreca.log.
    command = (
        "mpiexec "
        + papreca_path
        + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    )

    mpi_return = subprocess.run(command, shell=True)

    if mpi_return.returncode != 0:
        print("Error: PAPRECA did not finish successfully!")
        sys.exit(1)

    print("Papreca finished successfully...initiating diffusions test!")
    print(" ")

    try:
        papreca_positions = readPaprecaDiffusions("papreca.log")
        _, rows = readLastLammpsAtomFrame("lammps_full.dat")
        lammps_positions = getDiffusedPositions(rows)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        print("ERROR while parsing diffusion outputs:")
        print(str(exc))
        sys.exit(1)

    success = comparePositions(papreca_positions, lammps_positions)
    sys.exit(0 if success == 100.0 else 1)


if __name__ == "__main__":
    main()
