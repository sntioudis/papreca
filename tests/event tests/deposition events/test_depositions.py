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


def readPaprecaDepositions(file_path):
    """
    Read deposition site positions from PAPRECA's rank-0 structured log.

    papreca.log rows begin:
        step  Deposition  time  site_x  site_y  site_z  ...
    """
    positions = []
    steps = []

    with open(file_path, "r") as file:
        for line in file:
            fields = line.split()

            if len(fields) < 6 or fields[1:2] != ["Deposition"]:
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

    validatePaprecaSteps(steps, "Deposition")
    return positions


def readLastLammpsAtomFrame(file_path):
    """
    Read the final LAMMPS atom-dump frame.

    Returns rows as dictionaries keyed by the dump column names.
    """
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


def getMoleculeCOGs(rows, number_of_events):
    """Return COGs for deposited molecule IDs 1..number_of_events."""
    required = ("mol", "xu", "yu", "zu")

    for column in required:
        if rows and column not in rows[0]:
            raise RuntimeError(
                f"LAMMPS dump is missing required column {column!r}."
            )

    cogs = []

    for mol_id in range(1, number_of_events + 1):
        coords = []

        for row in rows:
            if int(row["mol"]) == mol_id:
                coords.append(
                    (
                        float(row["xu"]),
                        float(row["yu"]),
                        float(row["zu"]),
                    )
                )

        if not coords:
            raise RuntimeError(
                f"No atoms with deposited molecule ID {mol_id} were found."
            )

        xyz = np.asarray(coords, dtype=float)
        cogs.append(tuple(np.mean(xyz, axis=0)))

    return cogs


def comparePositions(papreca_positions, lammps_positions):
    print(" ")
    print("PARSED EVENT COUNTS")
    print("----------------------------------------------------------------")
    print("PAPRECA deposition count:", len(papreca_positions))
    print("LAMMPS molecule count:   ", len(lammps_positions))
    print("----------------------------------------------------------------")
    print(" ")

    if len(papreca_positions) != len(lammps_positions):
        print("ERROR: PAPRECA and LAMMPS event counts differ.")
        return 0.0

    successes = 0

    for i, (papreca_pos, lammps_pos) in enumerate(
        zip(papreca_positions, lammps_positions)
    ):
        match = np.allclose(
            papreca_pos,
            lammps_pos,
            atol=COORD_ATOL,
            rtol=COORD_RTOL,
        )

        print(
            f"Step {i}: PAPRECA deposition site = {papreca_pos}; "
            f"LAMMPS molecule COG = {lammps_pos}; match={match}"
        )

        if match:
            successes += 1

    success = 100.0 * successes / len(papreca_positions)

    print(" ")
    print("PRINTING TEST SUMMARY")
    print("----------------------------------------------------------------")
    print(
        f"Deposition test after comparing {len(papreca_positions)} "
        "deposition events..."
    )
    print("The test success rate was:", success)
    print("Coordinate atol:", COORD_ATOL)
    print("----------------------------------------------------------------")
    print(" ")

    return success


def main():
    if len(sys.argv) != 2:
        print(
            "Usage: python3(or python) test_depositions.py "
            "path/to/papreca/executable"
        )
        sys.exit(1)

    papreca_path = sys.argv[1]

    print("Running PAPRECA...")

    # Preserve merged MPI stdout for diagnostics, but do not use it as the
    # authoritative event chronology.
    command = (
        "mpiexec "
        + papreca_path
        + "/papreca -in in_kmc.lmp in_kmc.ppc > papreca_full.log"
    )

    mpi_return = subprocess.run(command, shell=True)

    if mpi_return.returncode != 0:
        print("Error: PAPRECA did not finish successfully!")
        sys.exit(1)

    print("Papreca finished successfully...initiating depositions test!")
    print(" ")

    try:
        papreca_positions = readPaprecaDepositions("papreca.log")
        _, rows = readLastLammpsAtomFrame("lammps_full.dat")
        lammps_positions = getMoleculeCOGs(
            rows,
            len(papreca_positions),
        )
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        print("ERROR while parsing deposition outputs:")
        print(str(exc))
        sys.exit(1)

    success = comparePositions(papreca_positions, lammps_positions)
    sys.exit(0 if success == 100.0 else 1)


if __name__ == "__main__":
    main()
