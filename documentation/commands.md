# %PAPRECA Commands

\anchor commands

This page contains information regarding %PAPRECA commands. For details related to LAMMPS commands please visit this [documentation page](https://docs.lammps.org/commands_list.html#).

> **Important note::**
> %PAPRECA will always ignore text to the right of a "#" character. You can use "#" characters to add comments to your input files.


\section FIX_papreca fix papreca command

\subsection FIX_papreca_syntax Syntax

```bash
fix papreca all papreca
```

\subsection FIX_papreca_description Description

LAMMPS fix command to be used in the LAMMPS input file. This command has to be utilized after defining the LAMMPS simulation box in your LAMMPS input file. 
The papreca fix for LAMMPS initializes and updates two neighbors lists (a half and a full) and facilitates predefined event detection and execution.
Certain kMC events (e.g., diffusion, deposition) require interference (i.e., collision) checks. Therefore, a full neighbor list has to be used with such events as it guarantees that
the collision can be detected when scanning through the neighbor list of either atom. On the other hand, the discovery of other kMC events (e.g., bond formation) is significantly more efficient when
half lists are used (since those lists each pair of atoms once). Effectively, this is why two neighbor lists are required for %PAPRECA.
Please refer to the [LAMMPS documentation](https://docs.lammps.org/Developer_par_neigh.html) for more information regarding neighbor lists.

> **Note 1:**
> If you plan to use [special_bonds](https://docs.lammps.org/special_bonds.html) in your simulation refrain from setting **ANY** of the special_bonds to zero. Setting a special_bond to zero eliminates the (1-2,1-3, or 1-4) neighbors from the neighbor lists. Please use a double number beyond the accuracy limits of a C++ double instead of zero (e.g., use "special_bonds lj 1e-100 1.0 1.0 coul 1e-100 1.0 1.0" in your input file instead of "special_bonds lj 0.0 1.0 1.0 coul 0.0 1.0 1.0" to include the 1-2 neighbors). Once again, this does not affect the computational efficiency of the MD stage but includes additional neighbor pairs in the neighbor list.

> **Note 2:**
> Consider using [neigh_modify](https://docs.lammps.org/neigh_modify.html) command with the "every 1", "delay 0", and "check yes" options in your LAMMPS input file (where applicable). This forces LAMMPS to update the neighbors lists on every timestep and increases the accuracy of your %PAPRECA run (in the expense of computational efficiency).

<hr>

\section KMC_steps kMC_steps command

\subsection KMCsteps_syntax Syntax

```bash
KMC_steps N
```

- N = Natural number denoting the total number of kMC stages.

\subsection KMCsteps_examples Example(s)

```bash
KMC_steps 10000 #Sets the total number of KMC stages to 10000
```

\subsection KMCsteps_description Description

Sets the total number of kMC stages for the %PAPRECA run.

> **Note:**
> This is a mandatory command. The %PAPRECA simulation will not start unless the total number of kMC stages is set.


<hr>


\section KMC_per_MD KMC_per_MD command

\subsection KMC_per_MD_syntax Syntax

```bash
KMC_per_MD N
```

- N = Natural number denoting the number of KMC stages per MD stage.

\subsection KMC_per_MD_examples Example(s)

```bash
KMC_per_MD 5 #This means that an MD stage will be simulated every 5 kMC stages.
```

\subsection KMC_per_MD_description Description

Sets the total number of kMC stages per MD stage for the %PAPRECA run.

> **Note:**
> This is a mandatory command. The %PAPRECA simulation will not start unless the total number of kMC stages per MD stage is set.


<hr>

\section KMC_per_longMD KMC_per_longMD command

\subsection KMC_per_longMD_syntax Syntax

```bash
KMC_per_longMD N
```

- N = Natural number denoting the number of KMC stages per long MD stage. This number has to be greater than KMC_per_MD as set in \ref KMC_per_MD.

\subsection KMC_per_longMD_examples Example(s)

```bash
KMC_per_longMD 5 #This means that a long MD stage will be simulated every 5 kMC stages.
```

\subsection KMC_per_MD_description Description

Sets the total number of kMC stages per long MD stage for the %PAPRECA run.


<hr>



\section random_seed random_seed command

\subsection random_seed_syntax Syntax

```bash
random_seed N
```

- N = Integer number, where 0 &le; N &le; 900000000.

\subsection random_seed_examples Example(s)

```bash
KMC_per_MD 19032024
```

\subsection random_seed_description Description

Defines the random seed to initialize the random number generator. %PAPRECA uses the following random number generator: RANMAR (the algorithm of Marsaglia, Zaman, Tsang)[1].

> **Note:**
> This is a mandatory command. The %PAPRECA simulation will not start unless the random seed is set.

\subsection random_seed_bibliography Bibliography

[1] James, F. "A review of pseudorandom number generators." Computer Physics Communications, vol. 60, 1990


<hr>


\section sigoptions sigmas_options command

\subsection sigoptions_syntax Syntax

```bash
sigmas_options style keyword values
```

- style = manual or LAMMPS.

```bash
keyword = mix
	mix value = geom or arithm
```

\subsection sigoptions_examples Example(s)

```bash
sigmas_options LAMMPS
sigmas_options manual mix geom
sigmas_options manual mix arithm
```

\subsection sigoptions_description Description

Defines the options related to the initialization of "sigma" (as in the Lennard-Jones potential, see [here](https://en.wikipedia.org/wiki/Lennard-Jones_potential) ) values for the %PAPRECA simulation. Note that, for events requiring collision tests (e.g., deposition events)
the sigma value defines the maximum acceptable distance between two species (i.e., if d &le; sigma %PAPRECA assumes that an interference exists and that atom collides) [1] , [2].

For style manual, the sigma values have to be manually declared in the %PAPRECA input file (using sequential \ref insigma) for all types.

For style LAMMPS, the sigma values are retrieved from LAMMPS. Note that, it is necessary that your [pair_style](https://docs.lammps.org/pair_style.html) of choice uses sigma values (e.g., [lj/cut/](https://docs.lammps.org/pair_lj.html) and that all sigma values are properly set in the LAMMPS input file.

The mix keyword performs geometric or arithmetic mixing on already defined sigma values. The mixing is performed after the last \ref insigma line in the %PAPRECA input (or right before the start of the %PAPRECA run if the \ref insigma has not been used). Hence, it can be very useful if you wish to solely
define the diagonal sigma terms. See this [LAMMPS documentation](https://docs.lammps.org/pair_modify.html) page for more information regarding the mixing styles.

> **Note:**
> The sigma values for all pairs of atom types (diagonal and cross-terms) have to be defined before the start of a %PAPRECA run.  The %PAPRECA run will start normally even if the sigmas of all pairs of atom types have not been defined. However, the %PAPRECA run will abort as soon as an interference check is performed for pairs of atom types of unknown sigma values.

\subsection sigoptions_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

<hr>

\section insigma init_sigma command

\subsection insigma_syntax Syntax

```bash
init_sigma M N K
```

- M = atom type of the first atom.
- M = atom type of the second atom.
- K = sigma value (in LAMMPS length units).

\subsection insigma_examples Example(s)

```bash
init_sigma 1 3 3.4
init_sigma 2 5 5.5
```

\subsection insigma_description Description

Defines a "sigma" value (as in the Lennard-Jones potential, see [here](https://en.wikipedia.org/wiki/Lennard-Jones_potential) ). Note that, for events requiring collision tests (e.g., deposition events)
the sigma value defines the maximum acceptable distance between two species (i.e., if d &le sigma %PAPRECA assumes that an interference exists and that atom collides) [1], [2].

> **Note 1:**
> The sigma values for all pairs of atom types (diagonal and cross-terms) have to be defined before the start of a %PAPRECA run.  The %PAPRECA run will start normally even if the sigmas of all pairs of atom types have not been defined. However, the %PAPRECA run will abort as soon as an interference check is performed for pairs of atom types of unknown sigma values.

> **Note 2:**
> A \ref sigoptions has to be present in your %PAPRECA input file before you use this command.

> **Note 3:**
> This command can be used as many times as necessary in order to define all the sigmas for the %PAPRECA simulation. When used multiple times, previously defined sigma values will be overwritten.

> **Note 4:**
> %PAPRECA will not check if the defined sigma corresponds to atom types that exist in the simulation (i.e., defined in the LAMMPS input file).

\subsection insigma_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421


<hr>


\section flutypes fluid_atomtypes command

\subsection flutypes_syntax Syntax

```bash
fluid_atomtypes N type1 type2 ... typeN
```

- N = Number of fluid atom types.
- typeX = atom type.

\subsection flutypes_examples Example(s)

```bash
fluid_atomtypes 4 1 2 3 4
fluid_atomtypes 2 4 9
```

\subsection flutypes_description Description

Defines the mobile atom types for the MD stage of the %PAPRECA simulation

> **Note 1:**
> This is a mandatory command. The %PAPRECA simulation will not start unless the fluid atom types have been defined.

> **Note 2:**
> %PAPRECA will not check if the defined mobile atom types are consistent with the mobile types of the MD (LAMMPS) stage.

<hr>


\section frotypes frozen_atomtypes command

\subsection frotypes_syntax Syntax

```bash
frozen_atomtypes N type1 type2 ... typeN
```

- N = Number of frozen atom types.
- typeX = atom type.

\subsection frotypes_examples Example(s)

```bash
frozen_atomtypes 4 1 2 3 4
frozen_atomtypes 2 4 9
```

\subsection frotypes_description Description

Defines the frozen atom types for the MD stage %PAPRECA simulation.

> **Note 1:**
> %PAPRECA will not check if the defined frozen atom types are consistent with the mobile types of the MD (LAMMPS) stage.

> **Note 2:**
> This command is not mandatory. However, you must use this command if certain atom types are defined as frozen in your LAMMPS input file (e.g., if you used "fix freeze_atoms frozen setforce 0.0 0.0 0.0" ). See the [LAMMPS documentation page](https://docs.lammps.org/fix_setforce.html) for more information.

\subsection frotypes_defaults Default

%PAPRECA will start with an empty frozen_atomtypes list if the frozen-atomtypes command is not used.


<hr>



\section time_end time_end command

\subsection time_end_syntax Syntax

```bash
time_end N
```

- N = Positive double number denoting the desired ending time of the simulation.

\subsection time_end_examples Example(s)

```bash
time_end 10
```

\subsection time_end_description Description

Sets a desired ending time for the %PAPRECA run, i.e., the %PAPRECA run will stop when the simulation time is equal to or greater than time_end.
This command does not replace \ref KMC_steps. It simply places a time limit on top of \ref KMC_steps.

\subsection time_end_default Default

time_end = std::numeric_limits< double >::max( ). This is the maximum limit of a double number in C++. Therefore, you will never reach that number.


<hr>



\section height_calculation height_calculation command

\subsection height_calculation_syntax Syntax

```bash
height_calculation style args
```

- style = mass_bins

```bash
mass_bins args = cutoff bin_width
	cutoff = mass percentage cutoff (greater than or equal to 0.0 and smaller than pr equal to 1.0).
	bin_width = width (in length units as defined in the LAMMPS file) of x-y bin.
```


\subsection height_calculation_examples Example(s)

```bash
height_calculation mass_bins 0.8 1.0
```

\subsection height_calculation_description Description

Sets up a height calculation for the %PAPRECA simulation. If this command is used, the height **across the z-coordinate** is calculated before every kMC stage. This command was designed to dynamically capture the height of thin films.

For style mass_bins the height calculation is performed as follows: firstly, the whole simulation box is divided into x-y segments. The x-y segments span from on side to another in the x- and y- directions
and have a thickness (in the z-direction) of bin_width. Secondly, the total mass of each x-y segment is calculated. Finally, a running sum of bin masses is calculated (starting from the lowermost x-y segment) and the film height is defined as the height of the first mass bin for which the running mass sum is greater that or equal to M &times; cutoff (where M is the total mass of the system at the current %PAPRECA step).

\subsection height_calculation_default Default

No height calculations are performed if the user does not include this command in the %PAPRECA input file.


<hr>


\section desorb desorption command

\subsection desorb_syntax Syntax

```bash
desorption height style keyword values
```

- height = deletion height (i.e., delete atoms whose **z-coordinate** is equal to or greater than height). Units are in length units (as defined in the LAMMPS input file).

- style = **gather_all** or **gather_local** or **LAMMPS_region**.

- (OPTIONAL AND NOT AVAILABLE FOR STYLE **LAMMPS_region**) keyword = max
```bash
max values = N
	N = integer number denoting the maximum number of atoms that can be deleted at once.
```

\subsection desorb_examples Example(s)

```bash
desorption 20 gather_local
desorption 40 gather_local max 200
desorption 50 gather_all
desorption 50 gather_all max 100
desorption 100 LAMMPS_region
```

\subsection desorb_description Description

This command does **NOT** set up a predefined desorption kMC event template. This command was created to assist film growth studies and is designed to delete all atoms whose **z-coordinate** is 
equal to or greater than height (as set in the command). 

For styles "gather_all" and "gather_local" results must be identical. However, the choice of style might affect the efficiency of the %PAPRECA run due to different implementations. See the relevant C++ function documentation
for more information about these 2 different approaches: deleteDesorbedAtoms(). Also, for styles "gather_all" and "gather_local" atoms bonded to deleted atoms are also deleted. For example, if the deletion height is set to 30 and an atom is above 30 (LAMMPS length units)
and bonded to another atom whose z-coordinate is 29, then, both atoms will be deleted. Bonded atoms are deleted to prevent "bond atoms missing from proc %d" errors (see [LAMMPS documentation page](https://www.afs.enea.it/software/lammps/doc17/html/Section_errors.html))
within the MD stage of %PAPRECA run. If the "max" keyword is used, then the maximum number of atoms that can be deleted at once becomes "N".

For style "LAMMPS_region" a wrapper (see deleteAtomsInBoxRegion()) around the [region](https://docs.lammps.org/region.html) and [delete_atoms](https://docs.lammps.org/delete_atoms.html) commands are used to delete atoms above the deletion height. Note that, unlike the "gather_all" and "gather_local" styles, the "LAMMPS_region" style will simply delete
all bonded interactions (i.e., bond, angles, dihedrals, and impropers) associated with the deleted atoms, but not the bonded atoms to deleted atoms. Consider This as it may lead to instabilities to do sudden system energy change.

> **Note:**
> Prior to using this command you must set up a height calculation for your %PAPRECA simulation (i.e., a \ref height_calculation must be present in your %PAPRECA above the current command line).

\subsection desorb_default Default

No atom deletions (desorptions) are performed if the user does not include this command in the %PAPRECA input file. Also, for styles "gather_all" and "gather_local", if the \ref desorb is used without the "max" keyword then any number of atoms can be deleted at once.

<hr>

\section minprior minimize_prior command

\subsection minprior_syntax Syntax

```bash
minimize_prior arg values
```

- (REQUIRED) arg = **no** or **yes**

```bash
no values = none
yes values = command
	command = valid LAMMPS minimization command (e.g., minimize).
```

\subsection minprior_examples Example(s)

```bash
minimize_prior yes minimize 1.0e-3 1.0e-5 100 1000
minimize_prior no
```

\subsection minprior_description Description

Defines a minimization command to be executed within the MD stage of the %PAPRECA run and before simulating the MD trajectory. This command might be helpful to relax the system to the closest Potential Energy Surface (PES) valley and avoid instabilities [1], [2] within the MD trajectory (e.g., "bond atoms missing from proc %d" errors (see [LAMMPS documentation page](https://www.afs.enea.it/software/lammps/doc17/html/Section_errors.html)).
Please see the relevant LAMMPS documentation page for the [minimize command](https://docs.lammps.org/minimize.html) for more information.

When the "yes" keyword is utilized, the user has to provide a valid LAMMPS minimization command. Note that, %PAPRECA will not check if the command is valid before the start of the simulation but will probably abort during runtime. Moreover,
note that any valid LAMMPS command can be passed to %PAPRECA with this command.

> **Note:**
> A LAMMPS minimization style has to be defined in the LAMMPS input file (see [min_style](https://docs.lammps.org/min_style.html) for more information regarding minimization styles).

\subsection minprior_default Default

No minimization prior to the MD trajectory is performed if this command is not included in the %PAPRECA input file.

\subsection minprior_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

<hr>

\section minafter minimize_after command

\subsection minafter_syntax Syntax

```bash
minimize_after arg values
```

- (REQUIRED) arg = **no** or **yes**

```bash
no values = none
yes values = command
	command = valid LAMMPS minimization command (e.g., minimize).
```

\subsection minafter_examples Example(s)

```bash
minimize_after yes minimize 1.0e-3 1.0e-5 100 1000
minimize_after no
```

\subsection minafter_description Description

Defines a minimization command to be executed within the MD stage of the %PAPRECA run and right after simulating the MD trajectory. This command might be helpful to relax the system to the closest Potential Energy Surface (PES) valley and avoid instabilities [1], [2] within the MD trajectory (e.g., "bond atoms missing from proc %d" errors (see [LAMMPS documentation page](https://www.afs.enea.it/software/lammps/doc17/html/Section_errors.html)).
Please see the relevant LAMMPS documentation page for the [minimize command](https://docs.lammps.org/minimize.html) for more information.

When the "yes" keyword is utilized, the user has to provide a valid LAMMPS minimization command. Note that, %PAPRECA will not check if the command is valid before the start of the simulation but will probably abort during runtime. Moreover,
note that any valid LAMMPS command can be passed to %PAPRECA with this command.

> **Note:**
> A LAMMPS minimization style has to be defined in the LAMMPS input file (see [min_style](https://docs.lammps.org/min_style.html) for more information regarding minimization styles).

\subsection minafter_default Default

No minimization after the MD trajectory is performed if this command is not included in the %PAPRECA input file.

\subsection minafter_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421


<hr>


\section trajdur trajectory_duration command

\subsection trajdur_syntax Syntax

```bash
trajectory_duration N
```

- N = double number denoting the duration of the MD stage.

\subsection trajdur_examples Example(s)

```bash
trajectory_duration 1000
trajectory_duration 0
```

\subsection trajdur_description Description

Sets the duration of the MD stage of the %PAPRECA run. This represents how many MD timesteps will be simulated by LAMMPS during each MD stage. The trajectory duration can be 0 (means that %PAPRECA will perform a pure off-lattice kMC run).

\subsection trajdur_default Default

trajectory_duration = -1.

<hr>

\section longtrajdur long_trajectory_duration command

\subsection trajdur_syntax Syntax

```bash
long_trajectory_duration N
```

- N = double number denoting the duration of the MD stage. This number has to be greater than the trajectory duration number as set in \ref trajdur.

\subsection longtrajdur_examples Example(s)

```bash
trajectory_duration 1000

trajectory_duration 0
```

\subsection longtrajdur_description Description

Sets the long duration of the MD stage of the %PAPRECA run. This represents how many MD timesteps will be simulated by LAMMPS during each MD stage.

\subsection trajdur_default Default

long_trajectory_duration = -1.

<hr>

\section nvelim nve_lim command

\subsection nvelim_syntax Syntax

```bash
nve_lim N xmax
```

- N = positive integer number denoting the number of LAMMPS MD time steps for nve/lim integration
- xmax = positive double number denoting the maximum distance (in LAMMPS distance units) an atom can move in one timestep

\subsection nvelim_examples Example(s)

```bash
nve_lim 10 0.1

nve_lim 20 0.01
```

\subsection nvelim_description Description

Applies an internal [nve/limit](https://docs.lammps.org/fix_nve_limit.html) integrator to atoms reacted through a create_BondBreak or create_BondForm commands in the kMC stage of Papreca. Bonded atoms of reacted atoms are also collected through a recursive bond-collection operation. The (parent and collected) reacted atoms are integrated via nve/limit for a total of "N" time steps and with a maximum allowed displacement of "xmax". This command can be useful to avoid instabilities due to energy release following reaction events.


\subsection nvelim_default Defaults

No nve/limit integration on any atoms.

<hr>


\section createBreak create_BondBreak command

\subsection createBreak_syntax Syntax

```bash
create_BondBreak atom1_type atom2_type bond_type arg values keyword values
```

- atom1_type = atom type of the first atom.
- atom2_type = atom type of the second atom.
- bond_type = bond type associated with the bond between atom types 1 and 2.

- (REQUIRED) arg = **rate_manual** or **rate_arrhenius**

```bash
rate_manual values = rate
	rate = bond-breaking rate in 1/s (Hz)
rate_arrhenius values = energy frequency temperature
	energy = activation energy in kcal/mol.
	frequency = attempt frequency in 1/s (Hz).
	temperature = temperature in K.
```
- (OPTIONAL) keyword = **catalyzed** or/and **limit**

```bash
catalyzed values = N type1 type2 ... typeN
	N = total number of catalyzing types
	type = atom type that catalyzes the bond-breaking event.
limit values = length_equil length_perc
	length_equil = User-defined assumed equilibrium bond length for that specific bond-breaking event.
	length_perc = Percentage over length_equil ( 0.0 < length_perc < 1.0 )
```

\subsection createBreak_examples Example(s)

```bash
create_BondBreak 4 5 7 rate_arrhenius 13.68 6.99e10 528.15 catalyzed 2 1 8 
create_BondBreak 1 2 5 rate_manual 1.0e13
```

\subsection createBreak_description Description

Create a predefined bond-break template (see PAPRECA::PredefinedReaction and PAPRECA::BondBreak) for the kMC stage of the %PAPRECA run. 
Note that, for the kMC event to function properly, you have to explicitly declare the [bond_style](https://docs.lammps.org/bond_style.html) (typically harmonic) along with the relevant [bond_coeff](https://docs.lammps.org/bond_coeff.html) in your LAMMPS input file.

On every kMC stage, %PAPRECA will attempt to remove bonds with a given probability (based on the chosen rate).
Bonds are removed from the system by calling the PAPRECA::deleteBond() LAMMPS wrapper function and the executeBondBreak() function of the papreca.cpp driver code.

You can provide the bond-breaking rate manually or input the activation energy, attempt frequency, and temperature of that kMC event to obtain the corresponding rate from the Arrhenius equation (see rates_calc.h rates_calc.cpp, and PAPRECA::getRateFromArrhenius() ).

When the catalyzed keyword is used, the bond-breaking event can be selected and executed (within the kMC stage of the %PAPRECA run) only if at least one atom of the specified atom type is present
in the full neighbor list of the parent atom (i.e., the atom on which the bond-breaking was discovered). Note that, %PAPRECA will not check if the provided catalyzing atom type is valid (i.e., exists in your simulation and has been defined in the LAMMPS input file).
Also, note that the bond-breaking probabilities will be influenced by any settings related to the building and updating of the neighbor lists. For example, if your pair_style cutoff is too small, then
fewer "catalyzing" types will be in the neighborhood of the parent atom.

When the limit keyword is used, bond-breaking events are only valid if the current distance between bonded atoms obeys the following inequality: (1-length_perc) * length_equil <= (1+length_perc) * length_equil. Careful, to avoid bond-missing errors it is suggested that the length_equil variable is set to the equilibrium bond length as defined in the LAMMPS input file.

> **Note 1:**
> In the current version each pair of atom types (e.g., type 1 and type 2) and their corresponding bond (e.g., bond type 5 for atom types 1 and 2) are allowed to be associated with only one predefined bond-breaking template.

> **Note 2:**
> %PAPRECA uses a dummy group with id=1 to perform bond deletions (see this LAMMPS wrapper function for more information PAPRECA::deleteBond()), which means that you should **NOT** use this bond id to define some another bonded interaction. Please see the "kmc.lmp" file located in ./Examples/Phosphate Film Growth from TCP on Fe110/ for an example demonstrating how you can define multiple bond types.

<hr>

\section createForm create_BondForm command

\subsection createForm_syntax Syntax

```bash
create_BondForm atom1_type atom2_type bond_dist delete_atoms lone_candidates same_mol arg values
```

- atom1_type = atom type of the first atom.
- atom2_type = atom type of the second atom.
- bond_type = bond type associated with the bond between atom types 1 and 2.
- bond_dist = bonding distance (in length units as defined in LAMMPS). The bonding event is valid if the distance between the two atoms is equal to or smaller than bond_dist.
- delete_atoms = **yes** (delete both atoms after the bond-formation event is executed) or **no**.
- lone_candidates = **yes** (the bond-formation event is valid only if bond candidates have no bonds with other atoms in the system) or **no**.
- same_mol = **yes** (the bond-formation event is valid only if both atoms are associated with different mol IDs) or **no**.

- (REQUIRED) arg = **rate_manual** or **rate_arrhenius**

```bash
rate_manual values = rate
	rate = bond-formation rate in 1/s (Hz)
rate_arrhenius values = energy frequency temperature
	energy = activation energy in kcal/mol.
	frequency = attempt frequency in 1/s (Hz).
	temperature = temperature in K.
```

- (OPTIONAL) keyword = **catalyzed**

```bash
catalyzed values = N type1 type2 ... typeN
	N = total number of catalyzing types
	type = atom type that catalyzes the bond-breaking event.
```

\subsection createForm_examples Example(s)

```bash
create_BondForm 5 4 8 3.7588562 no no no rate_manual 1.0e13
create_BondForm 4 4 12 3.4045958 yes no yes rate_manual 1.0e13
create_BondForm 5 6 9 3.7588562 no no no rate_arrhenius 13.0 1.0e13 500
```

\subsection createForm_description Description

Create a predefined bond-formation template (see PAPRECA::PredefinedBondForm and PAPRECA::BondForm) for the kMC stage of the %PAPRECA run. 
Note that, for the kMC event to function properly, you have to explicitly declare the [bond_style](https://docs.lammps.org/bond_style.html) (typically harmonic) along with the relevant [bond_coeff](https://docs.lammps.org/bond_coeff.html) in your LAMMPS input file.
%PAPRECA does not permit the formation of two bonds between exactly the same atoms. If you wish to create a "double" bond consider defining an additional bond type whose coefficients are representative of a "double" bond.

On every kMC stage, %PAPRECA will attempt to create the bond (as defined in the template) with a given probability (based on the chosen rate).
Bonds are created by calling PAPRECA::formBond() LAMMPS wrapper function and the executeBondForm() function of the papreca.cpp driver code.

You can provide the bond-formation rate manually or input the activation energy, attempt frequency, and temperature of that kMC event to obtain the corresponding rate from the Arrhenius equation (see rates_calc.h rates_calc.cpp, and PAPRECA::getRateFromArrhenius() ).

For more information regarding the catalyzed option please see \ref createBreak.

> **Note:**
> In the current version each pair of atom types (e.g., type 1 and type 2) and their corresponding bond (e.g., bond type 5 for atom types 1 and 2) are allowed to be associated with only one predefined bond-forming template.

<hr>

\section maxbonds species_maxbonds command

\subsection maxbonds_syntax Syntax

```bash
species_maxbonds N M
```

- N = atom type.
- M = maximum number of bonds

\subsection maxbonds_examples Example(s)

```bash
species_maxbonds 2 5 #Limits the maximum number of bonds for atom type 2 to 5.
```

\subsection maxbonds_description Description

Sets a maximum number of bonds for a specific atom type. This command can be useful to prevent "overbonding" in %PAPRECA studies involving bond-formation events.

\subsection maxbonds_default Default

The maximum number of bonds of a species is unlimited if the user does not include this command in the %PAPRECA input file.


<hr>

\section maxbondtypes species_maxbondtypes command

\subsection maxbondtypes_syntax Syntax

```bash
species_maxbondtypes N M K
```

- N = atom type.
- M = bond type.
- K = maximum number of bonds

\subsection maxbondtypes_examples Example(s)

```bash
species_maxbonds 7 3 2 #Limits the maximum number of bonds of type 3 for atom type 7 to 2.
```

\subsection maxbondtypes_description Description

Sets a maximum number of bonds of a specific bond type for a specific atom type. This command can be useful to prevent "overbonding" in %PAPRECA studies involving bond-formation events.

\subsection maxbondtypes_default Default

The maximum number of bonds of a specific bond type for a specific atom species is unlimited if the user does not include this command in the %PAPRECA input file.

<hr>

\section createDepo create_Deposition command

\subsection createDepo_syntax Syntax

```bash
create_Deposition parent_type depo_offset insertion_vel adsorbate_name arg values keyword values
```

- parent_type = atom type of parent atom (i.e., atom on which the event is detected).
- depo_offset = double number denoting the distance (in length units as in LAMMPS) between the parent atom and the center of mass (COM) of the deposited molecule.
- insertion_vel = double number denoting the velocity (in velocity units as in LAMMPS) of the inserted molecule.
- adsobate_name = name of adsorbate (this name has to be INDENTICAL to the molecule name as declared in the LAMMPS input file).

- (REQUIRED) arg = **rate_manual** or **rate_arrhenius** or **rate_hertz**

```bash
rate_manual values = rate
	rate = deposition rate in 1/s (Hz)
rate_arrhenius values = energy frequency temperature
	energy = activation energy in kcal/mol.
	frequency = attempt frequency in 1/s (Hz).
	temperature = temperature in K.
rate_hertz values = pressure area mass temperature
	pressure = partial pressure in Bar.
	area = adsorption site area in Angstroms^2.
	mass = molecule (particle) mass in g/mol.
	temperature = adsorption temperature in K.
```

- (OPTIONAL) keyword = sticking_coeff

```bash
sticking_coeff values = variable or constant
	constant value = double (fixed greater than zero and smaller than or equal to 1.0) number for the sticking coefficient of all deposition events involving the same adsorbate (i.e., defined by the same adsorbate_name).
	variable value = none
```

\subsection createDepo_examples Example(s)

```bash
create_Deposition 4 3.13315 0.0 mmmTCP rate_hertz 0.0001 7.640648 368.37 528.15 sticking_coeff variable
create_Deposition 1 3.503 0.001 particle rate_manual 1.0 sticking_coeff constant 1
```

\subsection createDepo_description Description

Create a predefined deposition template (see PAPRECA::PredefinedDeposition and PAPRECA::Deposition) for the kMC stage of the %PAPRECA run. 
This command will only work properly if you have previously defined a molecule template (i.e., you have prepared a separate molecule file and you have used the molecule command in your LAMMPS input file. See this part of the [LAMMPS documentation](https://docs.lammps.org/molecule.html) for more information.

On every kMC stage, %PAPRECA will attempt to insert a molecule in the system with a given probability (based on the chosen rate).
Molecules are inserted in the system by calling the PAPRECA::insertMolecule() LAMMPS wrapper function and the executeDeposition() function of the papreca.cpp driver code. Note that, the deposition candidate coordinates
coincide with the geometric center of the molecule (as defined in the LAMMPS input file).

For sticking_coeff = constant the user must select a constant sticking coefficient value. The sticking coefficient value will remain unchanged (as set) throughout the simulation.
Conversely, if sticking_coeff = variable, then the sticking coefficient is dynamically calculated by dividing the number of available (collision-free) deposition sites by the total number of sites (occupied and collision-free) [1] , [2].
Note that, if you have defined multiple deposition events with the same adsorbate_name (e.g., mmmTCP), then %PAPRECA will assume that the adsorption sites of all such events are identical and will assign/calculate an identical sticking coefficient.

You can provide the deposition rate manually or input the activation energy, attempt frequency, and temperature of that kMC event to obtain the corresponding rate from the Arrhenius equation (see rates_calc.h rates_calc.cpp, and PAPRECA::getRateFromArrhenius() ).
For predefined deposition events the rate can also be calculated from the kinetic theory of gases (Hertz-Knudsen equation). Please see 
this function PAPRECA::getDepoRateFromHertzKnudsen() and [here](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Surface_Science_(Nix)/02%3A_Adsorption_of_Molecules_on_Surfaces/2.03%3A_Kinetics_of_Adsorption) for a brief theoretical background.

> **Note:**
> In the current version each parent atom (e.g., type 1) is allowed to be associated with only one predefined deposition template. However, different atom types can be associated with the same adsorbate_name (in separate deposition templates). See the phosphates example located in (./Examples/Phosphate Film Growth from TCP on Fe110).

\subsection createDepo_defaults Default

If the sticking_coeff keyword is not used, then the predefined deposition event will be initialized with variable sticking coefficients.

\subsection createDepo_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

<hr>

\section depovecs random_depovecs command

\subsection depovecs_syntax Syntax

```bash
random_depovecs arg values
```

- (REQUIRED) arg = **no** or **yes**

```bash
no values = none.
yes values = none.
```

\subsection depovecs_examples Example(s)

```bash
random_depovecs no
random_depovecs yes
```

\subsection depovevecs_description Description

When random_depovecs are not activated (i.e., random_depovecs = no) the center of mass (COM) of any deposited molecule
is located directly above the parent atom (i.e., the atom on which searches for kMC events are performed) and at a distance of depo_offset (see \ref createDepo).

When random_depovecs are activated (i.e., random_depovecs = yes) the COM of any deposited molecule is located on the surface of the upper hemisphere (i.e., along the +z-direction)
of a sphere centered on the parent atom. The exact position of the COM is determined randomly (by drawing a random number, see getDepoPointCandidateCoords() for more information).

> **Note:**
> In the current version deposition candidates are always placed above the parent atom. This means that you cannot use this version to place a molecule underneath the parent atom. Hence, you cannot get a film growing towards -z.

\subsection depovecs_default Default

random_depovecs = no

<hr>

\section depoheights depoheights command

\subsection depoheights_syntax Syntax

```bash
depoheights height_scan height_reject
```

- height_scan = double number denoting the scan range (in LAMMPS length units) for all deposition events.
- height_reject = double number denoting the rejection height (in LAMMPS length units) for all deposition events.

\subsection depoheights_examples Example(s)

```bash
depoheights 5 20
```

\subsection depoheights_description Description

Set limits for discovering deposition events (see \ref createDepo) in the kMC stage of the %PAPRECA run. Note that, those limitations are applied to all deposition events defined in the simulation.
Specifically, $PAPRECA scans for deposition events between film_height - height_scan and film_height + height scan and rejects deposition candidates above film_height + height_reject.


> **Note:**
> Prior to using this command you must set up a height calculation for your %PAPRECA simulation (i.e., a \ref height_calculation must be present in your %PAPRECA above the current command line).

\subsection depoheights_default Default

Deposition scans are not limited to a distance above/below the current film height if the user does not include this command in the %PAPRECA input file.

<hr>

\section createDiff create_DiffusionHop command

\subsection createDiff_syntax Syntax

```bash
create_DiffusionHop parent_type diff_vel diff_dist displacive diffused_type arg values keyword values
```

- parent_type = atom type of parent atom (i.e., atom on which the event is detected).
- diff_vel = double number denoting the velocity (in velocity units as in LAMMPS) of the diffused atom.
- diff_dist = double number denoting the Eucledian distance (in length units as in LAMMPS) between the parent atom and the vacant site.
- displacive = **yes** (if the parent atom is displaced) or **no** if the parent atom remains intact and a new atom is placed on the vacant site.
- diffused_type = atom type of diffused atom.

- (REQUIRED) arg = **rate_manual** or **rate_arrhenius**

```bash
rate_manual values = rate
	rate = diffusion rate in 1/s (Hz)
rate_arrhenius values = energy frequency temperature
	energy = activation energy in kcal/mol.
	frequency = attempt frequency in 1/s (Hz).
	temperature = temperature in K.
```

- (OPTIONAL) keyword = custom

```bash
custom values = Fe_4PO4neib
	Fe_4PO4neib values = N Ptype
		N = has to be equal to 1
		Ptype = atom type number of the Phosphorus type in the simulation.
```

\subsection createDiff_examples Example(s)

```bash
create_DiffusionHop 1 0.0 4.14468 no 8 rate_arrhenius 11.53 1.4e13 528.15 custom Fe_4PO4neib 1 5
create_DiffusionHop 1 0.0 3 yes 1 rate_manual 1.0e13
```

\subsection createDiff_description Description

Create a predefined diffusion hop template (see PAPRECA::PredefinedDiffusionHop and PAPRECA::Diffusion) for the kMC stage of the %PAPRECA run. 
A "displacive" diffusion template (i.e., displacive = "yes" ), means that the parent atom moves from its initial position to the position of the vacant site.
Conversely, a "non-displacive" diffusion template (i.e., displacive = "no" ), means that the parent atom does not move and that a new atom is inserted on the vacant site. Note that,
in the current version of %PAPRECA, the inserted atom has zero charge and zero velocity.

On every kMC stage, %PAPRECA will attempt to move ("displacive" templates) or insert ("non-displacive" templates) an atom in the system with a given probability (based on the chosen rate).
Atoms are diffused in the system by calling the PAPRECA::diffuseAtom() LAMMPS wrapper function and the executeDiffusion() function of the papreca.cpp driver code.

If the custom template "Fe_4PO4neib" is used (see example above for syntax), then diffusion events require at least four PO4 structures in their neighborhood to be valid.
A PO4 structure consists of Phosphorus atom which is bonded (with an explicit bonded interaction) to four Oxygen atoms.
The "Fe_4PO4neib" custom template performs searches on the neighbor list of the parent atom. Also, note that the diffusion probabilities will be influenced by any settings related to the building and updating of the neighbor lists. For example, if your pair_style cutoff is too small, then
fewer PO4 structures will be in the neighborhood of the parent atom. The "Fe_4PO4neib" custom template was created to cover the needs of a very specific application related
to the formation and growth of thin film from tricresyl phosphate (TCP) molecules on an iron Fe110 surface [1] , [2].

You can provide the diffusion rate manually or input the activation energy, attempt frequency, and temperature of that kMC event to obtain the corresponding rate from the Arrhenius equation (see rates_calc.h rates_calc.cpp, and PAPRECA::getRateFromArrhenius() ).

> **Note:**
> In the current version each parent atom (e.g., type 1) is allowed to be associated with only one predefined diffusion template.

\subsection createDiff_bibliography Bibliography

[1] Ntioudis, S., et al. "PAPRECA: A parallel hybrid off-lattice kinetic Monte Carlo/molecular dynamics simulator", Journal of Open Source Software, 9(98), 6714 (2024). https://doi.org/10.21105/joss.06714

[2] Ntioudis, S., et al. "A hybrid off-lattice kinetic Monte Carlo/molecular dynamics method for amorphous thin film growth", Computational Materials Science, 229, 112421 (2023). https://doi.org/10.1016/j.commatsci.2023.112421

<hr>

\section diffvecs random_diffvecs command

\subsection diffvecs_syntax Syntax

```bash
random_diffvecs arg values
```

- (REQUIRED) arg = **no** or **yes**

```bash
no values = 2D or 3D
	2D = for diffusion sites located only directly above the parent atom.
	3D = for diffusion sites located above as well as below the parent atom.
yes values = 2D or 3D
```

\subsection diffvecs_examples Example(s)

```bash
random_diffvecs no
random_diffvecs no 2D
random_diffvecs yes 3D
```

\subsection diffvecs_description Description

When the "2D" option is used diffusion sites (i.e., locations towards which the atoms migrate) can
only be located above the parent atom (i.e., the atom on which searches for kMC events are performed). This means that atoms will solely diffuse towards the +z-direction.
When the "3D" option is used diffusion sites can be located above as well as below the parent atom.

When random_diffvecs are not activated (i.e., random_diffvecs = no) diffusion sites will be directly above/below (depending on whether the "2D" or "3D" options were used)
and at a distance of diff_dist (see \ref createDiff).

When random_diffvecs are activated (i.e., random_diffvecs  = yes ) diffusion sites will be located on the surface
of a sphere centered on the parent atom. The exact position of the COM is selected randomly (by drawing a random number, see getDiffPointCandidateCoords() for more information).

\subsection diffvecs_default Default

random_diffvecs = no

<hr>

\section createMonodes create_MonoatomicDesorption command

\subsection createMonodes_syntax Syntax

```bash
create_MonoatomicDesorption parent_type arg values
```

- parent_type = atom type of parent atom (i.e., atom on which the event is detected).

- (REQUIRED) arg = **rate_manual** or **rate_arrhenius**

```bash
rate_manual values = rate
	rate = monoatomic desorption rate in 1/s (Hz)
rate_arrhenius values = energy frequency temperature
	energy = activation energy in kcal/mol.
	frequency = attempt frequency in 1/s (Hz).
	temperature = temperature in K.
```

\subsection createMonodes_examples Example(s)

```bash
create_MonoatomicDesorption 2 rate_manual 1.0e13
create_MonoatomicDesorption 5 rate_arrhenius 30 1.0e13 1000
```

\subsection createMonodes_description Description

Create a monoatomic desorption template (see PAPRECA::PredefinedMonoatomicDesorption and PAPRECA::MonoatomicDesorption) for the kMC stage of the %PAPRECA run.
Note that, in the current %PAPRECA version monoatomic desorption can only occur if the candidate atom is lone (i.e., has no explicit bonds with any other system atom).

On every kMC stage, %PAPRECA will attempt to delete a single atom from the system with a given probability (based on the chosen rate).
Atoms are desorbed in the system by calling the PAPRECA::deletAtoms() LAMMPS wrapper function and the executeMonoatomicDesorption() function of the papreca.cpp driver code.

You can provide the monoatomic desorption rate manually or input the activation energy, attempt frequency, and temperature of that kMC event to obtain the corresponding rate from the Arrhenius equation (see rates_calc.h rates_calc.cpp, and PAPRECA::getRateFromArrhenius() ).

> **Note:**
> In the current version each parent atom (e.g., type 1) is allowed to be associated with only one monoatomic desorption template.


<hr>

\section heightvtime export_HeightVtime command

\subsection heightvtime_syntax Syntax

```bash
export_HeightVtime N
```

- N = integer number denoting the print frequency (i.e., %PAPRECA will print on the heightVtime.log file every N %PAPRECA steps).

\subsection heightVtime_examples Example(s)

```bash
export_HeightVtime 10
export_HeightVtime 1000
```

\subsection heightVtime_description Description

Generates a file named "heightVtime.log" in the current directory of the %PAPRECA run. The HeightVtime file comprises 2 columns. 
The first column lists the time (in seconds) of a %PAPRECA step. The second column lists the height of the system (along the z-coordinate) in units consistent with the LAMMPS length units.
%PAPRECA will append to the HeightVtime file every N steps.


\subsection heightVtime_default Default

If this command is not used in your %PAPRECA input file, no "heightVtime.log" file will be generated in your run %PAPRECA run directory.

<hr>

\section coverage export_SurfaceCoverage command

\subsection coverage_syntax Syntax

```bash
export_SurfaceCoverage N
```

- N = integer number denoting the print frequency (i.e., %PAPRECA will print on the surface_coverage.log file every N %PAPRECA steps).

\subsection coverage_examples Example(s)

```bash
export_SurfaceCoverage 10
export_SurfaceCoverage 1000
```

\subsection coverage_description Description

Generates a file named "surface_coverage.log" in the current directory of the %PAPRECA run. The SurfaceCoverage file comprises 2 columns. 
The first column lists the time (in seconds) of a %PAPRECA step. The second column lists the surface coverage, i.e., the number of occupied sites divided by the total number of sites in the system.
Note that, the surface coverage is linked to the sticking coefficient of the deposition event (see \ref createDepo).

%PAPRECA will append to the SurfaceCoverage file every N steps.

\subsection coverage_default Default

If this command is not used in your %PAPRECA input file, no "surface_coverage.log" file will be generated in your run %PAPRECA run directory.

<hr>

\section Edistributions export_ElementalDistributions command

\subsection Edistributions_syntax Syntax

```bash
export_ElementalDistributions N keyword values
```

- N = integer number denoting the print frequency (i.e., %PAPRECA will generate a distribution.log file every N %PAPRECA steps).

- (OPTIONAL) keyword = bin_width values

```bash
bin_width values = width
	width = width (across the z-direction) of the x-y slice used for the domain segmentation.
```

\subsection Edistributions_examples Example(s)

```bash
export_ElementalDistributions 10
export_ElementalDistributions 1000 bin_width 1.0
```

\subsection Edistributions_description Description

Generates a file named distributions.log every N steps, in the current directory of the %PAPRECA run. The first column of the distributions
lists the height (in length units as in LAMMPS). The remaining columns list the number of atoms of each type (in each height-bin).

For more information regarding the x-y bins and the domain segmentation please refer to \ref height_calculation.


\subsection Edistributions_default Default

If this command is not used in your %PAPRECA input file, no "distributions.log" files will be generated in your run %PAPRECA run directory.

<hr>

\section execution export_ExecutionTimes command

\subsection execution_syntax Syntax

```bash
export_ExecutionTimes N
```

- N = integer number denoting the print frequency (i.e., %PAPRECA will print on the execTimes.log file every N %PAPRECA steps).

\subsection execution_examples Example(s)

```bash
export_ExecutionTimes 10
export_ExecutionTimes 300
```

\subsection execution_description Description

Generates a file named "execTimes.log" in the current directory of the %PAPRECA run.
The execTimes file lists the minimum, average, and maximum wall times for the kMC, and MD stages. Note that the minimum, average, and maximum
times will be equal if you choose to run the simulation on a single processor.

%PAPRECA will append to the execTimes.log file every N steps.

> **Note:**
> The total walltime (at the bottom of the execTimes.log file) is the sum of all the average times. Note that, if you choose a print frequency (N) different than 1, the reported total wall time will be smaller than the actual walltime. Nevertheless, the total walltime will be printed (by LAMMPS) in the console at the end of the run.


\subsection execution_default Default

If this command is not used in your %PAPRECA input file, no "execTimes.log" file will be generated in your run %PAPRECA run directory.

<hr>

\section restart restart_freq command

\subsection restart_syntax Syntax

```bash
restart_freq N
```

- N = integer number denoting the restart dump frequency (i.e., %PAPRECA will dump a restart file every N %PAPRECA steps).

\subsection restart_examples Example(s)

```bash
restart_freq 10
restart_freq 10000
```
\subsection restart_description Description

Dump a LAMMPS restart file every N %PAPRECA steps, in the current directory. The LAMMPS restart file can be used to restart a %PAPRECA run. To restart a %PAPRECA run, create a new LAMMPS input 
where you read a dumped LAMMPS restart file (i.e., [read_restart](https://docs.lammps.org/restart.html) and start your run with that LAMMPS input file and the already used %PAPRECA input file.

See the PAPRECA::dumpRestart() function for a brief explanation regarding the reason why the dumping of restart files has to be controlled by %PAPRECA (and not directly from the LAMMPS input file).

\subsection restart_default Default

If this command is not used in your %PAPRECA input file, no restart files will be generated in your run %PAPRECA run directory.


