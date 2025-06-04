/*----------------------------------------------------------------------------------------
PAPRECA hybrid off-lattice kinetic Monte Carlo/Molecular dynamics simulator.
Copyright (C) 2024 Stavros Ntioudis, James P. Ewen, Daniele Dini, and C. Heath Turner

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
----------------------------------------------------------------------------------------*/

/// \file ///
///@brief General PAPRECA header containing all headers of the PAPRECA library.

#ifndef PAPRECA_H
#define PAPRECA_H

//Papreca Library headers
/// \cond
#include "papreca_config.h"
#include "papreca_error.h"

#include "rates_calc.h"
#include "geometry_calc.h"
#include "utilities.h"
#include "lammps_wrappers.h"
#include "mpi_wrappers.h"

#include "bond.h"
#include "debug.h"

#include "event.h"
#include "event_list.h"
#include "event_detect.h"
#include "event_select.h"
#include "event_execute.h"

#include "sim_clock.h"

#include "input_file.h"
#include "export_files.h"

#include "equilibration.h"
/// \endcond

#endif
