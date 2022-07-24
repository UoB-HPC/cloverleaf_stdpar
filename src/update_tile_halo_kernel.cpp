/*
 Crown Copyright 2012 AWE.

 This file is part of CloverLeaf.

 CloverLeaf is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your option)
 any later version.

 CloverLeaf is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 CloverLeaf. If not, see http://www.gnu.org/licenses/.
 */

#include "update_tile_halo_kernel.h"
#include "utils.hpp"

//   @brief Fortran kernel to update the external halo cells in a chunk.
//   @author Wayne Gaudin
//   @details Updates halo cells for the required fields at the required depth
//   for any halo cells that lie on an external boundary. The location and type
//   of data governs how this is carried out. External boundaries are always
//   reflective.


void update_tile_halo_l_kernel(
		int x_min, int x_max, int y_min, int y_max,
		clover::Buffer2D<double> &density0, clover::Buffer2D<double> &energy0,
		clover::Buffer2D<double> &pressure, clover::Buffer2D<double> &viscosity,
		clover::Buffer2D<double> &soundspeed, clover::Buffer2D<double> &density1,
		clover::Buffer2D<double> &energy1, clover::Buffer2D<double> &xvel0,
		clover::Buffer2D<double> &yvel0, clover::Buffer2D<double> &xvel1,
		clover::Buffer2D<double> &yvel1, clover::Buffer2D<double> &vol_flux_x,
		clover::Buffer2D<double> &vol_flux_y,
		clover::Buffer2D<double> &mass_flux_x,
		clover::Buffer2D<double> &mass_flux_y, int left_xmin, int left_xmax,
		int left_ymin, int left_ymax, clover::Buffer2D<double> &left_density0,
		clover::Buffer2D<double> &left_energy0,
		clover::Buffer2D<double> &left_pressure,
		clover::Buffer2D<double> &left_viscosity,
		clover::Buffer2D<double> &left_soundspeed,
		clover::Buffer2D<double> &left_density1,
		clover::Buffer2D<double> &left_energy1,
		clover::Buffer2D<double> &left_xvel0,
		clover::Buffer2D<double> &left_yvel0,
		clover::Buffer2D<double> &left_xvel1,
		clover::Buffer2D<double> &left_yvel1,
		clover::Buffer2D<double> &left_vol_flux_x,
		clover::Buffer2D<double> &left_vol_flux_y,
		clover::Buffer2D<double> &left_mass_flux_x,
		clover::Buffer2D<double> &left_mass_flux_y, const int fields[NUM_FIELDS],
		int depth) {
	// Density 0
	if (fields[field_density0] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &density0, &left_density0](int k) {
			for (int j = 0; j < depth; ++j) {
				density0(x_min - j, k) = left_density0(left_xmax + 1 - j, k);
			}
		}));
	}

	// Density 1
	if (fields[field_density1] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &density1, &left_density1](int k) {
			for (int j = 0; j < depth; ++j) {
				density1(x_min - j, k) = left_density1(left_xmax + 1 - j, k);
			}
		}));
	}

	// Energy 0
	if (fields[field_energy0] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &energy0, &left_energy0](int k) {
			for (int j = 0; j < depth; ++j) {
				energy0(x_min - j, k) = left_energy0(left_xmax + 1 - j, k);
			}
		}));
	}

	// Energy 1
	if (fields[field_energy1] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &energy1, &left_energy1](int k) {
			for (int j = 0; j < depth; ++j) {
				energy1(x_min - j, k) = left_energy1(left_xmax + 1 - j, k);
			}
		}));
	}

	// Pressure
	if (fields[field_pressure] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &pressure, &left_pressure](int k) {
			for (int j = 0; j < depth; ++j) {
				pressure(x_min - j, k) = left_pressure(left_xmax + 1 - j, k);
			}
		}));
	}

	// Viscosity
	if (fields[field_viscosity] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &viscosity, &left_viscosity](int k) {
			for (int j = 0; j < depth; ++j) {
				viscosity(x_min - j, k) = left_viscosity(left_xmax + 1 - j, k);
			}
		}));
	}

	// Soundspeed
	if (fields[field_soundspeed] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &soundspeed, &left_soundspeed](int k) {
			for (int j = 0; j < depth; ++j) {
				soundspeed(x_min - j, k) = left_soundspeed(left_xmax + 1 - j, k);
			}
		}));
	}

	// XVEL 0
	if (fields[field_xvel0] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &xvel0, &left_xvel0](int k) {
			for (int j = 0; j < depth; ++j) {
				xvel0(x_min - j, k) = left_xvel0(left_xmax + 1 - j, k);
			}
		}));
	}

	// XVEL 1
	if (fields[field_xvel1] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &xvel1, &left_xvel1](int k) {
			for (int j = 0; j < depth; ++j) {
				xvel1(x_min - j, k) = left_xvel1(left_xmax + 1 - j, k);
			}
		}));
	}

	// YVEL 0
	if (fields[field_yvel0] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &yvel0, &left_yvel0](int k) {
			for (int j = 0; j < depth; ++j) {
				yvel0(x_min - j, k) = left_yvel0(left_xmax + 1 - j, k);
			}
		}));
	}

	// YVEL 1
	if (fields[field_yvel1] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &yvel1, &left_yvel1](int k) {
			for (int j = 0; j < depth; ++j) {
				yvel1(x_min - j, k) = left_yvel1(left_xmax + 1 - j, k);
			}
		}));
	}

	// VOL_FLUX_X
	if (fields[field_vol_flux_x] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &vol_flux_x, &left_vol_flux_x](int k) {
			for (int j = 0; j < depth; ++j) {
				vol_flux_x(x_min - j, k) = left_vol_flux_x(left_xmax + 1 - j, k);
			}
		}));
	}

	// MASS_FLUX_X
	if (fields[field_mass_flux_x] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &mass_flux_x, &left_mass_flux_x](int k) {
			for (int j = 0; j < depth; ++j) {
				mass_flux_x(x_min - j, k) = left_mass_flux_x(left_xmax + 1 - j, k);
			}
		}));
	}

	// VOL_FLUX_Y
	if (fields[field_vol_flux_y] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &vol_flux_y, &left_vol_flux_y](int k) {
			for (int j = 0; j < depth; ++j) {
				vol_flux_y(x_min - j, k) = left_vol_flux_y(left_xmax + 1 - j, k);
			}
		}));
	}

	// MASS_FLUX_Y
	if (fields[field_mass_flux_y] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &mass_flux_y, &left_mass_flux_y](int k) {
			for (int j = 0; j < depth; ++j) {
				mass_flux_y(x_min - j, k) = left_mass_flux_y(left_xmax + 1 - j, k);
			}
		}));
	}
}

void update_tile_halo_r_kernel(
		int x_min, int x_max, int y_min, int y_max,
		clover::Buffer2D<double> &density0, clover::Buffer2D<double> &energy0,
		clover::Buffer2D<double> &pressure, clover::Buffer2D<double> &viscosity,
		clover::Buffer2D<double> &soundspeed, clover::Buffer2D<double> &density1,
		clover::Buffer2D<double> &energy1, clover::Buffer2D<double> &xvel0,
		clover::Buffer2D<double> &yvel0, clover::Buffer2D<double> &xvel1,
		clover::Buffer2D<double> &yvel1, clover::Buffer2D<double> &vol_flux_x,
		clover::Buffer2D<double> &vol_flux_y,
		clover::Buffer2D<double> &mass_flux_x,
		clover::Buffer2D<double> &mass_flux_y, int right_xmin, int right_xmax,
		int right_ymin, int right_ymax, clover::Buffer2D<double> &right_density0,
		clover::Buffer2D<double> &right_energy0,
		clover::Buffer2D<double> &right_pressure,
		clover::Buffer2D<double> &right_viscosity,
		clover::Buffer2D<double> &right_soundspeed,
		clover::Buffer2D<double> &right_density1,
		clover::Buffer2D<double> &right_energy1,
		clover::Buffer2D<double> &right_xvel0,
		clover::Buffer2D<double> &right_yvel0,
		clover::Buffer2D<double> &right_xvel1,
		clover::Buffer2D<double> &right_yvel1,
		clover::Buffer2D<double> &right_vol_flux_x,
		clover::Buffer2D<double> &right_vol_flux_y,
		clover::Buffer2D<double> &right_mass_flux_x,
		clover::Buffer2D<double> &right_mass_flux_y, const int fields[NUM_FIELDS],
		int depth) {
	// Density 0
	if (fields[field_density0] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &density0, &right_density0](int k) {
			for (int j = 0; j < depth; ++j) {
				density0(x_max + 2 + j, k) = right_density0(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Density 1
	if (fields[field_density1] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &density1, &right_density1](int k) {
			for (int j = 0; j < depth; ++j) {
				density1(x_max + 2 + j, k) = right_density1(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Energy 0
	if (fields[field_energy0] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &energy0, &right_energy0](int k) {
			for (int j = 0; j < depth; ++j) {
				energy0(x_max + 2 + j, k) = right_energy0(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Energy 1
	if (fields[field_energy1] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &energy1, &right_energy1](int k) {
			for (int j = 0; j < depth; ++j) {
				energy1(x_max + 2 + j, k) = right_energy1(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Pressure
	if (fields[field_pressure] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &pressure, &right_pressure](int k) {
			for (int j = 0; j < depth; ++j) {
				pressure(x_max + 2 + j, k) = right_pressure(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Viscosity
	if (fields[field_viscosity] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &viscosity, &right_viscosity](int k) {
			for (int j = 0; j < depth; ++j) {
				viscosity(x_max + 2 + j, k) = right_viscosity(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// Soundspeed
	if (fields[field_soundspeed] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &soundspeed, &right_soundspeed](int k) {
			for (int j = 0; j < depth; ++j) {
				soundspeed(x_max + 2 + j, k) = right_soundspeed(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// XVEL 0
	if (fields[field_xvel0] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &xvel0, &right_xvel0](int k) {
			for (int j = 0; j < depth; ++j) {
				xvel0(x_max + 1 + 2 + j, k) = right_xvel0(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// XVEL 1
	if (fields[field_xvel1] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &xvel1, &right_xvel1](int k) {
			for (int j = 0; j < depth; ++j) {
				xvel1(x_max + 1 + 2 + j, k) = right_xvel1(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// YVEL 0
	if (fields[field_yvel0] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &yvel0, &right_yvel0](int k) {
			for (int j = 0; j < depth; ++j) {
				yvel0(x_max + 1 + 2 + j, k) = right_yvel0(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// YVEL 1
	if (fields[field_yvel1] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &yvel1, &right_yvel1](int k) {
			for (int j = 0; j < depth; ++j) {
				yvel1(x_max + 1 + 2 + j, k) = right_yvel1(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// VOL_FLUX_X
	if (fields[field_vol_flux_x] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &vol_flux_x, &right_vol_flux_x](int k) {
			for (int j = 0; j < depth; ++j) {
				vol_flux_x(x_max + 1 + 2 + j, k) = right_vol_flux_x(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// MASS_FLUX_X
	if (fields[field_mass_flux_x] == 1) {
		// DO k=y_min-depth,y_max+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + depth + 2}, ([=, &mass_flux_x, &right_mass_flux_x](int k) {
			for (int j = 0; j < depth; ++j) {
				mass_flux_x(x_max + 1 + 2 + j, k) = right_mass_flux_x(right_xmin + 1 - 1 + 2 + j, k);
			}
		}));
	}

	// VOL_FLUX_Y
	if (fields[field_vol_flux_y] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &vol_flux_y, &right_vol_flux_y](int k) {
			for (int j = 0; j < depth; ++j) {
				vol_flux_y(x_max + 2 + j, k) = right_vol_flux_y(right_xmin - 1 + 2 + j, k);
			}
		}));
	}

	// MASS_FLUX_Y
	if (fields[field_mass_flux_y] == 1) {
		// DO k=y_min-depth,y_max+1+depth

		clover::par_ranged1(Range1d{y_min - depth + 1, y_max + 1 + depth + 2}, ([=, &mass_flux_y, &right_mass_flux_y](int k) {
			for (int j = 0; j < depth; ++j) {
				mass_flux_y(x_max + 2 + j, k) = right_mass_flux_y(right_xmin - 1 + 2 + j, k);
			}
		}));
	}
}

//  Top and bottom only do xmin -> xmax
//  This is because the corner ghosts will get communicated in the left right
//  communication

void update_tile_halo_t_kernel(
		int x_min, int x_max, int y_min, int y_max,
		clover::Buffer2D<double> &density0, clover::Buffer2D<double> &energy0,
		clover::Buffer2D<double> &pressure, clover::Buffer2D<double> &viscosity,
		clover::Buffer2D<double> &soundspeed, clover::Buffer2D<double> &density1,
		clover::Buffer2D<double> &energy1, clover::Buffer2D<double> &xvel0,
		clover::Buffer2D<double> &yvel0, clover::Buffer2D<double> &xvel1,
		clover::Buffer2D<double> &yvel1, clover::Buffer2D<double> &vol_flux_x,
		clover::Buffer2D<double> &vol_flux_y,
		clover::Buffer2D<double> &mass_flux_x,
		clover::Buffer2D<double> &mass_flux_y, int top_xmin, int top_xmax,
		int top_ymin, int top_ymax, clover::Buffer2D<double> &top_density0,
		clover::Buffer2D<double> &top_energy0,
		clover::Buffer2D<double> &top_pressure,
		clover::Buffer2D<double> &top_viscosity,
		clover::Buffer2D<double> &top_soundspeed,
		clover::Buffer2D<double> &top_density1,
		clover::Buffer2D<double> &top_energy1,
		clover::Buffer2D<double> &top_xvel0, clover::Buffer2D<double> &top_yvel0,
		clover::Buffer2D<double> &top_xvel1, clover::Buffer2D<double> &top_yvel1,
		clover::Buffer2D<double> &top_vol_flux_x,
		clover::Buffer2D<double> &top_vol_flux_y,
		clover::Buffer2D<double> &top_mass_flux_x,
		clover::Buffer2D<double> &top_mass_flux_y, const int fields[NUM_FIELDS],
		int depth) {
	// Density 0
	if (fields[field_density0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &density0, &top_density0](int j) {
				density0(j, y_max + 2 + k) = top_density0(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Density 1
	if (fields[field_density1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &density1, &top_density1](int j) {
				density1(j, y_max + 2 + k) = top_density1(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Energy 0
	if (fields[field_energy0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &energy0, &top_energy0](int j) {
				energy0(j, y_max + 2 + k) = top_energy0(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Energy 1
	if (fields[field_energy1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &energy1, &top_energy1](int j) {
				energy1(j, y_max + 2 + k) = top_energy1(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Pressure
	if (fields[field_pressure] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &pressure, &top_pressure](int j) {
				pressure(j, y_max + 2 + k) = top_pressure(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Viscocity
	if (fields[field_viscosity] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &viscosity, &top_viscosity](int j) {
				viscosity(j, y_max + 2 + k) = top_viscosity(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// Soundspeed
	if (fields[field_soundspeed] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &soundspeed, &top_soundspeed](int j) {
				soundspeed(j, y_max + 2 + k) = top_soundspeed(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// XVEL 0
	if (fields[field_xvel0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &xvel0, &top_xvel0](int j) {
				xvel0(j, y_max + 1 + 2 + k) = top_xvel0(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}

	// XVEL 1
	if (fields[field_xvel1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &xvel1, &top_xvel1](int j) {
				xvel1(j, y_max + 1 + 2 + k) = top_xvel1(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}

	// YVEL 0
	if (fields[field_yvel0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &yvel0, &top_yvel0](int j) {
				yvel0(j, y_max + 1 + 2 + k) = top_yvel0(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}

	// YVEL 1
	if (fields[field_yvel1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &yvel1, &top_yvel1](int j) {
				yvel1(j, y_max + 1 + 2 + k) = top_yvel1(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}

	// VOL_FLUX_X
	if (fields[field_vol_flux_x] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &vol_flux_x, &top_vol_flux_x](int j) {
				vol_flux_x(j, y_max + 2 + k) = top_vol_flux_x(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// MASS_FLUX_X
	if (fields[field_mass_flux_x] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &mass_flux_x, &top_mass_flux_x](int j) {
				mass_flux_x(j, y_max + 2 + k) = top_mass_flux_x(j, top_ymin - 1 + 2 + k);
			}));
		}
	}

	// VOL_FLUX_Y
	if (fields[field_vol_flux_y] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &vol_flux_y, &top_vol_flux_y](int j) {
				vol_flux_y(j, y_max + 1 + 2 + k) = top_vol_flux_y(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}

	// MASS_FLUX_Y
	if (fields[field_mass_flux_y] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &mass_flux_y, &top_mass_flux_y](int j) {
				mass_flux_y(j, y_max + 1 + 2 + k) = top_mass_flux_y(j, top_ymin + 1 - 1 + 2 + k);
			}));
		}
	}
}

void update_tile_halo_b_kernel(
		int x_min, int x_max, int y_min, int y_max,
		clover::Buffer2D<double> &density0, clover::Buffer2D<double> &energy0,
		clover::Buffer2D<double> &pressure, clover::Buffer2D<double> &viscosity,
		clover::Buffer2D<double> &soundspeed, clover::Buffer2D<double> &density1,
		clover::Buffer2D<double> &energy1, clover::Buffer2D<double> &xvel0,
		clover::Buffer2D<double> &yvel0, clover::Buffer2D<double> &xvel1,
		clover::Buffer2D<double> &yvel1, clover::Buffer2D<double> &vol_flux_x,
		clover::Buffer2D<double> &vol_flux_y,
		clover::Buffer2D<double> &mass_flux_x,
		clover::Buffer2D<double> &mass_flux_y, int bottom_xmin, int bottom_xmax,
		int bottom_ymin, int bottom_ymax,
		clover::Buffer2D<double> &bottom_density0,
		clover::Buffer2D<double> &bottom_energy0,
		clover::Buffer2D<double> &bottom_pressure,
		clover::Buffer2D<double> &bottom_viscosity,
		clover::Buffer2D<double> &bottom_soundspeed,
		clover::Buffer2D<double> &bottom_density1,
		clover::Buffer2D<double> &bottom_energy1,
		clover::Buffer2D<double> &bottom_xvel0,
		clover::Buffer2D<double> &bottom_yvel0,
		clover::Buffer2D<double> &bottom_xvel1,
		clover::Buffer2D<double> &bottom_yvel1,
		clover::Buffer2D<double> &bottom_vol_flux_x,
		clover::Buffer2D<double> &bottom_vol_flux_y,
		clover::Buffer2D<double> &bottom_mass_flux_x,
		clover::Buffer2D<double> &bottom_mass_flux_y, const int fields[NUM_FIELDS],
		int depth) {
	// Density 0
	if (fields[field_density0] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &density0, &bottom_density0](int j) {
				density0(j, y_min - k) = bottom_density0(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Density 1
	if (fields[field_density1] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &density1, &bottom_density1](int j) {
				density1(j, y_min - k) = bottom_density1(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Energy 0
	if (fields[field_energy0] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &energy0, &bottom_energy0](int j) {
				energy0(j, y_min - k) = bottom_energy0(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Energy 1
	if (fields[field_energy1] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &energy1, &bottom_energy1](int j) {
				energy1(j, y_min - k) = bottom_energy1(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Pressure
	if (fields[field_pressure] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &pressure, &bottom_pressure](int j) {
				pressure(j, y_min - k) = bottom_pressure(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Viscocity
	if (fields[field_viscosity] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &viscosity, &bottom_viscosity](int j) {
				viscosity(j, y_min - k) = bottom_viscosity(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// Soundspeed
	if (fields[field_soundspeed] == 1) {
		for (int k = 0; k < depth; ++k) {
			//  DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &soundspeed, &bottom_soundspeed](int j) {
				soundspeed(j, y_min - k) = bottom_soundspeed(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// XVEL 0
	if (fields[field_xvel0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &xvel0, &bottom_xvel0](int j) {
				xvel0(j, y_min - k) = bottom_xvel0(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// XVEL 1
	if (fields[field_xvel1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &xvel1, &bottom_xvel1](int j) {
				xvel1(j, y_min - k) = bottom_xvel1(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// YVEL 0
	if (fields[field_yvel0] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &yvel0, &bottom_yvel0](int j) {
				yvel0(j, y_min - k) = bottom_yvel0(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// YVEL 1
	if (fields[field_yvel1] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &yvel1, &bottom_yvel1](int j) {
				yvel1(j, y_min - k) = bottom_yvel1(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// VOL_FLUX_X
	if (fields[field_vol_flux_x] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &vol_flux_x, &bottom_vol_flux_x](int j) {
				vol_flux_x(j, y_min - k) = bottom_vol_flux_x(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// MASS_FLUX_X
	if (fields[field_mass_flux_x] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+1+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + 1 + depth + 2}, ([=, &mass_flux_x, &bottom_mass_flux_x](int j) {
				mass_flux_x(j, y_min - k) = bottom_mass_flux_x(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// VOL_FLUX_Y
	if (fields[field_vol_flux_y] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &vol_flux_y, &bottom_vol_flux_y](int j) {
				vol_flux_y(j, y_min - k) = bottom_vol_flux_y(j, bottom_ymax + 1 - k);
			}));
		}
	}

	// MASS_FLUX_Y
	if (fields[field_mass_flux_y] == 1) {
		for (int k = 0; k < depth; ++k) {
			// DO j=x_min-depth, x_max+depth

			clover::par_ranged1(Range1d{x_min - depth + 1, x_max + depth + 2}, ([=, &mass_flux_y, &bottom_mass_flux_y](int j) {
				mass_flux_y(j, y_min - k) = bottom_mass_flux_y(j, bottom_ymax + 1 - k);
			}));
		}
	}
}


