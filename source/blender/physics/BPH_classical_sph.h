/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) Blender Foundation
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __BPH_SPH_CLASSICAL_H__
#define __BPH_SPH_CLASSICAL_H__

struct BVHTree;
struct BVHTree_RangeQuery;
struct ParticleSystem;
struct ParticleKey;
struct ParticleSimulationData;
struct ParticleData;
struct SPHData;

#define SPH_NEIGHBORS 512
typedef struct SPHNeighbor {
  ParticleSystem *psys;
  int index;
} SPHNeighbor;

typedef struct SPHRangeData {
  SPHNeighbor neighbors[SPH_NEIGHBORS];
  int tot_neighbors;

  float* data;

  ParticleSystem *npsys;
  ParticleData *pa;

  float h;
  float mass;
  float massfac;
  int use_size;
} SPHRangeData;

/* General SPH functions */
void BPH_sph_evaluate_func(BVHTree*, ParticleSystem**, float*, SPHRangeData*, float, BVHTree_RangeQuery);
void BPH_sph_particle_courant(SPHData*, SPHRangeData*);
//void BPH_sph_integrate(ParticleSimulationData*, ParticleData*, float, SPHData*);

/* Classical SPH only functions */
void BPH_sphclassical_density_accum_cb(void*, int, float);
void BPH_sphclassical_neighbour_accum_cb(void*, int, float);
void BPH_sphclassical_force_cb(void*, ParticleKey*, float*, float*);
void BPH_sphclassical_calc_dens(ParticleData*, float, SPHData*);

#endif
