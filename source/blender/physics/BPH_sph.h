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
 * Classical SPH                                                                                                                                                             
 * Copyright 2011-2012 AutoCRC
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __BPH_SPH_H__
#define __BPH_SPH_H__

struct Object;
struct Scene;
struct BVHTree;
struct BVHTree_RangeQuery;
struct ParticleSettings;
struct ParticleSystem;
struct ParticleKey;
struct ParticleSimulationData;
struct ParticleData;
//struct SPHData;

#define SPH_NEIGHBORS 512
typedef struct SPHNeighbor {
  ParticleSystem *psys;
  int index;
} SPHNeighbor;

typedef struct SPHParams {
	float density;
	float near_density;
	float pressure;
	float near_pressure;
} SPHParams;

typedef struct SPHRangeData {
  SPHNeighbor neighbors[SPH_NEIGHBORS];
  int tot_neighbors;

  SPHParams params;

  ParticleSystem *npsys;
  ParticleData *pa;

  float h;
  float mass;
  float massfac;
  int use_size;
} SPHRangeData;

typedef struct SPHData {
  ParticleSystem *psys[10];
  ParticleData *pa;
  float mass;
  struct EdgeHash *eh;
  float *gravity;
  float hfac;

  /* Some parameters are needed for the equation of state. These are set in the
   * init callback. */
  float rest_density;
  float stiffness;
  float stiffness_near_fac;

  /* Average distance to neighbours (other particles in the support domain),
   * for calculating the Courant number (adaptive time step). */
  int pass;
  float element_size;
  float flow[3];

  /* Integrator callbacks. This allows different SPH implementations. */
  void (*init) (struct SPHData *sphdata);
  void (*force_cb) (void *sphdata_v, ParticleKey *state, float *force, float *impulse);
  void (*density_cb) (void *rangedata_v, int index, float squared_dist);
  void (*equation_of_state) (struct SPHData *sphdata, SPHParams *params);
} SPHData;

/* General SPH functions */
void BPH_sph_unsplit_particle(struct ParticleSimulationData *sim, float cfra);
void BPH_sph_split_particle(struct ParticleSimulationData *sim, int index, float cfra);
void BPH_sph_planar_split(struct ParticleSimulationData *sim, int index, float cfra);

/* DDR SPH */
void BPH_sphDDR_step(struct ParticleSimulationData *sim, float dtime, float cfra);

/* Classical SPH only functions */
void BPH_sphclassical_step(struct ParticleSimulationData *sim, float dtime, float cfra);

/* Adaptive resolution */
void BPH_sph_refiners_init(struct ListBase **refiners, struct ParticleSystem *psys);
void BPH_sph_refiners_end(struct ListBase **refiners);
void sphclassical_refiners_add(struct ListBase **refiners);
void BPH_sph_refiners_remove(struct ListBase **refiners);

void psys_sph_init(struct ParticleSimulationData *sim, SPHData *sphdata);
void psys_sph_finalise(SPHData *sphdata);
void psys_sph_sample(struct BVHTree *tree, SPHData *sphdata, float co[3], SPHParams *params);

#endif
