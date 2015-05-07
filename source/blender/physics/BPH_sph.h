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

typedef struct SPHData {
  ParticleSystem *psys[10];
  ParticleData *pa;
  float mass;
  struct EdgeHash *eh;
  float *gravity;
  float hfac;
  /* Average distance to neighbours (other particles in the support domain),
   * for calculating the Courant number (adaptive time step). */
  int pass;
  float element_size;
  float flow[3];

  /* Integrator callbacks. This allows different SPH implementations. */
  void (*force_cb) (void *sphdata_v, ParticleKey *state, float *force, float *impulse);
  void (*density_cb) (void *rangedata_v, int index, float squared_dist);
} SPHData;

/* General SPH functions */
void BPH_sph_unsplit_particle(ParticleSimulationData *sim, int index, float cfra);
void BPH_sph_split_particle(ParticleSimulationData *sim, int index, float cfra);

/* DDR SPH */
void BPH_sphDDR_step(ParticleSimulationData *sim, float dtime, float cfra);

/* Classical SPH only functions */
void BPH_sphclassical_step(ParticleSimulationData *sim, float dtime, float cfra);

#endif
