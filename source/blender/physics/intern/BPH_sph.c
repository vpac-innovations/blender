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

/** \file blender/physics/intern/BPH_classical_sph.c
 *  ingroup bph
 */

#include "MEM_guardedalloc.h"

#include "DNA_particle_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BLI_math.h"
#include "BLI_threads.h"
#include "BLI_kdopbvh.h"
#include "BLI_edgehash.h"

#include "BKE_particle.h"

#include "BPH_sph.h"

static ThreadRWMutex psys_bvhtree_rwlock = BLI_RWLOCK_INITIALIZER;

#define PSYS_FLUID_SPRINGS_INITIAL_SIZE 256

/* Calculate the speed of the particle relative to the local scale of the
 * simulation. This should be called once per particle during a simulation
 * step, after the velocity has been updated. element_size defines the scale of                                                                        
 * the simulation, and is typically the distance to neighboring particles. */
static void BPH_sph_update_courant_num(ParticleSimulationData *sim, ParticleData *pa,
				       float dtime, SPHData *sphdata)
{
  float relative_vel[3];
  float speed;

  sub_v3_v3v3(relative_vel, pa->prev_state.vel, sphdata->flow);
  speed = len_v3(relative_vel);
  if (sim->courant_num < speed * dtime / sphdata->element_size)
    sim->courant_num = speed * dtime / sphdata->element_size;
}

static ParticleSpring *BPH_sph_spring_add(ParticleSystem *psys, ParticleSpring *spring)
{
  /* Are more refs required? */
  if (psys->alloc_fluidsprings == 0 || psys->fluid_springs == NULL) {
    psys->alloc_fluidsprings = PSYS_FLUID_SPRINGS_INITIAL_SIZE;
    psys->fluid_springs = (ParticleSpring*)MEM_callocN(psys->alloc_fluidsprings * sizeof(ParticleSpring), "Particle Fluid Springs");
  }
  else if (psys->tot_fluidsprings == psys->alloc_fluidsprings) {
    /* Double the number of refs allocated */
    psys->alloc_fluidsprings *= 2;
    psys->fluid_springs = (ParticleSpring*)MEM_reallocN(psys->fluid_springs, psys->alloc_fluidsprings * sizeof(ParticleSpring));
  }

  memcpy(psys->fluid_springs + psys->tot_fluidsprings, spring, sizeof(ParticleSpring));
  psys->tot_fluidsprings++;

  return psys->fluid_springs + psys->tot_fluidsprings - 1;
}


static EdgeHash *BPH_sph_springhash_build(ParticleSystem *psys)
{
  EdgeHash *springhash = NULL;
  ParticleSpring *spring;
  int i = 0;

  springhash = BLI_edgehash_new_ex(__func__, psys->tot_fluidsprings);

  for (i=0, spring=psys->fluid_springs; i<psys->tot_fluidsprings; i++, spring++)
    BLI_edgehash_insert(springhash, spring->particle_index[0], spring->particle_index[1], SET_INT_IN_POINTER(i+1));

  return springhash;
}

static void BPH_sph_spring_delete(ParticleSystem *psys, int j)
{
	if (j != psys->tot_fluidsprings - 1)
		psys->fluid_springs[j] = psys->fluid_springs[psys->tot_fluidsprings - 1];

	psys->tot_fluidsprings--;

	if (psys->tot_fluidsprings < psys->alloc_fluidsprings/2 && psys->alloc_fluidsprings > PSYS_FLUID_SPRINGS_INITIAL_SIZE) {
		psys->alloc_fluidsprings /= 2;
		psys->fluid_springs = (ParticleSpring*)MEM_reallocN(psys->fluid_springs,  psys->alloc_fluidsprings * sizeof(ParticleSpring));
	}
}

static void BPH_sph_springs_modify(ParticleSystem *psys, float dtime)
{
	SPHFluidSettings *fluid = psys->part->fluid;
	ParticleData *pa1, *pa2;
	ParticleSpring *spring = psys->fluid_springs;
	
	float h, d, Rij[3], rij, Lij;
	int i;

	float yield_ratio = fluid->yield_ratio;
	float plasticity = fluid->plasticity_constant;
	/* scale things according to dtime */
	float timefix = 25.f * dtime;

	if ((fluid->flag & SPH_VISCOELASTIC_SPRINGS)==0 || fluid->spring_k == 0.f)
		return;

	/* Loop through the springs */
	for (i=0; i<psys->tot_fluidsprings; i++, spring++) {
		pa1 = psys->particles + spring->particle_index[0];
		pa2 = psys->particles + spring->particle_index[1];

		sub_v3_v3v3(Rij, pa2->prev_state.co, pa1->prev_state.co);
		rij = normalize_v3(Rij);

		/* adjust rest length */
		Lij = spring->rest_length;
		d = yield_ratio * timefix * Lij;

		if (rij > Lij + d) // Stretch
			spring->rest_length += plasticity * (rij - Lij - d) * timefix;
		else if (rij < Lij - d) // Compress
			spring->rest_length -= plasticity * (Lij - d - rij) * timefix;

		h = 4.f*pa1->size;

		if (spring->rest_length > h)
			spring->delete_flag = 1;
	}

	/* Loop through springs backwaqrds - for efficient delete function */
	for (i=psys->tot_fluidsprings-1; i >= 0; i--) {
		if (psys->fluid_springs[i].delete_flag)
			BPH_sph_spring_delete(psys, i);
	}
}

static void BPH_sph_evaluate_func(BVHTree *tree, ParticleSystem **psys, float co[3], SPHRangeData *pfr, float interaction_radius, BVHTree_RangeQuery callback)
{
  int i;

  pfr->tot_neighbors = 0;

  for (i=0; i < 10 && psys[i]; i++) {
    pfr->npsys    = psys[i];
    pfr->massfac  = psys[i]->part->mass / pfr->mass;
    pfr->use_size = psys[i]->part->flag & PART_SIZEMASS;

    if (tree) {
      BLI_bvhtree_range_query(tree, co, interaction_radius, callback, pfr);
      break;
    }
    else {
      BLI_rw_mutex_lock(&psys_bvhtree_rwlock, THREAD_LOCK_READ);
      BLI_bvhtree_range_query(psys[i]->bvhtree, co, interaction_radius, callback, pfr);
      BLI_rw_mutex_unlock(&psys_bvhtree_rwlock);
    }
  }
}

static void BPH_sph_particle_courant(SPHData *sphdata, SPHRangeData *pfr)
{
  ParticleData *pa, *npa;
  int i;
  float flow[3], offset[3], dist;

  zero_v3(flow);

  dist = 0.0f;
  if (pfr->tot_neighbors > 0) {
    pa = pfr->pa;
    for (i=0; i < pfr->tot_neighbors; i++) {
      npa = pfr->neighbors[i].psys->particles + pfr->neighbors[i].index;
      sub_v3_v3v3(offset, pa->prev_state.co, npa->prev_state.co);
      dist += len_v3(offset);
      add_v3_v3(flow, npa->prev_state.vel);
    }
    dist += sphdata->psys[0]->part->fluid->radius; // TODO: remove this? - z0r            
    sphdata->element_size = dist / pfr->tot_neighbors;
    mul_v3_v3fl(sphdata->flow, flow, 1.0f / pfr->tot_neighbors);
  }
  else {
    sphdata->element_size = FLT_MAX;
    copy_v3_v3(sphdata->flow, flow);
  }
}

static void BPH_sph_density_accum_cb(void *userdata, int index, float squared_dist)
{
  SPHRangeData *pfr = (SPHRangeData *)userdata;
  ParticleData *npa = pfr->npsys->particles + index;
  float q;
  float dist;

  if (npa == pfr->pa || squared_dist < FLT_EPSILON)
    return;

  /* Ugh! One particle has too many neighbors! If some aren't taken into                                                                                                                               
   * account, the forces will be biased by the tree search order. This                                                                                                                                 
   * effectively adds enery to the system, and results in a churning motion.                                                                                                                           
   * But, we have to stop somewhere, and it's not the end of the world.                                                                                                                                
   *  - jahka and z0r                                                                                                                                                                                  
   */
  if (pfr->tot_neighbors >= SPH_NEIGHBORS)
    return;

  pfr->neighbors[pfr->tot_neighbors].index = index;
  pfr->neighbors[pfr->tot_neighbors].psys = pfr->npsys;
  pfr->tot_neighbors++;

  dist = sqrtf(squared_dist);
  q = (1.f - dist/pfr->h) * pfr->massfac;

  if (pfr->use_size)
    q *= npa->size;

  pfr->data[0] += q*q;
  pfr->data[1] += q*q*q;
}

static void BPH_sph_force_cb(void *sphdata_v, ParticleKey *state, float *force, float *UNUSED(impulse))
{
	SPHData *sphdata = (SPHData *)sphdata_v;
	ParticleSystem **psys = sphdata->psys;
	ParticleData *pa = sphdata->pa;
	SPHFluidSettings *fluid = psys[0]->part->fluid;
	ParticleSpring *spring = NULL;
	SPHRangeData pfr;
	SPHNeighbor *pfn;
	float *gravity = sphdata->gravity;
	EdgeHash *springhash = sphdata->eh;

	float q, u, rij, dv[3];
	float pressure, near_pressure;

	float visc = fluid->viscosity_omega;
	float stiff_visc = fluid->viscosity_beta * (fluid->flag & SPH_FAC_VISCOSITY ? fluid->viscosity_omega : 1.f);

	float inv_mass = 1.0f / sphdata->mass;
	float spring_constant = fluid->spring_k;

	/* 4.0 seems to be a pretty good value */
	float interaction_radius = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * pa->size : 1.0f);
	float h = interaction_radius * sphdata->hfac;
	float rest_density = fluid->rest_density * (fluid->flag & SPH_FAC_DENSITY ? 4.77f : 1.f); /* 4.77 is an experimentally determined density factor */
	float rest_length = fluid->rest_length * (fluid->flag & SPH_FAC_REST_LENGTH ? 2.588f * pa->size : 1.f);

	float stiffness = fluid->stiffness_k;
	float stiffness_near_fac = fluid->stiffness_knear * (fluid->flag & SPH_FAC_REPULSION ? fluid->stiffness_k : 1.f);

	ParticleData *npa;
	float vec[3];
	float vel[3];
	float co[3];
	float data[2];
	float density, near_density;

	int i, spring_index, index = pa - psys[0]->particles;

	data[0] = data[1] = 0;
	pfr.data = data;
	pfr.h = h;
	pfr.pa = pa;
	pfr.mass = sphdata->mass;

	BPH_sph_evaluate_func( NULL, psys, state->co, &pfr, interaction_radius, BPH_sph_density_accum_cb);

	density = data[0];
	near_density = data[1];

	pressure =  stiffness * (density - rest_density);
	near_pressure = stiffness_near_fac * near_density;

	pfn = pfr.neighbors;
	for (i=0; i<pfr.tot_neighbors; i++, pfn++) {
		npa = pfn->psys->particles + pfn->index;

		madd_v3_v3v3fl(co, npa->prev_state.co, npa->prev_state.vel, state->time);

		sub_v3_v3v3(vec, co, state->co);
		rij = normalize_v3(vec);

		q = (1.f - rij/h) * pfn->psys->part->mass * inv_mass;

		if (pfn->psys->part->flag & PART_SIZEMASS)
			q *= npa->size;

		copy_v3_v3(vel, npa->prev_state.vel);

		/* Double Density Relaxation */
		madd_v3_v3fl(force, vec, -(pressure + near_pressure*q)*q);

		/* Viscosity */
		if (visc > 0.f  || stiff_visc > 0.f) {
			sub_v3_v3v3(dv, vel, state->vel);
			u = dot_v3v3(vec, dv);

			if (u < 0.f && visc > 0.f)
				madd_v3_v3fl(force, vec, 0.5f * q * visc * u );

			if (u > 0.f && stiff_visc > 0.f)
				madd_v3_v3fl(force, vec, 0.5f * q * stiff_visc * u );
		}

		if (spring_constant > 0.f) {
			/* Viscoelastic spring force */
			if (pfn->psys == psys[0] && fluid->flag & SPH_VISCOELASTIC_SPRINGS && springhash) {
				/* BLI_edgehash_lookup appears to be thread-safe. - z0r */
				spring_index = GET_INT_FROM_POINTER(BLI_edgehash_lookup(springhash, index, pfn->index));

				if (spring_index) {
					spring = psys[0]->fluid_springs + spring_index - 1;

					madd_v3_v3fl(force, vec, -10.f * spring_constant * (1.f - rij/h) * (spring->rest_length - rij));
				}
				else if (fluid->spring_frames == 0 || (pa->prev_state.time-pa->time) <= fluid->spring_frames) {
					ParticleSpring temp_spring;
					temp_spring.particle_index[0] = index;
					temp_spring.particle_index[1] = pfn->index;
					temp_spring.rest_length = (fluid->flag & SPH_CURRENT_REST_LENGTH) ? rij : rest_length;
					temp_spring.delete_flag = 0;

					/* sph_spring_add is not thread-safe. - z0r */
#pragma omp critical
					BPH_sph_spring_add(psys[0], &temp_spring);
				}
			}
			else {/* PART_SPRING_HOOKES - Hooke's spring force */
				madd_v3_v3fl(force, vec, -10.f * spring_constant * (1.f - rij/h) * (rest_length - rij));
			}
		}
	}
	
	/* Artificial buoyancy force in negative gravity direction  */
	if (fluid->buoyancy > 0.f && gravity)
		madd_v3_v3fl(force, gravity, fluid->buoyancy * (density-rest_density));

	if (sphdata->pass == 0 && psys[0]->part->time_flag & PART_TIME_AUTOSF)
		BPH_sph_particle_courant(sphdata, &pfr);
	sphdata->pass++;
}

static void BPH_sphclassical_density_accum_cb(void *userdata, int index, float UNUSED(squared_dist))
{
  SPHRangeData *pfr = (SPHRangeData *)userdata;
  ParticleData *npa = pfr->npsys->particles + index;
  float q;
  float qfac = 21.0f / (256.f * (float)M_PI);
  float rij, rij_h;
  float vec[3];

  /* Exclude particles that are more than 2h away. Can't use squared_dist here                  
   * because it is not accurate enough. Use current state, i.e. the output of                   
   * basic_integrate() - z0r */
  sub_v3_v3v3(vec, npa->state.co, pfr->pa->state.co);
  rij = len_v3(vec);
  rij_h = rij / pfr->h;
  if (rij_h > 2.0f)
  return;

  /* Smoothing factor. Utilise the Wendland kernel. gnuplot:                                    
   *     q1(x) = (2.0 - x)**4 * ( 1.0 + 2.0 * x)                                                
   *     plot [0:2] q1(x) */

  
  q  = qfac / pow3f(pfr->h) * pow4f(2.0f - rij_h) * ( 1.0f + 2.0f * rij_h);
  q *= pfr->npsys->part->mass;
  
  if (pfr->use_size)
    q *= pfr->pa->size;

  pfr->data[0] += q;
  pfr->data[1] += q / npa->sphdensity;
}

static void BPH_sphclassical_neighbour_accum_cb(void *userdata, int index, float UNUSED(squared_dist))
{
  SPHRangeData *pfr = (SPHRangeData *)userdata;
  ParticleData *npa = pfr->npsys->particles + index;
  float rij, rij_h;
  float vec[3];

  if (pfr->tot_neighbors >= SPH_NEIGHBORS)
    return;

  /* Exclude particles that are more than 2h away. Can't use squared_dist here                  
   * because it is not accurate enough. Use current state, i.e. the output of                   
   * basic_integrate() - z0r */
  sub_v3_v3v3(vec, npa->state.co, pfr->pa->state.co);
  rij = len_v3(vec);
  rij_h = rij / pfr->h;
  if (rij_h > 2.0f)
    return;

  pfr->neighbors[pfr->tot_neighbors].index = index;
  pfr->neighbors[pfr->tot_neighbors].psys = pfr->npsys;
  pfr->tot_neighbors++;
}

static void BPH_sphclassical_force_cb(void *sphdata_v, ParticleKey *state, float *force, float *UNUSED(impulse))
{
	SPHData *sphdata = (SPHData *)sphdata_v;
	ParticleSystem **psys = sphdata->psys;
	ParticleData *pa = sphdata->pa;
	SPHFluidSettings *fluid = psys[0]->part->fluid;
	SPHRangeData pfr;
	SPHNeighbor *pfn;
	float *gravity = sphdata->gravity;

	float dq, u, rij, dv[3];
	float pressure, npressure;

	float visc = fluid->viscosity_omega;

	float interaction_radius;
	float h, hinv;
	/* 4.77 is an experimentally determined density factor */
	float rest_density = fluid->rest_density * (fluid->flag & SPH_FAC_DENSITY ? 4.77f : 1.0f);

	// Use speed of sound squared
	float stiffness = pow2f(fluid->stiffness_k);

	ParticleData *npa;
	float vec[3];
	float co[3];
	float pressureTerm;

	int i;

	float qfac2 = 42.0f / (256.0f * (float)M_PI);
	float rij_h;

	/* 4.0 here is to be consistent with previous formulation/interface */
	interaction_radius = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * pa->size : 1.0f);
	h = interaction_radius * sphdata->hfac;
	hinv = 1.0f / h;

	pfr.h = h;
	pfr.pa = pa;

	BPH_sph_evaluate_func(NULL, psys, state->co, &pfr, interaction_radius, BPH_sphclassical_neighbour_accum_cb);
	pressure =  stiffness * (pow7f(pa->sphdensity / rest_density) - 1.0f);

	/* multiply by mass so that we return a force, not accel */
	qfac2 *= sphdata->mass / pow3f(pfr.h);

	pfn = pfr.neighbors;
	for (i = 0; i < pfr.tot_neighbors; i++, pfn++) {
		npa = pfn->psys->particles + pfn->index;
		if (npa == pa) {
			/* we do not contribute to ourselves */
			continue;
		}

		/* Find vector to neighbor. Exclude particles that are more than 2h
		 * away. Can't use current state here because it may have changed on
		 * another thread - so do own mini integration. Unlike basic_integrate,
		 * SPH integration depends on neighboring particles. - z0r */
		madd_v3_v3v3fl(co, npa->prev_state.co, npa->prev_state.vel, state->time);
		sub_v3_v3v3(vec, co, state->co);
		rij = normalize_v3(vec);
		rij_h = rij / pfr.h;
		if (rij_h > 2.0f)
			continue;

		npressure = stiffness * (pow7f(npa->sphdensity / rest_density) - 1.0f);

		/* First derivative of smoothing factor. Utilise the Wendland kernel.
		 * gnuplot:
		 *     q2(x) = 2.0 * (2.0 - x)**4 - 4.0 * (2.0 - x)**3 * (1.0 + 2.0 * x)
		 *     plot [0:2] q2(x)
		 * Particles > 2h away are excluded above. */
		dq = qfac2 * (pow4f(2.0f - rij_h) - 2.0f * pow3f(2.0f - rij_h) * (1.0f + 2.0f * rij_h)  );
		dq *= sphdata->mass;

		if (pfn->psys->part->flag & PART_SIZEMASS)
			dq *= npa->size;

		pressureTerm = pressure / pow2f(pa->sphdensity) + npressure / pow2f(npa->sphdensity);

		/* Note that 'minus' is removed, because vec = vecBA, not vecAB.
		 * This applies to the viscosity calculation below, too. */
		madd_v3_v3fl(force, vec, pressureTerm * dq);

		/* Viscosity */
		if (visc > 0.0f) {
			sub_v3_v3v3(dv, npa->prev_state.vel, pa->prev_state.vel);
			u = dot_v3v3(vec, dv);
			/* Apply parameters */
			u *= -dq * hinv * visc / (0.5f * npa->sphdensity + 0.5f * pa->sphdensity);
			madd_v3_v3fl(force, vec, u);
		}
	}

	/* Artificial buoyancy force in negative gravity direction  */
	if (fluid->buoyancy > 0.f && gravity)
		madd_v3_v3fl(force, gravity, fluid->buoyancy * (pa->sphdensity - rest_density));

	if (sphdata->pass == 0 && psys[0]->part->time_flag & PART_TIME_AUTOSF)
		BPH_sph_particle_courant(sphdata, &pfr);
	sphdata->pass++;
}

void BPH_psys_sph_init(ParticleSimulationData *sim, SPHData *sphdata)
{
	ParticleTarget *pt;
	int i;

	// Add other coupled particle systems.
	sphdata->psys[0] = sim->psys;
	for (i=1, pt=sim->psys->targets.first; i<10; i++, pt=(pt?pt->next:NULL))
		sphdata->psys[i] = pt ? psys_get_target_system(sim->ob, pt) : NULL;

	if (psys_uses_gravity(sim))
		sphdata->gravity = sim->scene->physics_settings.gravity;
	else
		sphdata->gravity = NULL;
	sphdata->eh = BPH_sph_springhash_build(sim->psys);

	// These per-particle values should be overridden later, but just for
	// completeness we give them default values now.
	sphdata->pa = NULL;
	sphdata->mass = 1.0f;

	if (sim->psys->part->fluid->solver == SPH_SOLVER_DDR) {
		sphdata->force_cb = BPH_sph_force_cb;
		sphdata->density_cb = BPH_sph_density_accum_cb;
		sphdata->hfac = 1.0f;
	}
	else {
		/* SPH_SOLVER_CLASSICAL */
		sphdata->force_cb = BPH_sphclassical_force_cb;
		sphdata->density_cb = BPH_sphclassical_density_accum_cb;
		sphdata->hfac = 0.5f;
	}

}

void BPH_psys_sph_finalise(SPHData *sphdata)
{
	if (sphdata->eh) {
		BLI_edgehash_free(sphdata->eh, NULL);
		sphdata->eh = NULL;
	}
}

/* Sample the density field at a point in space. */
void BPH_psys_sph_density(BVHTree *tree, SPHData *sphdata, float co[3], float vars[2])
{
  ParticleSystem **psys = sphdata->psys;
  SPHFluidSettings *fluid = psys[0]->part->fluid;
  /* 4.0 seems to be a pretty good value */
  float interaction_radius  = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * psys[0]->part->size : 1.0f);
  SPHRangeData pfr;
  float density[2];

  density[0] = density[1] = 0.0f;
  pfr.data = density;
  pfr.h = interaction_radius * sphdata->hfac;
  pfr.mass = sphdata->mass;

  BPH_sph_evaluate_func(tree, psys, co, &pfr, interaction_radius, sphdata->density_cb);

  vars[0] = pfr.data[0];
  vars[1] = pfr.data[1];
}

static void BPH_sphclassical_calc_dens(ParticleData *pa, float UNUSED(dfra), SPHData *sphdata)
{
  ParticleSystem **psys = sphdata->psys;
  SPHFluidSettings *fluid = psys[0]->part->fluid;
  /* 4.0 seems to be a pretty good value */
  float interaction_radius  = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * psys[0]->part->size : 1.0f);
  SPHRangeData pfr;
  float data[2];

  data[0] = 0;
  data[1] = 0;
  pfr.data = data;
  pfr.h = interaction_radius * sphdata->hfac;
  pfr.pa = pa;
  pfr.mass = sphdata->mass;

  BPH_sph_evaluate_func( NULL, psys, pa->state.co, &pfr, interaction_radius, BPH_sphclassical_density_accum_cb);
  pa->sphdensity = MIN2(MAX2(data[0], fluid->rest_density * 0.9f), fluid->rest_density * 1.1f);

}

static void BPH_sph_integrate(ParticleSimulationData *sim, 
			      ParticleData           *pa, 
			      float                   dfra, 
			      SPHData                *sphdata)
{
  ParticleSettings *part = sim->psys->part;
  // float timestep = psys_get_timestep(sim); // UNUSED                                                                                                                                                
  float pa_mass = part->mass * (part->flag & PART_SIZEMASS ? pa->size : 1.f);
  float dtime = dfra*psys_get_timestep(sim);
  // int steps = 1; // UNUSED                                                                                                                                                                          
  float effector_acceleration[3];

  sphdata->pa = pa;
  sphdata->mass = pa_mass;
  sphdata->pass = 0;
  //sphdata.element_size and sphdata.flow are set in the callback.                                                                                                     
  /* restore previous state and treat gravity & effectors as external acceleration*/
  sub_v3_v3v3(effector_acceleration, pa->state.vel, pa->prev_state.vel);
  mul_v3_fl(effector_acceleration, 1.f/dtime);

  copy_particle_key(&pa->state, &pa->prev_state, 0);
  
  integrate_particle(part, pa, dtime, effector_acceleration, sphdata->force_cb, sphdata);
}

void BPH_sphDDR_step(ParticleSimulationData *sim, ParticleData *pa, float cfra)
{
  SPHData sphdata;
  ParticleSettings *part = sim->psys->part;
  ParticleSystem *psys = sim->psys;
  float timestep, dtime;
  int p;

  timestep = psys_get_timestep(sim);
  dtime = pa->state.time*timestep;
  BPH_psys_sph_init(sim, &sphdata);

  /* Apply SPH forces using double-density relaxation algorithm                                                                                                                
   * (Clavat et. al.) */
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    /* do global forces & effectors */
    basic_integrate(sim, p, pa->state.time, cfra);

    /* actual fluids calculations */
    BPH_sph_integrate(sim, pa, pa->state.time, &sphdata);

    if (sim->colliders)
      collision_check(sim, p, pa->state.time, cfra);

    /* SPH particles are not physical particles, just interpolation                                                                                                      
     * particles,  thus rotation has not a direct sense for them */
    basic_rotate(part, pa, pa->state.time, timestep);

#pragma omp critical
    if (part->time_flag & PART_TIME_AUTOSF)
      BPH_sph_update_courant_num(sim, pa, dtime, &sphdata);
  }

  BPH_sph_springs_modify(psys, timestep);
  BPH_psys_sph_finalise(&sphdata);
}

void BPH_sphclassical_step(ParticleSimulationData *sim, ParticleData *pa, float cfra)
{
  SPHData sphdata;
  ParticleSettings *part = sim->psys->part;
  ParticleSystem *psys = sim->psys;
  float timestep, dtime;
  int p;

  timestep = psys_get_timestep(sim);
  dtime = pa->state.time*timestep;
  BPH_psys_sph_init(sim, &sphdata);

#pragma omp parallel for private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    basic_integrate(sim, p, pa->state.time, cfra);
  }
  
  /* Calculate summation density */
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    BPH_sphclassical_calc_dens(pa, pa->state.time, &sphdata);
  }
  
  /* Do global forces and effectors */  
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    BPH_sph_integrate(sim, pa, pa->state.time, &sphdata);
    
    if (sim->colliders)
      collision_check(sim, p, pa->state.time, cfra);
    
    /* SPH particles are not physical particles, just interpolation
     * particles, thus rotation has not a direct sense for them */
    basic_rotate(part, pa, pa->state.time, timestep);

#pragma omp critical
    if (part->time_flag & PART_TIME_AUTOSF)
      BPH_sph_update_courant_num(sim, pa, dtime, &sphdata);
  }

  BPH_psys_sph_finalise(&sphdata);
}
