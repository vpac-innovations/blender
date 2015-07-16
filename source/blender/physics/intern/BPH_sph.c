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
#include "DNA_object_refiner.h"
#include "DNA_scene_types.h"

#include "BLI_math.h"
#include "BLI_threads.h"
#include "BLI_kdopbvh.h"
#include "BLI_edgehash.h"
#include "BLI_listbase.h"

#include "BKE_particle.h"
#include "BKE_collision.h"
#include "BKE_bvhutils.h"

#include "BPH_sph.h"

static ThreadRWMutex psys_bvhtree_rwlock = BLI_RWLOCK_INITIALIZER;

#define PSYS_FLUID_SPRINGS_INITIAL_SIZE 256

/* Calculate the speed of the particle relative to the local scale of the
 * simulation. This should be called once per particle during a simulation
 * step, after the velocity has been updated. element_size defines the scale of                                                                        
 * the simulation, and is typically the distance to neighboring particles. */
static void sph_update_courant_num(ParticleSimulationData *sim, ParticleData *pa, float dtime, SPHData *sphdata)
{
  float relative_vel[3];
  float speed;

  sub_v3_v3v3(relative_vel, pa->prev_state.vel, sphdata->flow);
  speed = len_v3(relative_vel);
  if (sim->courant_num < speed * dtime / sphdata->element_size)
    sim->courant_num = speed * dtime / sphdata->element_size;
}

static ParticleSpring* sph_spring_add(ParticleSystem *psys, ParticleSpring *spring)
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

static EdgeHash* sph_springhash_build(ParticleSystem *psys)
{
  EdgeHash *springhash = NULL;
  ParticleSpring *spring;
  int i = 0;

  springhash = BLI_edgehash_new_ex(__func__, psys->tot_fluidsprings);

  for (i=0, spring=psys->fluid_springs; i<psys->tot_fluidsprings; i++, spring++)
    BLI_edgehash_insert(springhash, spring->particle_index[0], spring->particle_index[1], SET_INT_IN_POINTER(i+1));

  return springhash;
}

static void sph_spring_delete(ParticleSystem *psys, int j)
{
	if (j != psys->tot_fluidsprings - 1)
		psys->fluid_springs[j] = psys->fluid_springs[psys->tot_fluidsprings - 1];

	psys->tot_fluidsprings--;

	if (psys->tot_fluidsprings < psys->alloc_fluidsprings/2 && psys->alloc_fluidsprings > PSYS_FLUID_SPRINGS_INITIAL_SIZE) {
		psys->alloc_fluidsprings /= 2;
		psys->fluid_springs = (ParticleSpring*)MEM_reallocN(psys->fluid_springs,  psys->alloc_fluidsprings * sizeof(ParticleSpring));
	}
}

static void sph_springs_modify(ParticleSystem *psys, float dtime)
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

	/* Loop through springs backwards - for efficient delete function */
	for (i=psys->tot_fluidsprings-1; i >= 0; i--) {
		if (psys->fluid_springs[i].delete_flag)
			sph_spring_delete(psys, i);
	}
}

static void sph_evaluate_func(BVHTree *tree, ParticleSystem **psys, float co[3], SPHRangeData *pfr, float interaction_radius, BVHTree_RangeQuery callback)
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

static void sph_particle_courant(SPHData *sphdata, SPHRangeData *pfr)
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

static void sph_init(SPHData *sphdata) {
	SPHFluidSettings *fluid = sphdata->psys[0]->part->fluid;

	/* 4.77 is an experimentally determined density factor */
	sphdata->rest_density = fluid->rest_density * (fluid->flag & SPH_FAC_DENSITY ? 4.77f : 1.f);
	sphdata->stiffness = fluid->stiffness_k;
	sphdata->stiffness_near_fac = fluid->stiffness_knear * (fluid->flag & SPH_FAC_REPULSION ? fluid->stiffness_k : 1.f);
}

static void sph_density_accum_cb(void *userdata, int index, float squared_dist)
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

  pfr->params.density += pow2f(q);
  pfr->params.near_density += pow3f(q);
}

static void sph_equation_of_state(SPHData *sphdata, SPHParams *params) {
	params->pressure =  sphdata->stiffness * (params->density - sphdata->rest_density);
	params->near_pressure = sphdata->stiffness_near_fac * params->near_density;
}

static void sph_force_cb(void *sphdata_v, ParticleKey *state, float *force, float *UNUSED(impulse))
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

	float visc = fluid->viscosity_omega;
	float stiff_visc = fluid->viscosity_beta * (fluid->flag & SPH_FAC_VISCOSITY ? fluid->viscosity_omega : 1.f);

	float inv_mass = 1.0f / sphdata->mass;
	float spring_constant = fluid->spring_k;

	/* 4.0 seems to be a pretty good value */
	float interaction_radius = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * pa->size : 1.0f);
	float h = interaction_radius * sphdata->hfac;
	float rest_length = fluid->rest_length * (fluid->flag & SPH_FAC_REST_LENGTH ? 2.588f * pa->size : 1.f);

	ParticleData *npa;
	float vec[3];
	float vel[3];
	float co[3];

	int i, spring_index, index = pa - psys[0]->particles;

	pfr.params.density = pfr.params.near_density = 0.0f;
	pfr.h = h;
	pfr.pa = pa;
	pfr.mass = sphdata->mass;

	sph_evaluate_func( NULL, psys, state->co, &pfr, interaction_radius, sph_density_accum_cb);
	sph_equation_of_state(sphdata, &pfr.params);

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
		madd_v3_v3fl(force, vec, -(pfr.params.pressure + pfr.params.near_pressure*q)*q);

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
					sph_spring_add(psys[0], &temp_spring);
				}
			}
			else {/* PART_SPRING_HOOKES - Hooke's spring force */
				madd_v3_v3fl(force, vec, -10.f * spring_constant * (1.f - rij/h) * (rest_length - rij));
			}
		}
	}
	
	/* Artificial buoyancy force in negative gravity direction  */
	if (fluid->buoyancy > 0.f && gravity)
		madd_v3_v3fl(force, gravity, fluid->buoyancy * (pfr.params.density - sphdata->rest_density));

	if (sphdata->pass == 0 && psys[0]->part->time_flag & PART_TIME_AUTOSF)
		sph_particle_courant(sphdata, &pfr);
	sphdata->pass++;
}

static void sphclassical_init(SPHData *sphdata) {
	SPHFluidSettings *fluid = sphdata->psys[0]->part->fluid;

	/* 4.77 is an experimentally determined density factor */
	sphdata->rest_density = fluid->rest_density * (fluid->flag & SPH_FAC_DENSITY ? 4.77f : 1.0f);
	/* Use speed of sound squared */
	sphdata->stiffness = pow2f(fluid->stiffness_k);
}

static void sphclassical_density_accum_cb(void *userdata, int index, float squared_dist)
{
  SPHRangeData *pfr = (SPHRangeData *)userdata;
  ParticleData *npa = pfr->npsys->particles + index;
  float massfac = npa->sphmassfac;
  float q;
  float qfac = 21.0f / (256.f * (float)M_PI);
  float rij, rij_h;
  float vec[3];

  /* Exclude particles that are more than 2h away. Can't use squared_dist here                  
   * because it is not accurate enough. Use current state, i.e. the output of                   
   * basic_integrate() - z0r */
  if (pfr->pa == NULL) {
    rij = sqrtf(squared_dist);
  } else {
    sub_v3_v3v3(vec, npa->state.co, pfr->pa->state.co);
    rij = len_v3(vec);
  }
  rij_h = rij / (pfr->h * npa->sphalpha);
  if (rij_h > 2.0f)
  return;

  /* Smoothing factor. Use a 5th order Wendland kernel.
   * http://arxiv.org/abs/1204.2471
   * The `pow3` is because we're working in 3D space.
   *
   * gnuplot:
   *     q1(x) = (2.0 - x)**4 * ( 1.0 + 2.0 * x)
   *     plot [0:2] q1(x) */
  q  = qfac / pow3f(pfr->h * npa->sphalpha) * pow4f(2.0f - rij_h) * ( 1.0f + 2.0f * rij_h);
  q *= pfr->npsys->part->mass * massfac;
  
  if (pfr->use_size)
    q *= pfr->pa->size;

  pfr->params.density += q;
}

static void sphclassical_neighbour_accum_cb(void *userdata, int index, float squared_dist)
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
  if (pfr->pa == NULL) {
    rij = sqrtf(squared_dist);
  } else {
    sub_v3_v3v3(vec, npa->state.co, pfr->pa->state.co);
    rij = len_v3(vec);
  }
  rij_h = rij / pfr->h;
  if (rij_h > 2.0f)
    return;

  pfr->neighbors[pfr->tot_neighbors].index = index;
  pfr->neighbors[pfr->tot_neighbors].psys = pfr->npsys;
  pfr->tot_neighbors++;
}

static void sphclassical_equation_of_state(SPHData *sphdata, SPHParams *params) {
	params->pressure = sphdata->stiffness * (pow7f(params->density / sphdata->rest_density) - 1.0f);
}

static void sphclassical_force_cb(void *sphdata_v, ParticleKey *state, float *force, float *UNUSED(impulse))
{
	SPHData *sphdata = (SPHData *)sphdata_v;
	ParticleSystem **psys = sphdata->psys;
	ParticleData *pa = sphdata->pa;
	SPHFluidSettings *fluid = psys[0]->part->fluid;
	SPHRangeData pfr;
	SPHNeighbor *pfn;
	SPHParams nparams;
	float *gravity = sphdata->gravity;

	float dq, u, rij, dv[3];

	float visc = fluid->viscosity_omega;

	float interaction_radius;
	float h, hinv;

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

	sph_evaluate_func(NULL, psys, state->co, &pfr, interaction_radius, sphclassical_neighbour_accum_cb);

	/* Calculate pressure from density. Density is calculated in a preceding
	 * loop, so it has to be fetched here. */
	pfr.params.density = pa->sphdensity;
	sphclassical_equation_of_state(sphdata, &pfr.params);

	/* multiply by mass so that we return a force, not accel */
	qfac2 *= sphdata->mass * pa->sphmassfac;

	pfn = pfr.neighbors;
	for (i = 0; i < pfr.tot_neighbors; i++, pfn++) {
		npa = pfn->psys->particles + pfn->index;
		if (npa == pa) {
			/* we do not contribute to ourselves */
			continue;
		}

		/*Symmetrize kernel for particle splitting */
		pfr.h = (h/2) * (pa->sphalpha + npa->sphalpha);

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

		/* Equation of state: convert density of neighbor to pressure. */
		nparams.density = npa->sphdensity;
		sphclassical_equation_of_state(sphdata, &nparams);

		/* First derivative of smoothing factor. Utilise the Wendland kernel
		 * (see above).
		 *
		 * gnuplot:
		 *     dq(x) = (2.0 - x)**4 - 2.0 * (2.0 - x)**3 * (1.0 + 2.0 * x)
		 *     plot [0:2] dq(x)
		 * comparison with smoothing factor q1 (see above) in gnuplot:
		 *     plot [0:2] q(x), dq(x)
		 * Particles > 2h away are excluded above. */
		dq = (qfac2 / pow3f(pfr.h)) * (pow4f(2.0f - rij_h) - 2.0f * pow3f(2.0f - rij_h) * (1.0f + 2.0f * rij_h)  );
		dq *= sphdata->mass * npa->sphmassfac;

		if (pfn->psys->part->flag & PART_SIZEMASS)
			dq *= npa->size;

		pressureTerm = (pfr.params.pressure / pow2f(pfr.params.density) +
		                nparams.pressure / pow2f(nparams.density));

		/* Apply pressure of neighbor (scalar) as a force (vector). The
		 * total force is found by summing over all neighboring particles.
		 * Note that 'minus' is removed, because vec = vecBA, not vecAB.
		 * This applies to the viscosity calculation below, too.
		 *
		 * We multiply the pressure term by the gradient of the smoothing
		 * kernel `dq` to find the pressure gradient. */
		madd_v3_v3fl(force, vec, pressureTerm * dq);

		/* Viscosity */
		if (visc > 0.0f) {
			sub_v3_v3v3(dv, npa->prev_state.vel, pa->prev_state.vel);
			u = dot_v3v3(vec, dv);
			/* Apply parameters */
			u *= -dq * hinv * visc / (0.5f * nparams.density + 0.5f * pfr.params.density);
			madd_v3_v3fl(force, vec, u);
		}
	}

	/* Artificial buoyancy force in negative gravity direction  */
	if (fluid->buoyancy > 0.f && gravity)
		madd_v3_v3fl(force, gravity, fluid->buoyancy * (pfr.params.density - sphdata->rest_density));

	if (sphdata->pass == 0 && psys[0]->part->time_flag & PART_TIME_AUTOSF)
		sph_particle_courant(sphdata, &pfr);
	sphdata->pass++;
}

void psys_sph_init(ParticleSimulationData *sim, SPHData *sphdata)
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
	sphdata->eh = sph_springhash_build(sim->psys);

	// These per-particle values should be overridden later, but just for
	// completeness we give them default values now.
	sphdata->pa = NULL;
	sphdata->mass = 1.0f;

	if (sim->psys->part->fluid->solver == SPH_SOLVER_DDR) {
		sphdata->init = sph_init;
		sphdata->force_cb = sph_force_cb;
		sphdata->density_cb = sph_density_accum_cb;
		sphdata->equation_of_state = sph_equation_of_state;
		sphdata->hfac = 1.0f;
	}
	else {
		/* SPH_SOLVER_CLASSICAL */
		sphdata->init = sphclassical_init;
		sphdata->force_cb = sphclassical_force_cb;
		sphdata->density_cb = sphclassical_density_accum_cb;
		sphdata->equation_of_state = sphclassical_equation_of_state;
		sphdata->hfac = 0.5f;
	}

	sphdata->init(sphdata);
}

void psys_sph_finalise(SPHData *sphdata)
{
	if (sphdata->eh) {
		BLI_edgehash_free(sphdata->eh, NULL);
		sphdata->eh = NULL;
	}
}

void psys_deadpars_add(DeadParticles* deadpars, int index)
{
	int *new_data;

	/* If there is space in dead particle array add new index to the end. */
	if(deadpars->size < deadpars->capacity){
		*(deadpars->data + deadpars->size) = index;
		deadpars->size++;
		return;
	}

	/* Otherwise reallocate dead particle array then add new index to the end. */
	/* Allocate new memory */
	new_data = MEM_callocN((deadpars->capacity + 1000)*sizeof(int), "dead particles");

	/* Copy old data to newly allocated memory. */
	memcpy(new_data, deadpars->data, deadpars->capacity * sizeof(int));

	/* Free old memory. */
	MEM_freeN(deadpars->data);

	/* Set pointer to new memory. */
	deadpars->data = new_data;

	/* Add new index to end, update size and capacity. */
	*(deadpars->data + deadpars->size) = index;
	deadpars->size++;
	deadpars->capacity += 1000;
}

/* Sample the density field at a point in space. */
void psys_sph_sample(BVHTree *tree, SPHData *sphdata, float co[3], SPHParams *params)
{
  ParticleSystem **psys = sphdata->psys;
  SPHFluidSettings *fluid = psys[0]->part->fluid;
  /* 4.0 seems to be a pretty good value */
  float interaction_radius  = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * psys[0]->part->size : 1.0f);
  SPHRangeData pfr;

  pfr.params.density = pfr.params.near_density = 0.0f;
  pfr.params.pressure = pfr.params.near_pressure = 0.0f;
  pfr.h = interaction_radius * sphdata->hfac;
  pfr.mass = sphdata->mass;
  pfr.pa = NULL;

  sph_evaluate_func(tree, psys, co, &pfr, interaction_radius, sphdata->density_cb);
  sphdata->equation_of_state(sphdata, &pfr.params);
  memcpy(params, &pfr.params, sizeof(SPHParams));
}

static void sphclassical_calc_dens(ParticleData *pa, float UNUSED(dfra), SPHData *sphdata)
{
  ParticleSystem **psys = sphdata->psys;
  SPHFluidSettings *fluid = psys[0]->part->fluid;
  /* 4.0 seems to be a pretty good value */
  float interaction_radius  = fluid->radius * (fluid->flag & SPH_FAC_RADIUS ? 4.0f * psys[0]->part->size : 1.0f);
  SPHRangeData pfr;

  pfr.params.density = 0.0f;
  pfr.h = interaction_radius * sphdata->hfac;
  pfr.pa = pa;
  pfr.mass = sphdata->mass;

  sph_evaluate_func( NULL, psys, pa->state.co, &pfr, interaction_radius, sphclassical_density_accum_cb);
  pa->sphdensity = MIN2(MAX2(pfr.params.density, fluid->rest_density * 0.9f), fluid->rest_density * 1.1f);
}

static void sph_integrate(ParticleSimulationData *sim, ParticleData *pa, float dfra, SPHData *sphdata)
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

static int split_through_wall_test(ParticleSimulationData *sim, ParticleData *pa, BVHTreeRayHit *hit)
{
	ColliderCache *coll;
	ListBase *colliders = sim->colliders;
	BVHTreeFromMesh treeData = {NULL};
	Object *current;
	float ray_start[3], ray_end[3], ray_dir[3], old_dist;

	if(BLI_listbase_is_empty(colliders))
		return 0;

	copy_v3_v3(ray_start, pa->prev_state.co);
	copy_v3_v3(ray_end, pa->state.co);

	sub_v3_v3v3(ray_dir, ray_end, ray_start);
	hit->index = -1;
	hit->dist = len_v3(ray_dir);

	/* Iterate over colliders and check for intersect */
	for(coll=colliders->first; coll; coll=coll->next){
		current = coll->ob;
		old_dist = hit->dist;

		bvhtree_from_mesh_faces(&treeData, current->derivedFinal, 0.0f, 4, 6);

		/* Ray cast. */
		BLI_bvhtree_ray_cast(treeData.tree, ray_start, ray_dir, pa->size, hit, treeData.raycast_callback, &treeData);

		/* Throw out new hit distance if previous one was shorter. */
		if(old_dist < hit->dist)
			hit->dist = old_dist;

		free_bvhtree_from_mesh(&treeData);
	}
	return hit->index >= 0;
}

static int nearest_split(ParticleSimulationData *sim, SPHRangeData *pfr)
{
	/* Since particles array is not being resized, use nearest
	 * neighbour callback to efficiently find closest neighbour
	 * particle */
	ParticleSystem *psys = sim->psys;
	ParticleData *pa = pfr->pa, *npa, test_pa;
	BVHTreeRayHit hit;
	SPHNeighbor *pfn;
	float min_dist = 2.f * pfr->h;
	float dist, vec[3];
	int index=0, p;

	pfn = pfr->neighbors;
	for(p=0; p < pfr->tot_neighbors; p++, pfn++){
		npa = pfn->psys->particles + pfn->index;
		if(pa->sphmassfac+npa->sphmassfac >= 1.05f || npa->alive != PARS_ALIVE ||
		   ((npa->state.co[2] < 0.030) && (npa->state.co[2] > -0.030) &&
		    (npa->state.co[1] < 0.030 && npa->state.co[1] > -0.030) &&
		    (npa->state.co[0] < 0.030 && npa->state.co[0] > -0.030)) ||
			npa == pa){
			continue;
		}
/*
		copy_particle_key(&test_pa.state, &npa->state, 0);
		copy_particle_key(&test_pa.prev_state, &pa->state, 0);
		if(split_through_wall_test(sim, &test_pa, &hit))
			continue;
*/
		sub_v3_v3v3(vec, pa->state.co, npa->state.co);
		dist = len_v3(vec);
		if(dist < min_dist){
			min_dist = dist;
			index = pfn->index;
		}
	}

	return index;
}

void BPH_sph_unsplit_particle(ParticleSimulationData *sim, float cfra)
{
	SPHData sphdata;
	ParticleSystem *psys = sim->psys;
	ParticleData *npa;
	SPHRangeData pfr;
	float qfac = qfac = 21.0f / (256.f * (float)M_PI);
	float interaction_radius = psys->part->fluid->radius;
	float h = interaction_radius * 0.5f;
	float mass = psys->part->mass;
	float fac1, fac2, wma, wmb, rij_h, old_massfac, vec[3], old_co[3], old_vel[3];
	int index_n, i;
	PARTICLE_P;

	psys_sph_init(sim, &sphdata);

	LOOP_DYNAMIC_PARTICLES{
		if((pa->alive != PARS_ALIVE || pa->sphmassfac > 0.91f)){
			continue;
		}

		/* Check if particle is within range of a refiner */
		if(((pa->state.co[2] < 0.030) && (pa->state.co[2] > -0.030) &&
	     (pa->state.co[1] < 0.030 && pa->state.co[1] > -0.030) &&
	     (pa->state.co[0] < 0.030 && pa->state.co[0] > -0.030)))
			continue;
		/*if(((pa->state.co[2] < 0.10) && (pa->state.co[2] > -0.60) &&
		    (pa->state.co[1] < 2.000 && pa->state.co[1] > -2.000) &&
		    (pa->state.co[0] < 2.000 && pa->state.co[0] > -2.000)) ||
		   ((pa->state.co[2] < -3.90) && (pa->state.co[2] > -4.60) &&
		    (pa->state.co[1] < 2.000 && pa->state.co[1] > -2.000) &&
		    (pa->state.co[0] < 2.000 && pa->state.co[0] > -2.000))){
			continue;
		}*/
		pfr.h = h;
		pfr.pa = pa;

		/* Find nearest neighbours for current particle */
		sph_evaluate_func(NULL, sphdata.psys, pa->state.co, &pfr, interaction_radius, sphclassical_neighbour_accum_cb);

		/* Find index of nearest split particle */
		index_n = nearest_split(sim, &pfr);
		if(index_n == 0)
			continue;

		npa = psys->particles+index_n;
		old_massfac = pa->sphmassfac;
		copy_v3_v3(old_co, pa->state.co);
		copy_v3_v3(old_vel, pa->state.vel);

		/* Merge particles, kill unused particles*/
		/* -- Set massfac for merged particle. */
		pa->sphmassfac += npa->sphmassfac;

		fac1 = (npa->sphmassfac)/(pa->sphmassfac);
		fac2 = old_massfac/(pa->sphmassfac);

		/* -- Set position for merged particle */
		madd_v3_v3flv3fl(pa->state.co, npa->state.co, fac1, old_co, fac2);

		/* -- Set velocity for merged particle */
		madd_v3_v3flv3fl(pa->state.vel, npa->state.vel, fac1, old_vel, fac2);

		/* -- Set sphalpha for merged particle */
		sub_v3_v3v3(vec, pa->state.co, old_co);
		rij_h = len_v3(vec)/(h * pa->sphalpha);
		wma = old_massfac * mass * qfac / pow3f(h * pa->sphalpha) * pow4f(2.0f-rij_h) * (1.0f + 2.0f * rij_h);

		sub_v3_v3v3(vec, pa->state.co, npa->state.co);
		rij_h = len_v3(vec)/(h * npa->sphalpha);
		wmb = npa->sphmassfac * mass * qfac / pow3f(h * npa->sphalpha) * pow4f(2.0f-rij_h) * (1.0f + 2.0f * rij_h);

		fac1 = (21.f * pa->sphmassfac * mass)/(16.f * (float)M_PI * (wma + wmb));
		pa->sphalpha = pow(fac1 , 1.f/3.f)/h;

		/* -- Set state variables
			  Not checking for sphmassfac == 1 to
			  allow for precision errors.          */
		if(pa->sphmassfac >= 0.95f){
			pa->split = PARS_UNSPLIT;
			pa->sphmassfac = 1.f;
			pa->sphalpha = 1.f;
		}
		/* -- Kill other particle. */
		npa->dietime = cfra + 0.001/((float)(psys->part->subframes + 1));
		npa->alive = PARS_DEAD;
		psys_deadpars_add(&psys->deadpars, index_n);
	}
}

static void split_positions(ParticleSimulationData *sim, ParticleData *pa, int num)
{
	ParticleData test_pa;
	BVHTreeRayHit hit;
	float h = sim->psys->part->fluid->radius * 0.5f;
	float eps = 0.7f * h;
	float factor = eps * sqrt(3.f) / 3.f;
	int i;

	memcpy(&test_pa, pa, sizeof(ParticleData));
	for(i = 0; i < 3; i++){
		test_pa.prev_state.co[i] = test_pa.state.co[i];
	}

	switch(num){
		case 1:
			/* Test for split through wall */
			test_pa.state.co[0] += factor;
			test_pa.state.co[1] += factor;
			test_pa.state.co[2] += factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] += factor;
			pa->state.co[1] += factor;
			pa->state.co[2] += factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] += factor;
			pa->prev_state.co[1] += factor;
			pa->prev_state.co[2] += factor;
			break;
		case 2:
			/* Test for split through wall */
			test_pa.state.co[0] += factor;
			test_pa.state.co[1] += factor;
			test_pa.state.co[2] -= factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] += factor;
			pa->state.co[1] += factor;
			pa->state.co[2] -= factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] += factor;
			pa->prev_state.co[1] += factor;
			pa->prev_state.co[2] -= factor;
			break;
		case 3:
			test_pa.state.co[0] += factor;
			test_pa.state.co[1] -= factor;
			test_pa.state.co[2] -= factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] += factor;
			pa->state.co[1] -= factor;
			pa->state.co[2] -= factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] += factor;
			pa->prev_state.co[1] -= factor;
			pa->prev_state.co[2] -= factor;
			break;
		case 4:
			test_pa.state.co[0] -= factor;
			test_pa.state.co[1] -= factor;
			test_pa.state.co[2] -= factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] -= factor;
			pa->state.co[1] -= factor;
			pa->state.co[2] -= factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] -= factor;
			pa->prev_state.co[1] -= factor;
			pa->prev_state.co[2] -= factor;
			break;
		case 5:
			test_pa.state.co[0] -= factor;
			test_pa.state.co[1] += factor;
			test_pa.state.co[2] += factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] -= factor;
			pa->state.co[1] += factor;
			pa->state.co[2] += factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] -= factor;
			pa->prev_state.co[1] += factor;
			pa->prev_state.co[2] += factor;
			break;
		case 6:
			test_pa.state.co[0] -= factor;
			test_pa.state.co[1] -= factor;
			test_pa.state.co[2] += factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] -= factor;
			pa->state.co[1] -= factor;
			pa->state.co[2] += factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] -= factor;
			pa->prev_state.co[1] -= factor;
			pa->prev_state.co[2] += factor;
			break;
		case 7:
			test_pa.state.co[0] += factor;
			test_pa.state.co[1] -= factor;
			test_pa.state.co[2] += factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] += factor;
			pa->state.co[1] -= factor;
			pa->state.co[2] += factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] += factor;
			pa->prev_state.co[1] -= factor;
			pa->prev_state.co[2] += factor;
			break;
		case 8:
			test_pa.state.co[0] -= factor;
			test_pa.state.co[1] += factor;
			test_pa.state.co[2] -= factor;

			if(sim->colliders){
				i = split_through_wall_test(sim, &test_pa, &hit);
				if(i)
					factor = (hit.dist-1.f*pa->size)*sqrt(3.f) / 3.f;
			}

			pa->state.co[0] -= factor;
			pa->state.co[1] += factor;
			pa->state.co[2] -= factor;
			/* Not sure if this is appropriate.*/
			pa->prev_state.co[0] -= factor;
			pa->prev_state.co[1] += factor;
			pa->prev_state.co[2] -= factor;
			break;
		default:
			break;
	}
}

void BPH_sph_refiners_init(ListBase **refiners, ParticleSystem *psys)
{
	SPHRefiner *refiner= MEM_callocN(sizeof(SPHRefiner), "SPHRefiner");

	if(*refiners == NULL)
		*refiners = MEM_callocN(sizeof(ListBase), "refiners list");

	/* Add refiner(s) */
	refiner->co[0] = 0.0f;
	refiner->co[1] = 0.0f;
	refiner->co[2] = 0.0f;
	refiner->radius = 0.02f;

	BLI_addtail(*refiners, refiner);
}

void BPH_sph_refiners_end(ListBase **refiners)
{
	if (*refiners) {
		SPHRefiner *sref = (*refiners)->first;

		BLI_freelistN(*refiners);
		MEM_freeN(*refiners);
		*refiners = NULL;
	}
}

int check_refiners(ListBase *refiners, ParticleData *pa)
{
	SPHRefiner *sref;
	float vec[3], dist;

	if(refiners) for(sref = refiners->first; sref; sref=sref->next) {
		sub_v3_v3v3(vec, pa->state.co, sref->co);
		dist = len_v3(vec);
		if (dist < sref->radius)
			return 1;
	}

	return 0;
}

void BPH_sph_split_particle(ParticleSimulationData *sim, int index, float cfra)
{
	ParticleSystem *psys = sim->psys;
	ParticleSettings *part = psys->part;
	ParticleData *pa, *new_pa;
	int oldtotpart = psys->totpart;
	int newtotpart;
	int newparticles = 8-psys->deadpars.size;
	int i;

	pa = psys->particles+index;

	/* Check if particle is within a refinement zone */
	if(!check_refiners(psys->refiners, pa))
		return;

	if(pa->split == PARS_UNSPLIT){
		pa->split = PARS_SPLIT;
		pa->sphalpha = 0.75f;
		pa->sphmassfac = 0.2f;

		if(newparticles > 0){
			/* Re-allocate particles array. */
			newtotpart = oldtotpart+newparticles;
			realloc_particles(sim, newtotpart);
			pa = psys->particles+index;

			/* Make copies of parent particle at end of particles array. */
			for(i = 0; i < newparticles; i++){
				new_pa = psys->particles+oldtotpart+i;
				memcpy(new_pa, pa, sizeof(ParticleData));
				new_pa->sphmassfac = 0.1f;

				/* Set position for new particle. */
				split_positions(sim, new_pa, i+1);

				/* Set birth time. Offset to avoid particle reset. */
				psys -> particles[oldtotpart+i].time = cfra - 0.001/((float)(part->subframes + 1));
			}
			/* Update ParticleSettings->totpart. */
			psys->part->totpart = newtotpart;
			psys->totadded += newparticles;
		}
		else
			newparticles = 0;

		for(i = 0; i < 8 - newparticles; i++){
			new_pa = psys->particles+psys->deadpars.data[psys->deadpars.size-1-i];
			memcpy(new_pa, pa, sizeof(ParticleData));
			new_pa->sphmassfac = 0.1f;

			/* Set position for new particle. */
			split_positions(sim, new_pa, newparticles+i+1);

			/* Set birth time. Offset to avoid particle reset. */
			new_pa->alive = PARS_ALIVE;
			new_pa->time = cfra - 0.001/((float)(part->subframes + 1));
			//psys -> particles[oldtotpart+i].time = cfra - 0.001/((float)(part->subframes + 1));
		}
		psys->deadpars.size -= 8 - newparticles;
	}
}

void BPH_sph_planar_split(ParticleSimulationData *sim, int index, float cfra)
{
	ParticleSystem *psys = sim->psys;
	ParticleSettings *part = psys->part;
	ParticleData *pa, *new_pa;
	int oldtotpart = psys->totpart;
	int newtotpart = oldtotpart+2;
	int i;

	pa = psys->particles+index;

	/* Split particles in predefined box */
	if((pa->state.co[2] > 0.040 || pa->state.co[2] < -0.020) ||
	   (pa->state.co[1] > 0.020 || pa->state.co[1] < -0.020) ||
	   (pa->state.co[0] > 0.020 || pa->state.co[0] < -0.020))
	   return;

	if(pa->split == PARS_UNSPLIT){
		/* Re-allocate particles array */
		realloc_particles(sim, newtotpart);
		pa = psys->particles+index;

		pa->split = PARS_SPLIT;
		pa->sphalpha = pow(1.f/3.f, 1.f/3.f);
		pa->sphmassfac = 1.f/3.f;

		/* Make copies of parent particle at end of particles array */
		new_pa = psys->particles+oldtotpart;
		memcpy(new_pa, pa, sizeof(ParticleData));

		new_pa->state.co[0] -= (1.f/6.f)*new_pa->size;
		new_pa->state.co[1] -= (sqrt(3.f)/18.f)*new_pa->size;
		new_pa->prev_state.co[0] -= (1.f/6.f)*new_pa->size;
		new_pa->prev_state.co[1] -= (sqrt(3.f)/18.f)*new_pa->size;

		new_pa = psys->particles+oldtotpart+1;
		memcpy(new_pa, pa, sizeof(ParticleData));

		new_pa->state.co[0] += (1.f/6.f)*new_pa->size;
		new_pa->state.co[1] -= (sqrt(3.f)/18.f)*new_pa->size;
		new_pa->prev_state.co[0] += (1.f/6.f)*new_pa->size;
		new_pa->prev_state.co[1] -= (sqrt(3.f)/18.f)*new_pa->size;

		pa->state.co[1] += (sqrt(3.f)/9.f)*new_pa->size;
		pa->prev_state.co[1] += (sqrt(3.f)/9.f)*new_pa->size;

		/* Set position for new particle */
//		split_positions(sim, new_pa, i+1);

		/* Set birth time. Offset to avoid particle reset. Is this robust though?*/
		psys -> particles[oldtotpart].time = cfra - 0.001/((float)(part->subframes + 1));
		psys -> particles[oldtotpart+1].time = cfra - 0.001/((float)(part->subframes + 1));
		/* Update ParticleSettings->totpart.
				  ParticleSystem->totpart? */
		psys->part->totpart = newtotpart;
		psys->totadded += 2;
	}
}

void BPH_sphDDR_step(ParticleSimulationData *sim, float dtime, float cfra)
{
  SPHData sphdata;
  ParticleSettings *part = sim->psys->part;
  ParticleSystem *psys = sim->psys;
  float timestep;
  PARTICLE_P;

  timestep = psys_get_timestep(sim);
  psys_sph_init(sim, &sphdata);

  /* Apply SPH forces using double-density relaxation algorithm                                                                                                                
   * (Clavat et. al.) */
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    /* do global forces & effectors */
    basic_integrate(sim, p, pa->state.time, cfra);

    /* actual fluids calculations */
    sph_integrate(sim, pa, pa->state.time, &sphdata);

    if (sim->colliders)
      collision_check(sim, p, pa->state.time, cfra);

    /* SPH particles are not physical particles, just interpolation                                                                                                      
     * particles,  thus rotation has not a direct sense for them */
    basic_rotate(part, pa, pa->state.time, timestep);

#pragma omp critical
    if (part->time_flag & PART_TIME_AUTOSF)
      sph_update_courant_num(sim, pa, dtime, &sphdata);
  }

  sph_springs_modify(psys, timestep);
  psys_sph_finalise(&sphdata);
}

void BPH_sphclassical_step(ParticleSimulationData *sim, float dtime, float cfra)
{
  SPHData sphdata;
  ParticleSettings *part = sim->psys->part;
  ParticleSystem *psys = sim->psys;
  float timestep;
  PARTICLE_P;

  timestep = psys_get_timestep(sim);
  psys_sph_init(sim, &sphdata);

#pragma omp parallel for private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    basic_integrate(sim, p, pa->state.time, cfra);
  }
  
  /* Calculate summation density */
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    sphclassical_calc_dens(pa, pa->state.time, &sphdata);
  }
  
  /* Do global forces and effectors */  
#pragma omp parallel for firstprivate (sphdata) private (pa) schedule(dynamic,5)
  LOOP_DYNAMIC_PARTICLES {
    sph_integrate(sim, pa, pa->state.time, &sphdata);
    
    if (sim->colliders)
      collision_check(sim, p, pa->state.time, cfra);
    
    /* SPH particles are not physical particles, just interpolation
     * particles, thus rotation has not a direct sense for them */
    basic_rotate(part, pa, pa->state.time, timestep);

#pragma omp critical
    if (part->time_flag & PART_TIME_AUTOSF)
      sph_update_courant_num(sim, pa, dtime, &sphdata);
  }

  psys_sph_finalise(&sphdata);
}
