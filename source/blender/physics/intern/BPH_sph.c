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
#include "BLI_listbase.h"

#include "BKE_particle.h"
#include "BKE_collision.h"
#include "BKE_bvhutils.h"

#include "BPH_sph.h"

static ThreadRWMutex psys_bvhtree_rwlock = BLI_RWLOCK_INITIALIZER;

#define PSYS_FLUID_SPRINGS_INITIAL_SIZE 256

static int split_through_wall_test(ParticleSimulationData *sim, ParticleData *pa, BVHTreeRayHit *hit)
{
	ColliderCache *coll;
	ListBase *colliders = sim->colliders;
	BVHTreeFromMesh treeData = {NULL};
	Object *current;
	float ray_start[3], ray_end[3], ray_dir[3];

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

		/* Create bvhtree from current object. Check last 3 args are appropriate. */
		/* 0.0f is epsilon parameter, if 0.0f, over-written with FLT_EPSILON	  */
		/* Next parameter is 'tree_type' not sure what that controls or how to	  */
		/* choose it appropriately for our purpose.								  */
		/* See BLI_bvhtree_new() for consequences of 'axis' parameter. (last par) */
		/* Not exactly sure what that controls either.							  */
		bvhtree_from_mesh_faces(&treeData, current->derivedFinal, 0.0f, 4, 6);

		/* Ray cast. Not sure what comes out of this at the moment. */
		BLI_bvhtree_ray_cast(treeData.tree, ray_start, ray_dir, pa->size, &hit, treeData.raycast_callback, &treeData);

		/* Free pointer to re-use? No idea if this works */
		/* Free's memory allocated by treeData, then sets sizeof(treeData) sized  */
		/* block of memory pointed by treeData to zero							  */
		/* bvhtree_from_mesh_faces() relies on an if(tree == NULL) so will that   */
		/* still work?															  */
		free_bvhtree_from_mesh(&treeData);
	}

	return hit->index >= 0;
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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

			i = split_through_wall_test(sim, &test_pa, &hit);
			if(i)
				factor = (hit.dist-2.f*pa->size)*sqrt(3.f) / 3.f;

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

void BPH_sph_split_particle(ParticleSimulationData *sim, int index, float cfra)
{
	ParticleSystem *psys = sim->psys;
	ParticleSettings *part = psys->part;
	ParticleData *pa, *new_pa;
	int oldtotpart = psys->totpart;
	int newtotpart = oldtotpart+8;
	int i;

	pa = psys->particles+index;

	/* Split particles in predefined box */
	if((pa->state.co[2] >= -0.5 || pa->state.co[2] <= -0.6) ||
	   (pa->state.co[1] <= 0.4 || pa->state.co[1] >= 0.7) ||
	   (pa->state.co[0] >= -0.5 || pa->state.co[0] <= -0.6))
		return;

	if(pa->split == PARS_UNSPLIT){
		/* Re-allocate particles array */
		realloc_particles(sim, newtotpart);
		pa = psys->particles+index;

		pa->split = PARS_SPLIT;
		pa->sphalpha = 0.75f;
		pa->sphmassfac = 0.253311f;

		/* Make copies of parent particle at end of particles array */
		for(i = 0; i < 8; i++){
			new_pa = psys->particles+oldtotpart+i;
			memcpy(new_pa, pa, sizeof(ParticleData));
			new_pa->sphmassfac = 0.0933361f;

			/* Set position for new particle */
			split_positions(sim, new_pa, i+1);

			/* Set birth time. Offset to avoid particle reset. Is this robust though?*/
			psys -> particles[oldtotpart+i].time = cfra - 0.001/((float)(part->subframes + 1));
		}
		/* Update ParticleSettings->totpart.
				  ParticleSystem->totpart? */
		psys->part->totpart = newtotpart;
		psys->totsplit += 1;
	}
}

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
