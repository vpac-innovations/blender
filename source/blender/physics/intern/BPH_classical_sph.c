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

/** \file blender/physics/intern/BPH_classical_sph.c
 *  \ingroup bph
 */

#include <stdio.h>

#include "DNA_particle_types.h"

#include "BLI_math.h"
#include "BLI_threads.h"
#include "BLI_kdopbvh.h"

#include "BKE_particle.h"

#include "BPH_classical_sph.h"

static ThreadRWMutex psys_bvhtree_rwlock = BLI_RWLOCK_INITIALIZER;

void BPH_sphclassical_density_accum_cb(void *userdata, 
				       int   index, 
				       float UNUSED(squared_dist))
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

void BPH_sphclassical_neighbour_accum_cb(void *userdata, 
					 int   index, 
					 float UNUSED(squared_dist))
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
void BPH_sphclassical_force_cb(void        *sphdata_v, 
			       ParticleKey *state, 
			       float       *force, 
			       float       *UNUSED(impulse))
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
		dq = qfac2 * (2.0f * pow4f(2.0f - rij_h) - 4.0f * pow3f(2.0f - rij_h) * (1.0f + 2.0f * rij_h)  );

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

void BPH_sphclassical_calc_dens(ParticleData *pa, 
				float         UNUSED(dfra), 
				SPHData       *sphdata)
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

void BPH_sph_evaluate_func(BVHTree           *tree,
			   ParticleSystem   **psys,
			   float              co[3],
			   SPHRangeData      *pfr,
			   float              interaction_radius, 
			   BVHTree_RangeQuery callback)
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

void BPH_sph_particle_courant(SPHData      *sphdata, 
			      SPHRangeData *pfr)
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
/*void BPH_sph_integrate(ParticleSimulationData *sim, 
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
  /* restore previous state and treat gravity & effectors as external acceleration*//*
  sub_v3_v3v3(effector_acceleration, pa->state.vel, pa->prev_state.vel);
  mul_v3_fl(effector_acceleration, 1.f/dtime);

  copy_particle_key(&pa->state, &pa->prev_state, 0);

  integrate_particle(part, pa, dtime, effector_acceleration, sphdata->force_cb, sphdata);
  }*/
