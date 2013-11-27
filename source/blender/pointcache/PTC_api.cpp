/*
 * Copyright 2013, Blender Foundation.
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
 */

#include "PTC_api.h"

#ifdef WITH_ALEMBIC

#include "reader.h"
#include "writer.h"
#include "export.h"
#include "particles.h"

#include "util_path.h"

extern "C" {
#include "DNA_object_types.h"
#include "DNA_particle_types.h"
#include "DNA_pointcache_types.h"

#include "RNA_access.h"
}

using namespace PTC;

void PTC_writer_free(PTCWriter *_writer)
{
	PTC::Writer *writer = (PTC::Writer *)_writer;
	delete writer;
}

void PTC_write_sample(struct PTCWriter *_writer)
{
	PTC::Writer *writer = (PTC::Writer *)_writer;
	writer->write_sample();
}

void PTC_bake(struct Main *bmain, struct Scene *scene, struct PTCWriter *_writer, int start_frame, int end_frame,
              short *stop, short *do_update, float *progress)
{
	PTC::Writer *writer = (PTC::Writer *)_writer;
	PTC::Exporter exporter(bmain, scene, stop, do_update, progress);
	exporter.bake(writer, start_frame, end_frame);
}


void PTC_reader_free(PTCReader *_reader)
{
	PTC::Reader *reader = (PTC::Reader *)_reader;
	delete reader;
}

void PTC_read_sample(struct PTCReader *_reader)
{
	PTC::Reader *reader = (PTC::Reader *)_reader;
	reader->read_sample();
}

/* get writer/reader from RNA type */
PTCWriter *PTC_writer_from_rna(Scene *scene, PointerRNA *ptr)
{
	if (RNA_struct_is_a(ptr->type, &RNA_ParticleSystem)) {
		Object *ob = (Object *)ptr->id.data;
		ParticleSystem *psys = (ParticleSystem *)ptr->data;
		return PTC_writer_particles(scene, ob, psys);
	}
	return NULL;
}

PTCReader *PTC_reader_from_rna(Scene *scene, PointerRNA *ptr)
{
	if (RNA_struct_is_a(ptr->type, &RNA_ParticleSystem)) {
		Object *ob = (Object *)ptr->id.data;
		ParticleSystem *psys = (ParticleSystem *)ptr->data;
		return PTC_reader_particles(scene, ob, psys);
	}
	return NULL;
}

/* Particles */
PTCWriter *PTC_writer_particles(Scene *scene, Object *ob, ParticleSystem *psys)
{
	PointCache *cache = psys->pointcache;
	if (!cache)
		return NULL;
	std::string filename = ptc_archive_path(cache->name, cache->index, cache->path, &ob->id,
	                                        true, true,
	                                        cache->flag & PTCACHE_EXTERNAL,
	                                        cache->flag & PTCACHE_IGNORE_LIBPATH);
	
	PTC::ParticlesWriter *writer = new PTC::ParticlesWriter(filename, scene, ob, psys);
	return (PTCWriter *)writer;
}

PTCReader *PTC_reader_particles(Scene *scene, Object *ob, ParticleSystem *psys)
{
	PointCache *cache = psys->pointcache;
	if (!cache)
		return NULL;
	std::string filename = ptc_archive_path(cache->name, cache->index, cache->path, &ob->id,
	                                        true, true,
	                                        cache->flag & PTCACHE_EXTERNAL,
	                                        cache->flag & PTCACHE_IGNORE_LIBPATH);

	PTC::ParticlesReader *reader = new PTC::ParticlesReader(filename, scene, ob, psys);
	return (PTCReader *)reader;
}

#else

void PTC_writer_free(PTCWriter *_writer)
{
}

void PTC_write(struct PTCWriter *_writer)
{
}

PTCWriter *PTC_writer_create_particles(const char *filename, struct Object *ob, struct ParticleSystem *psys)
{
	return NULL;
}

#endif