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

#ifndef PTC_API_H
#define PTC_API_H

#ifdef __cplusplus
extern "C" {
#endif

struct Main;
struct Scene;
struct Object;
struct ParticleSystem;
struct PointerRNA;

struct PTCWriter;
struct PTCReader;

void PTC_writer_free(struct PTCWriter *writer);
void PTC_write_sample(struct PTCWriter *writer);

void PTC_reader_free(struct PTCReader *reader);
void PTC_read_sample(struct PTCReader *reader);

void PTC_bake(struct Main *bmain, struct Scene *scene, struct PTCWriter *writer, int start_frame, int end_frame,
              short *stop, short *do_update, float *progress);


/* get writer/reader from RNA type */
struct PTCWriter *PTC_writer_from_rna(struct Scene *scene, struct PointerRNA *ptr);
struct PTCReader *PTC_reader_from_rna(struct Scene *scene, struct PointerRNA *ptr);

/* Particles */
struct PTCWriter *PTC_writer_particles(struct Scene *scene, struct Object *ob, struct ParticleSystem *psys);
struct PTCReader *PTC_reader_particles(struct Scene *scene, struct Object *ob, struct ParticleSystem *psys);

#ifdef __cplusplus
} /* extern C */
#endif

#endif  /* PTC_API_H */