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
 * Contributor(s): Blender Foundation (2015), Sean Loh
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/makesrna/intern/rna_object_refine.c
 *  \ingroup RNA
 */

#include <stdlib.h>

#include "DNA_object_types.h"
#include "DNA_object_refiner.h"
#include "DNA_particle_types.h"
#include "DNA_scene_types.h"

#include "RNA_define.h"

#include "rna_internal.h"

#include "WM_api.h"
#include "WM_types.h"

#ifdef RNA_RUNTIME

#else
static void rna_def_field(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;

	static EnumPropertyItem field_type_items[] = {
	    {0, "NONE", 0, "None", ""},
		{REFINE_OBJ, "REFINER", ICON_PARTICLES, "Refiner", "Refiner for classical SPH adaptive resolution"},
	    {0, NULL, 0, NULL, NULL}
	};
/*
	srna = RNA_def_struct(brna, "FieldSettings", NULL);
	RNA_def_struct_sdna(srna, "PartRefine");
	RNA_def_struct_path_func(srna, "rna_FieldSettings_path");
	RNA_def_struct_ui_text(srna, "Field Settings", "Field settings for an object in physics simulation");
	RNA_def_struct_ui_icon(srna, ICON_PHYSICS);

	/* Enums *//*

	prop = RNA_def_property(srna, "refiner", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "flag");
	RNA_def_property_enum_items(prop, field_type_items);
	RNA_def_property_enum_funcs(prop, NULL, "rna_FieldSettings_type_set", NULL);
	RNA_def_property_ui_text(prop, "Refiner", "Is object a refiner");
	RNA_def_property_update(prop, 0, "rna_FieldSettings_dependency_update");*/
}

void RNA_def_object_refine(BlenderRNA *brna)
{
	rna_def_field(brna);
}

#endif
