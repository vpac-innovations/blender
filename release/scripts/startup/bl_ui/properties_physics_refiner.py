# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>
import bpy
from bpy.types import Panel

from bl_ui.properties_physics_common import (
        basic_refiner_settings_ui,
        #basic_refiner_falloff_ui,
        )


class PhysicButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "physics"

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return (context.object) and (not rd.use_game_engine)


class PHYSICS_PT_refiner(PhysicButtonsPanel, Panel):
    bl_label = "Refiner"

    @classmethod
    def poll(cls, context):
        ob = context.object
        rd = context.scene.render
        return (not rd.use_game_engine) and (ob.refiner) and (ob.refiner.type != 'NONE')

    def draw(self, context):
        layout = self.layout

        ob = context.object
        refiner = ob.refiner

        split = layout.split(percentage=0.2)
        split.label(text="Type:")

        split.prop(refiner, "type", text="")

        basic_refiner_settings_ui(self, context, refiner)
"""
        if refiner.type == 'POINT':
            split = layout.split(percentage=0.2)
            split.label(text="Shape:")
            split.prop(refiner, "shape", text="")
        elif refiner.type == 'SURFACE':
            return  # not implemented yet
        else:
            return  # nothing to draw

        split = layout.split()"""
"""
        else:
            basic_refiner_settings_ui(self, context, refiner)
        elif field.type == 'REFINER':
            col = split.column()
            col.prop(field, "guide_minimum")
            col.prop(field, "guide_free")
            col.prop(field, "falloff_power")
            col.prop(field, "use_guide_path_add")
            col.prop(field, "use_guide_path_weight")

            col = split.column()
            col.label(text="Clumping:")
            col.prop(field, "guide_clump_amount")
            col.prop(field, "guide_clump_shape")

            row = layout.row()
            row.prop(field, "use_max_distance")
            sub = row.row()
            sub.active = field.use_max_distance
            sub.prop(field, "distance_max")

            layout.separator()

            layout.prop(field, "guide_kink_type")
            if field.guide_kink_type != 'NONE':
                layout.prop(field, "guide_kink_axis")

                split = layout.split()

                col = split.column()
                col.prop(field, "guide_kink_frequency")
                col.prop(field, "guide_kink_shape")

                col = split.column()
                col.prop(field, "guide_kink_amplitude")

        elif field.type == 'TEXTURE':
            col = split.column()
            col.prop(field, "strength")
            col.prop(field, "texture_mode", text="")
            col.prop(field, "texture_nabla")

            col = split.column()
            col.prop(field, "use_object_coords")
            col.prop(field, "use_2d_force")
        elif field.type == 'SMOKE_FLOW':
            col = split.column()
            col.prop(field, "strength")
            col.prop(field, "flow")
            col = split.column()
            col.label(text="Domain Object:")
            col.prop(field, "source_object", "")
            col.prop(field, "use_smoke_density")
        else:
            basic_force_field_settings_ui(self, context, field)

        if field.type not in {'NONE', 'GUIDE'}:

            layout.label(text="Falloff:")
            layout.prop(field, "falloff_type", expand=True)

            basic_force_field_falloff_ui(self, context, field)

            if field.falloff_type == 'CONE':
                layout.separator()

                split = layout.split(percentage=0.35)

                col = split.column()
                col.label(text="Angular:")
                col.prop(field, "use_radial_min", text="Use Minimum")
                col.prop(field, "use_radial_max", text="Use Maximum")

                col = split.column()
                col.prop(field, "radial_falloff", text="Power")

                sub = col.column()
                sub.active = field.use_radial_min
                sub.prop(field, "radial_min", text="Angle")

                sub = col.column()
                sub.active = field.use_radial_max
                sub.prop(field, "radial_max", text="Angle")

            elif field.falloff_type == 'TUBE':
                layout.separator()

                split = layout.split(percentage=0.35)

                col = split.column()
                col.label(text="Radial:")
                col.prop(field, "use_radial_min", text="Use Minimum")
                col.prop(field, "use_radial_max", text="Use Maximum")

                col = split.column()
                col.prop(field, "radial_falloff", text="Power")

                sub = col.column()
                sub.active = field.use_radial_min
                sub.prop(field, "radial_min", text="Distance")

                sub = col.column()
                sub.active = field.use_radial_max
                sub.prop(field, "radial_max", text="Distance")
"""
if __name__ == "__main__":  # only for live edit.
    bpy.utils.register_module(__name__)
