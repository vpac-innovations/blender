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

        split = layout.split(percentage=0.3)
        split.label(text="Type:")
        split.prop(refiner, "type", text="")

        split = layout.split(percentage=0.3)
        split.label(text="Split ratio:")
        split.prop(refiner, "ratio", text="")

        split = layout.split(percentage=0.3)
        split.label(text="No. of splits:")
        split.prop(refiner, "nsplit", text="")

        split = layout.split()
        split.prop(refiner, "use_falloff")

        basic_refiner_settings_ui(self, context, refiner)

if __name__ == "__main__":  # only for live edit.
    bpy.utils.register_module(__name__)
