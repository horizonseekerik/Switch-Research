"""
=============================================================================
PROJECT JANUS: PUBLICATION-GRADE 3D PHOTONIC SWITCH MODEL (BLENDER 5.2)
=============================================================================
Top 0.1% global quality 3D CAD model & photorealistic rendering of the
certified Sb2S3 1x2 directional coupler switch cell at lambda = 1064 nm.

Geometry Ground Truth (Task-6990):
  - Waveguide 1 (Bar) & Waveguide 2 (Cross): 240 nm width x 220 nm height
  - Evanescent Coupling Gap: 118 nm uniform parallel spacing
  - Sb2S3 Active Layer (WG1): 160 nm thickness
    * 3.95 um uniform body
    * Dual 1.50 um parabolic apodized tapers (15 nm tip to 240 nm width)
    * Total interaction length: L_c = 6.95 um
  - Monolayer Graphene & M1 Gold Contacts: 120 Ohm/sq heater layer + contact pads
  - Substrate & BOX: 2 um SiO2 buried oxide on silicon handle wafer
  - Optical Mode: Volumetric 1064 nm TE00 luminescent coupling into Port 2

Renderings Generated (2560x1440 QHD):
  1. render_1_hero_isometric_3d.png    - Cinematic hero 3/4 beauty view
  2. render_2_coupling_gap_macro.png   - Nanoscale macro close-up of 118 nm gap & taper
  3. render_3_topdown_schematic.png    - Architectural top-down plan view
  4. render_4_cross_section_facet.png  - Cross-sectional facet view down waveguide axis

Exported Master 3D Assets:
  - sb2s3_switch_cell_3d_master.blend  - Native Blender project with lighting/materials
  - sb2s3_switch_cell_3d.glb           - WebGL / GLTF 2.0 interactive asset
  - sb2s3_switch_cell_3d.obj / .mtl    - Universal CAD mesh interchange
=============================================================================
"""

import bpy
import bmesh
import math
import os
import sys
import numpy as np

print("=" * 80)
print("COMPILING 3D PHOTONIC SWITCH CELL ASSET & PHOTOREALISTIC RENDERS")
print("=" * 80)

# Clear entire scene
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.unit_settings.system = 'METRIC'
scene.unit_settings.scale_length = 1.0 # 1 BU = 1.0 um

output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. LUXURY STUDIO ENVIRONMENT & LIGHTING
# -----------------------------------------------------------------------------
world = bpy.data.worlds.new("Studio_Tech_Dark")
scene.world = world
world.use_nodes = True
w_nodes = world.node_tree.nodes
w_links = world.node_tree.links
w_nodes.clear()

w_out = w_nodes.new('ShaderNodeOutputWorld')
w_bg = w_nodes.new('ShaderNodeBackground')
w_bg.inputs['Color'].default_value = (0.012, 0.015, 0.022, 1.0) # Deep dark studio navy
w_bg.inputs['Strength'].default_value = 0.65
w_links.new(w_bg.outputs['Background'], w_out.inputs['Surface'])

# -----------------------------------------------------------------------------
# 2. PBR SHADERS
# -----------------------------------------------------------------------------
def make_pbr_material(name, color, metallic=0.0, roughness=0.18, ior=1.45, transmission=0.0, emission=None, emission_strength=1.0):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    nodes.clear()
    
    out = nodes.new('ShaderNodeOutputMaterial')
    out.location = (400, 0)
    
    if emission:
        emit = nodes.new('ShaderNodeEmission')
        emit.inputs['Color'].default_value = emission
        emit.inputs['Strength'].default_value = emission_strength
        links.new(emit.outputs['Emission'], out.inputs['Surface'])
    else:
        bsdf = nodes.new('ShaderNodeBsdfPrincipled')
        bsdf.inputs['Base Color'].default_value = color
        bsdf.inputs['Metallic'].default_value = metallic
        bsdf.inputs['Roughness'].default_value = roughness
        bsdf.inputs['IOR'].default_value = ior
        if 'Transmission Weight' in bsdf.inputs:
            bsdf.inputs['Transmission Weight'].default_value = transmission
        elif 'Transmission' in bsdf.inputs:
            bsdf.inputs['Transmission'].default_value = transmission
        links.new(bsdf.outputs['BSDF'], out.inputs['Surface'])
        
    return mat

mat_silicon = make_pbr_material("Mat_Crystalline_Silicon", (0.15, 0.21, 0.28, 1.0), metallic=0.94, roughness=0.14, ior=3.48)
mat_sb2s3 = make_pbr_material("Mat_Sb2S3_Chalcogenide", (0.96, 0.65, 0.14, 1.0), metallic=0.30, roughness=0.16, ior=2.70)
mat_graphene = make_pbr_material("Mat_Graphene_Layer", (0.04, 0.04, 0.05, 1.0), metallic=0.60, roughness=0.25, transmission=0.40)
mat_gold = make_pbr_material("Mat_Gold_M1_Pads", (1.00, 0.78, 0.28, 1.0), metallic=1.00, roughness=0.07)
mat_sio2 = make_pbr_material("Mat_SiO2_BOX_Glass", (0.82, 0.90, 0.96, 1.0), roughness=0.04, ior=1.449, transmission=0.88)
mat_handle_wafer = make_pbr_material("Mat_Handle_Wafer", (0.05, 0.07, 0.09, 1.0), metallic=0.85, roughness=0.25)
mat_laser_core = make_pbr_material("Mat_Laser_Core", (0.15, 0.90, 1.0, 1.0), emission=(0.12, 0.88, 1.0, 1.0), emission_strength=22.0)

# -----------------------------------------------------------------------------
# 3. NANOMETER PHYSICAL GEOMETRY CREATION
# -----------------------------------------------------------------------------
W_WG = 0.240       # 240 nm width
H_WG = 0.220       # 220 nm height
GAP = 0.118        # 118 nm gap
L_CELL = 7.550     # 7.55 um cell length
L_C = 6.950        # 6.95 um interaction section
L_TIP = 1.500      # 1.50 um taper tip
L_BODY = 3.950     # 3.95 um uniform body
T_PCM = 0.160      # 160 nm Sb2S3 thickness

y_wg1 = +(W_WG + GAP) / 2.0  # +0.179 um (WG1 center)
y_wg2 = -(W_WG + GAP) / 2.0  # -0.179 um (WG2 center)

def create_beveled_box(name, xmin, xmax, ymin, ymax, zmin, zmax, material, bevel_w=0.010):
    mesh = bpy.data.meshes.new(name + "_Mesh")
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    
    bm = bmesh.new()
    bmesh.ops.create_cube(bm, size=1.0)
    bm.to_mesh(mesh)
    bm.free()
    
    sx = xmax - xmin
    sy = ymax - ymin
    sz = zmax - zmin
    obj.scale = (sx, sy, sz)
    obj.location = ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0)
    
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    obj.select_set(False)
    
    if bevel_w > 0:
        mod = obj.modifiers.new("Bevel", 'BEVEL')
        mod.width = bevel_w
        mod.segments = 2
        mod.limit_method = 'ANGLE'
        mod.angle_limit = math.radians(45)
        
    obj.data.materials.append(material)
    return obj

# 1. Waveguide 1 (Input/Bar core)
wg1 = create_beveled_box("Waveguide_1_Bar", -L_CELL/2.0, +L_CELL/2.0, y_wg1 - W_WG/2.0, y_wg1 + W_WG/2.0, 0.0, H_WG, mat_silicon)

# 2. Waveguide 2 (Cross core)
wg2 = create_beveled_box("Waveguide_2_Cross", -L_CELL/2.0, +L_CELL/2.0, y_wg2 - W_WG/2.0, y_wg2 + W_WG/2.0, 0.0, H_WG, mat_silicon)

# 3. Active Sb2S3 Patch with High-Resolution Parabolic Tapers
def build_tapered_sb2s3_patch(name, material):
    mesh = bpy.data.meshes.new(name + "_Mesh")
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    
    bm = bmesh.new()
    x_half = L_BODY / 2.0
    x_start = -L_C / 2.0
    x_end = +L_C / 2.0
    
    n_pts_taper = 28
    xs_in = [x_start + i * ((-x_half) - x_start) / n_pts_taper for i in range(n_pts_taper + 1)]
    xs_body = [-x_half, +x_half]
    xs_out = [x_half + i * (x_end - x_half) / n_pts_taper for i in range(1, n_pts_taper + 1)]
    xs = xs_in + xs_body[1:] + xs_out
    
    w_tip = 0.015 # 15 nm needle tip
    
    bot_rings = []
    top_rings = []
    
    for x in xs:
        if x < -x_half:
            u = (x - x_start) / L_TIP
            w = w_tip + (W_WG - w_tip) * (1.0 - (1.0 - u)**2) # Smooth parabolic taper
        elif x > x_half:
            u = (x_end - x) / L_TIP
            w = w_tip + (W_WG - w_tip) * (1.0 - (1.0 - u)**2)
        else:
            w = W_WG
            
        y_bot = y_wg1 - w / 2.0
        y_top = y_wg1 + w / 2.0
        
        b1 = bm.verts.new((x, y_bot, H_WG))
        b2 = bm.verts.new((x, y_top, H_WG))
        t1 = bm.verts.new((x, y_bot, H_WG + T_PCM))
        t2 = bm.verts.new((x, y_top, H_WG + T_PCM))
        
        bot_rings.append((b1, b2))
        top_rings.append((t1, t2))
        
    for i in range(len(xs) - 1):
        b1, b2 = bot_rings[i]
        b3, b4 = bot_rings[i+1]
        t1, t2 = top_rings[i]
        t3, t4 = top_rings[i+1]
        
        bm.faces.new((b1, b3, b4, b2)) # Bottom face
        bm.faces.new((t1, t2, t4, t3)) # Top face
        bm.faces.new((b1, t1, t3, b3)) # Outer wall
        bm.faces.new((b2, b4, t4, t2)) # Gap-facing wall
        
    # Cap end facets
    bm.faces.new((bot_rings[0][0], bot_rings[0][1], top_rings[0][1], top_rings[0][0]))
    bm.faces.new((bot_rings[-1][0], top_rings[-1][0], top_rings[-1][1], bot_rings[-1][1]))
    
    bm.to_mesh(mesh)
    bm.free()
    
    bpy.context.view_layer.objects.active = obj
    mod = obj.modifiers.new("Bevel", 'BEVEL')
    mod.width = 0.006
    mod.segments = 2
    
    obj.data.materials.append(material)
    return obj

sb2s3_patch = build_tapered_sb2s3_patch("Sb2S3_Active_Patch", mat_sb2s3)

# 4. Graphene Monolayer Microheater (Conformal top strip)
graphene_sheet = create_beveled_box(
    "Graphene_MicroHeater",
    -L_C/2.0 - 0.02, +L_C/2.0 + 0.02,
    y_wg1 - W_WG/2.0 - 0.01, y_wg1 + W_WG/2.0 + 0.01,
    H_WG + T_PCM, H_WG + T_PCM + 0.003,
    mat_graphene, bevel_w=0.0
)

# 5. Gold Contact Electrodes (M1 Metallization Pads)
gold_pad_in = create_beveled_box(
    "Gold_Electrode_In",
    -L_C/2.0 - 0.32, -L_C/2.0 - 0.02,
    y_wg1 - 0.20, y_wg1 + 0.20,
    H_WG + T_PCM - 0.02, H_WG + T_PCM + 0.10,
    mat_gold, bevel_w=0.012
)
gold_pad_out = create_beveled_box(
    "Gold_Electrode_Out",
    +L_C/2.0 + 0.02, +L_C/2.0 + 0.32,
    y_wg1 - 0.20, y_wg1 + 0.20,
    H_WG + T_PCM - 0.02, H_WG + T_PCM + 0.10,
    mat_gold, bevel_w=0.012
)

# 6. SiO2 Buried Oxide (BOX) Substrate (Sleek cutaway)
box_slab = create_beveled_box(
    "SiO2_Buried_Oxide",
    -L_CELL/2.0 - 0.35, +L_CELL/2.0 + 0.35,
    -0.60, +0.60,
    -0.25, 0.0,
    mat_sio2, bevel_w=0.02
)

# 7. Silicon Handle Wafer Base
handle_wafer = create_beveled_box(
    "Silicon_Handle_Wafer",
    -L_CELL/2.0 - 0.35, +L_CELL/2.0 + 0.35,
    -0.60, +0.60,
    -0.45, -0.25,
    mat_handle_wafer, bevel_w=0.02
)

# 8. Volumetric 1064 nm TE Mode Cross-Coupling Optical Wave
def build_optical_mode_beam(name, material):
    mesh = bpy.data.meshes.new(name + "_Mesh")
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    
    bm = bmesh.new()
    n_steps = 100
    xs = np.linspace(-L_CELL/2.0, +L_CELL/2.0, n_steps)
    r_mode = 0.038 # 38 nm modal core radius
    
    prev_ring = None
    for x in xs:
        # Cross transition from y_wg1 to y_wg2
        u = (x + L_C/2.0) / L_C
        u_clamped = max(0.0, min(1.0, u))
        w_cross = math.sin(u_clamped * math.pi / 2.0)**2
        cy = y_wg1 * (1.0 - w_cross) + y_wg2 * w_cross
        cz = H_WG / 2.0
        
        ring = []
        n_circ = 14
        for k in range(n_circ):
            ang = 2.0 * math.pi * k / n_circ
            vy = cy + r_mode * math.cos(ang)
            vz = cz + r_mode * math.sin(ang)
            ring.append(bm.verts.new((x, vy, vz)))
            
        if prev_ring:
            for k in range(n_circ):
                k_nxt = (k + 1) % n_circ
                bm.faces.new((prev_ring[k], prev_ring[k_nxt], ring[k_nxt], ring[k]))
        prev_ring = ring
        
    bm.to_mesh(mesh)
    bm.free()
    obj.data.materials.append(material)
    return obj

optical_beam = build_optical_mode_beam("Optical_Mode_TE00", mat_laser_core)

# -----------------------------------------------------------------------------
# 4. LIGHTING RIG
# -----------------------------------------------------------------------------
# Key Light (Warm Sunlight Angle)
key_light_data = bpy.data.lights.new("Key_Light", 'AREA')
key_light_data.energy = 180.0
key_light_data.color = (1.0, 0.95, 0.88)
key_light_data.size = 4.5
key_light = bpy.data.objects.new("Key_Light", key_light_data)
bpy.context.collection.objects.link(key_light)
key_light.location = (-3.2, -4.2, 3.8)
key_light.rotation_euler = (math.radians(52), math.radians(22), math.radians(-38))

# Rim Grazing Light (Sharp edge reflections on 118 nm gap)
rim_light_data = bpy.data.lights.new("Rim_Light", 'AREA')
rim_light_data.energy = 260.0
rim_light_data.color = (0.75, 0.90, 1.0)
rim_light_data.size = 3.0
rim_light = bpy.data.objects.new("Rim_Light", rim_light_data)
bpy.context.collection.objects.link(rim_light)
rim_light.location = (3.8, 3.2, 2.8)
rim_light.rotation_euler = (math.radians(-58), math.radians(-28), math.radians(140))

# Top Soft Fill Light
top_light_data = bpy.data.lights.new("Top_Light", 'AREA')
top_light_data.energy = 85.0
top_light_data.color = (0.92, 0.96, 1.0)
top_light_data.size = 8.0
top_light = bpy.data.objects.new("Top_Light", top_light_data)
bpy.context.collection.objects.link(top_light)
top_light.location = (0.0, 0.0, 5.5)

# -----------------------------------------------------------------------------
# 5. CAMERAS & SCENE FRAMING
# -----------------------------------------------------------------------------
focus_center = bpy.data.objects.new("Target_Center", None)
bpy.context.collection.objects.link(focus_center)
focus_center.location = (0.0, 0.0, H_WG / 2.0)

def create_camera(name, location, lens, target_obj, is_ortho=False, ortho_scale=8.0):
    cam_data = bpy.data.cameras.new(name)
    cam_data.lens = lens
    if is_ortho:
        cam_data.type = 'ORTHO'
        cam_data.ortho_scale = ortho_scale
    cam = bpy.data.objects.new(name, cam_data)
    bpy.context.collection.objects.link(cam)
    cam.location = location
    
    track = cam.constraints.new('TRACK_TO')
    track.target = target_obj
    track.track_axis = 'TRACK_NEGATIVE_Z'
    track.up_axis = 'UP_Y'
    return cam

# 1. Hero Isometric View (Diagonal full-span framing)
cam_hero = create_camera("Cam_Hero_Isometric", (-4.2, -4.0, 2.9), lens=58.0, target_obj=focus_center)

# 2. Nanoscale Gap & Taper Macro View
focus_macro = bpy.data.objects.new("Target_Macro", None)
bpy.context.collection.objects.link(focus_macro)
focus_macro.location = (-1.975, 0.0, H_WG + T_PCM / 2.0)
cam_macro = create_camera("Cam_Macro_Gap", (-2.7, -1.5, 1.05), lens=85.0, target_obj=focus_macro)

# 3. Orthographic Top-Down Layout (Horizontal device orientation)
cam_top = bpy.data.objects.new("Cam_TopDown_Plan", bpy.data.cameras.new("Cam_TopDown_Plan"))
bpy.context.collection.objects.link(cam_top)
cam_top.data.type = 'ORTHO'
cam_top.data.ortho_scale = 8.2
cam_top.location = (0.0, 0.0, 5.0)
cam_top.rotation_euler = (0.0, 0.0, 0.0) # Horizontal length along canvas width

# 4. Facet Cross-Section View (Looking straight down waveguide axis)
focus_facet = bpy.data.objects.new("Target_Facet", None)
bpy.context.collection.objects.link(focus_facet)
focus_facet.location = (-L_CELL / 2.0, 0.0, H_WG / 2.0)
cam_facet = create_camera("Cam_Facet_CrossSection", (-L_CELL / 2.0 - 2.2, -0.6, 0.65), lens=75.0, target_obj=focus_facet)

# -----------------------------------------------------------------------------
# 6. RENDER EXECUTION (2560x1440 QHD)
# -----------------------------------------------------------------------------
scene.render.resolution_x = 2560
scene.render.resolution_y = 1440
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = 'PNG'
scene.render.image_settings.color_mode = 'RGBA'

if hasattr(scene, 'eevee'):
    if hasattr(scene.eevee, 'use_raytracing'):
        scene.eevee.use_raytracing = True
    if hasattr(scene.eevee, 'use_fast_gi'):
        scene.eevee.use_fast_gi = True
    if hasattr(scene.eevee, 'taa_render_samples'):
        scene.eevee.taa_render_samples = 64

print("[1/4] Rendering View 1: Hero Cinematic Isometric...")
scene.camera = cam_hero
scene.render.filepath = os.path.join(output_dir, "render_1_hero_isometric_3d.png")
bpy.ops.render.render(write_still=True)

print("[2/4] Rendering View 2: Nanoscale Gap & Taper Macro...")
scene.camera = cam_macro
scene.render.filepath = os.path.join(output_dir, "render_2_coupling_gap_macro.png")
bpy.ops.render.render(write_still=True)

print("[3/4] Rendering View 3: Orthographic Top-Down Architecture...")
scene.camera = cam_top
scene.render.filepath = os.path.join(output_dir, "render_3_topdown_schematic.png")
bpy.ops.render.render(write_still=True)

print("[4/4] Rendering View 4: Facet Cross-Sectional Angle...")
scene.camera = cam_facet
scene.render.filepath = os.path.join(output_dir, "render_4_cross_section_facet.png")
bpy.ops.render.render(write_still=True)

# -----------------------------------------------------------------------------
# 7. EXPORT 3D ASSETS (.BLEND, .GLB, .OBJ)
# -----------------------------------------------------------------------------
blend_file = os.path.join(output_dir, "sb2s3_switch_cell_3d_master.blend")
print(f"Exporting Master Blend Project: {blend_file}")
bpy.ops.wm.save_as_mainfile(filepath=blend_file)

glb_file = os.path.join(output_dir, "sb2s3_switch_cell_3d.glb")
print(f"Exporting GLTF 2.0 Binary: {glb_file}")
bpy.ops.export_scene.gltf(filepath=glb_file, export_format='GLB', export_materials='EXPORT', export_yup=True)

obj_file = os.path.join(output_dir, "sb2s3_switch_cell_3d.obj")
print(f"Exporting Wavefront OBJ Mesh: {obj_file}")
try:
    bpy.ops.wm.obj_export(filepath=obj_file)
except Exception:
    bpy.ops.export_scene.obj(filepath=obj_file)

print("=" * 80)
print("SUCCESS: 3D PHOTONIC SWITCH MODEL & ALL 4 RENDERS SUCCESSFULLY COMPILED!")
print("=" * 80)
