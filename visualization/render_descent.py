"""
render_descent.py - Render Vis 1 (a* descent channels) + Vis 2 (sigma surface).

Run:
  blender --background --python render_descent.py

Loads viz_outputs/descent_a.csv, descent_b.csv. Each goes in its own Collection
with its own coordinate scaling. Six snapshots per viz to snapshots_descent/.
Saves descent_3d.blend.
"""
import csv
from pathlib import Path

import bpy
import numpy as np
from mathutils import Vector

VIZ_DIR = Path("C:/Collatz/visualization")
DATA_DIR = VIZ_DIR / "viz_outputs"
SNAP_DIR = VIZ_DIR / "snapshots_descent"
SNAP_DIR.mkdir(exist_ok=True)

# Vis 1 (descent channels by a*): x=log_val, y=step, z=a_star_idx
V1_X_SCALE = 4.0
V1_Y_SCALE = 1.0
V1_Z_SCALE = 20.0  # widely separate the 6 z-channels
V1_RADIUS = 0.18

# Vis 2 (sigma cloud): x=log_n, y=sigma, z=log(peak/n)
V2_X_SCALE = 8.0
V2_Y_SCALE = 0.25
V2_Z_SCALE = 9.0
V2_RADIUS = 0.40

# Vis 2 lives at +200 BU in X to keep clouds visually separable in same blend
V2_OFFSET_X = 220.0

# a* index -> color (0..6, matches generator output)
A_PALETTE = np.array([
    [1.00, 1.00, 1.00, 1.0],  # 0: white      (unused, a*=1 has 0 classes)
    [1.00, 0.30, 0.30, 1.0],  # 1: red        (a* = 3)
    [1.00, 0.65, 0.15, 1.0],  # 2: orange     (a* = 9)
    [0.95, 0.85, 0.10, 1.0],  # 3: gold       (a* = 27) - darker than yellow for visibility
    [0.30, 0.95, 0.40, 1.0],  # 4: green      (a* = 81)
    [0.20, 0.65, 1.00, 1.0],  # 5: blue       (a* = 243)
    [0.85, 0.35, 1.00, 1.0],  # 6: magenta    (a* = 729)
], dtype=np.float32)
GREY = np.array([0.5, 0.5, 0.5, 1.0], dtype=np.float32)


def csv_to_arrays(csv_path, cache_path):
    if cache_path.exists():
        d = np.load(cache_path)
        return d["x"], d["y"], d["z"], d["a_idx"]
    print(f"  parsing {csv_path.name} ...", flush=True)
    xs, ys, zs, ais = [], [], [], []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            zs.append(float(row["z"]))
            ais.append(int(row["a_star_idx"]))
    arr_x = np.asarray(xs, dtype=np.float32)
    arr_y = np.asarray(ys, dtype=np.float32)
    arr_z = np.asarray(zs, dtype=np.float32)
    arr_ai = np.asarray(ais, dtype=np.int8)
    np.savez(cache_path, x=arr_x, y=arr_y, z=arr_z, a_idx=arr_ai)
    print(f"  cached -> {cache_path.name} ({len(arr_x):,} rows)", flush=True)
    return arr_x, arr_y, arr_z, arr_ai


def make_point_cloud_object(name, x, y, z, a_idx,
                            x_scale, y_scale, z_scale,
                            x_offset=0.0, y_offset=0.0, z_offset=0.0):
    n = len(x)
    coords = np.empty((n, 3), dtype=np.float32)
    coords[:, 0] = x * x_scale + x_offset
    coords[:, 1] = y * y_scale + y_offset
    coords[:, 2] = z * z_scale + z_offset

    a_clip = np.clip(a_idx.astype(np.int32), 0, len(A_PALETTE) - 1)
    fallback_mask = (a_idx < 0)
    colors = A_PALETTE[a_clip]
    colors[fallback_mask] = GREY

    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(n)
    mesh.vertices.foreach_set("co", coords.flatten())
    mesh.update()
    col_attr = mesh.color_attributes.new(name="Col", type="FLOAT_COLOR", domain="POINT")
    col_attr.data.foreach_set("color", colors.flatten())

    obj = bpy.data.objects.new(name, mesh)
    return obj


def setup_geo_nodes_points(obj, point_radius):
    mat = bpy.data.materials.new(name=f"{obj.name}_mat")
    mat.use_nodes = True
    mt = mat.node_tree
    for n in list(mt.nodes):
        mt.nodes.remove(n)
    attr = mt.nodes.new("ShaderNodeAttribute"); attr.attribute_name = "Col"; attr.location = (-400, 0)
    emit = mt.nodes.new("ShaderNodeEmission"); emit.inputs["Strength"].default_value = 1.0; emit.location = (-100, 0)
    out = mt.nodes.new("ShaderNodeOutputMaterial"); out.location = (200, 0)
    mt.links.new(attr.outputs["Color"], emit.inputs["Color"])
    mt.links.new(emit.outputs["Emission"], out.inputs["Surface"])

    mod = obj.modifiers.new(name="PointCloud", type="NODES")
    ng = bpy.data.node_groups.new(name=f"{obj.name}_geo", type="GeometryNodeTree")
    ng.interface.new_socket(name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    ng.interface.new_socket(name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    nodes = ng.nodes; links = ng.links
    inp = nodes.new("NodeGroupInput"); inp.location = (-600, 0)
    onp = nodes.new("NodeGroupOutput"); onp.location = (600, 0)
    m2p = nodes.new("GeometryNodeMeshToPoints"); m2p.location = (-300, 0)
    m2p.inputs["Radius"].default_value = point_radius
    sm = nodes.new("GeometryNodeSetMaterial"); sm.location = (200, 0)
    sm.inputs["Material"].default_value = mat
    links.new(inp.outputs["Geometry"], m2p.inputs["Mesh"])
    links.new(m2p.outputs["Points"], sm.inputs["Geometry"])
    links.new(sm.outputs["Geometry"], onp.inputs["Geometry"])
    mod.node_group = ng


def setup_scene():
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for col in list(bpy.data.collections):
        bpy.data.collections.remove(col)
    scene = bpy.context.scene
    scene.render.engine = "BLENDER_EEVEE"
    scene.render.resolution_x = 1600
    scene.render.resolution_y = 1000
    scene.render.image_settings.file_format = "PNG"
    if hasattr(scene, "eevee"):
        for attr in ("use_bloom", "use_gtao", "use_ssr"):
            if hasattr(scene.eevee, attr):
                setattr(scene.eevee, attr, False)
    scene.view_settings.view_transform = "Standard"
    bg = scene.world or bpy.data.worlds.new("World")
    scene.world = bg
    bg.use_nodes = True
    bgn = bg.node_tree.nodes.get("Background")
    if bgn:
        bgn.inputs["Color"].default_value = (0.02, 0.02, 0.03, 1.0)


def make_camera(name, location, look_at, lens=35):
    cd = bpy.data.cameras.new(name); cd.lens = lens
    cam = bpy.data.objects.new(name, cd)
    cam.location = location
    cam.rotation_euler = (Vector(look_at) - Vector(location)).to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.collection.objects.link(cam)
    return cam


def render_to(filepath, cam, hide_collections=()):
    bpy.context.scene.camera = cam
    for c in bpy.data.collections:
        c.hide_render = (c.name in hide_collections)
    bpy.context.scene.render.filepath = str(filepath)
    bpy.ops.render.render(write_still=True)


def main():
    setup_scene()

    print("[load] Vis 1 descent (a* channels)", flush=True)
    x1, y1, z1, ai1 = csv_to_arrays(DATA_DIR / "descent_a.csv", DATA_DIR / "descent_a.npz")
    print("[load] Vis 2 (sigma surface)", flush=True)
    x2, y2, z2, ai2 = csv_to_arrays(DATA_DIR / "descent_b.csv", DATA_DIR / "descent_b.npz")

    coll_1 = bpy.data.collections.new("Vis1_Descent_aStar")
    coll_2 = bpy.data.collections.new("Vis2_Sigma_Surface")
    bpy.context.scene.collection.children.link(coll_1)
    bpy.context.scene.collection.children.link(coll_2)

    print(f"[mesh] Vis 1 ({len(x1):,} verts)", flush=True)
    obj1 = make_point_cloud_object("Vis1_Descent", x1, y1, z1, ai1,
                                   V1_X_SCALE, V1_Y_SCALE, V1_Z_SCALE)
    coll_1.objects.link(obj1)
    setup_geo_nodes_points(obj1, V1_RADIUS)

    print(f"[mesh] Vis 2 ({len(x2):,} verts)", flush=True)
    obj2 = make_point_cloud_object("Vis2_Sigma", x2, y2, z2, ai2,
                                   V2_X_SCALE, V2_Y_SCALE, V2_Z_SCALE,
                                   x_offset=V2_OFFSET_X)
    coll_2.objects.link(obj2)
    setup_geo_nodes_points(obj2, V2_RADIUS)

    # Vis 1 bounds: x ~ 0..100, y ~ 0..100, z ~ 0..120
    v1_center = np.array([50, 50, 60])
    cam_v1_iso = make_camera("cam_v1_iso", (200, -150, 200), v1_center)
    cam_v1_front = make_camera("cam_v1_front", (50, -200, 60), v1_center)  # see x-z (descent shape per channel)
    cam_v1_side = make_camera("cam_v1_side", (220, 50, 60), v1_center)     # see y-z (channel separation along step)
    cam_v1_top = make_camera("cam_v1_top", (50, 50, 240), v1_center)       # see x-y (descent in log_value-step)

    # Vis 2 bounds (after offset): x ~ V2_OFFSET..V2_OFFSET+115, y ~ 0..117, z ~ 0..103
    v2_center = np.array([V2_OFFSET_X + 55, 60, 50])
    cam_v2_iso = make_camera("cam_v2_iso", (V2_OFFSET_X + 200, -150, 200), v2_center)
    cam_v2_front = make_camera("cam_v2_front", (V2_OFFSET_X + 55, -200, 50), v2_center)
    cam_v2_side = make_camera("cam_v2_side", (V2_OFFSET_X + 220, 60, 50), v2_center)
    cam_v2_top = make_camera("cam_v2_top", (V2_OFFSET_X + 55, 60, 220), v2_center)

    blend_path = VIZ_DIR / "descent_3d.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    print(f"[save] {blend_path}", flush=True)

    snaps = [
        ("v1_01_iso.png",   cam_v1_iso,   ("Vis2_Sigma_Surface",)),
        ("v1_02_front.png", cam_v1_front, ("Vis2_Sigma_Surface",)),
        ("v1_03_side.png",  cam_v1_side,  ("Vis2_Sigma_Surface",)),
        ("v1_04_top.png",   cam_v1_top,   ("Vis2_Sigma_Surface",)),
        ("v2_01_iso.png",   cam_v2_iso,   ("Vis1_Descent_aStar",)),
        ("v2_02_front.png", cam_v2_front, ("Vis1_Descent_aStar",)),
        ("v2_03_side.png",  cam_v2_side,  ("Vis1_Descent_aStar",)),
        ("v2_04_top.png",   cam_v2_top,   ("Vis1_Descent_aStar",)),
    ]
    for fname, cam, hide in snaps:
        out = SNAP_DIR / fname
        print(f"[render] {fname}  cam={cam.name}", flush=True)
        render_to(out, cam, hide_collections=hide)

    print("\n[done] all renders complete", flush=True)


if __name__ == "__main__":
    main()
