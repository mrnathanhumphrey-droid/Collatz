"""
render_viz.py - Build qx+1 phase-space 3D visualization in Blender.

Run:
  blender --background --python render_viz.py

Loads vis_a.csv, vis_b.csv, vis_b_div.csv from C:/Collatz/visualization/.
Builds three Collections, each with one Mesh whose vertices are the data points.
Geometry Nodes converts vertices to render-points with vertex-color emission.
Caches CSV->NPZ on first run for fast subsequent rebuilds.
Renders 7 snapshots from informative angles to snapshots/.
Saves qx1_3d.blend.
"""
import csv
import math
import sys
from pathlib import Path

import bpy
import numpy as np
from mathutils import Vector

VIZ_DIR = Path("C:/Collatz/visualization")
SNAP_DIR = VIZ_DIR / "snapshots"
SNAP_DIR.mkdir(exist_ok=True)

# Coordinate scaling (shared so the three clouds align in z)
X_SCALE = 2.5   # log-value 0..43 -> 0..107 BU
Z_SCALE = 4.0   # q 3..25       -> 12..100 BU

# Per-viz Y scaling (different units, per brief)
Y_SCALE_A = 0.4   # step idx 0..200 -> 0..80
Y_SCALE_B = 0.25  # sigma 0..381    -> 0..95
Y_SCALE_BD = 2.5  # log final 0..43 -> 0..107
SENTINEL_Y_B = -20.0  # sentinel parking lane below the data plane

POINT_RADIUS = 0.25

# a★ power -> RGBA (perceptual spectrum, 0..6)
A_PALETTE = np.array([
    [1.00, 1.00, 1.00, 1.0],  # 0: white   (a★ = 1)
    [1.00, 0.20, 0.20, 1.0],  # 1: red     (a★ = q)
    [1.00, 0.55, 0.15, 1.0],  # 2: orange  (a★ = q^2)
    [1.00, 0.95, 0.20, 1.0],  # 3: yellow  (a★ = q^3)
    [0.25, 0.95, 0.30, 1.0],  # 4: green   (a★ = q^4)
    [0.25, 0.55, 1.00, 1.0],  # 5: blue    (a★ = q^5)
    [0.75, 0.30, 1.00, 1.0],  # 6: purple  (a★ = q^6)
], dtype=np.float32)
GREY = np.array([0.5, 0.5, 0.5, 1.0], dtype=np.float32)


def csv_to_arrays(csv_path, cache_path):
    """Load CSV columns. Cache as .npz on first run."""
    if cache_path.exists():
        d = np.load(cache_path)
        return d["x"], d["y"], d["z"], d["a_pow"], d["status_code"]

    print(f"  parsing {csv_path.name} ...", flush=True)
    status_map = {"converged": 0, "divergent": 1, "timeout": 2}
    xs, ys, zs, aps, scs = [], [], [], [], []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            zs.append(int(row["z"]))
            aps.append(int(row["a_star_power"]))
            scs.append(status_map[row["status"]])
    arr_x = np.asarray(xs, dtype=np.float32)
    arr_y = np.asarray(ys, dtype=np.float32)
    arr_z = np.asarray(zs, dtype=np.int16)
    arr_ap = np.asarray(aps, dtype=np.int8)
    arr_sc = np.asarray(scs, dtype=np.int8)
    np.savez(cache_path, x=arr_x, y=arr_y, z=arr_z, a_pow=arr_ap, status_code=arr_sc)
    print(f"  cached -> {cache_path.name}  ({len(arr_x):,} rows)", flush=True)
    return arr_x, arr_y, arr_z, arr_ap, arr_sc


def make_point_cloud_object(name, x, y, z, a_pow, y_scale, sentinel_y=None):
    """Create a mesh-as-vertices object with a vertex color attribute."""
    n = len(x)
    y_mapped = y.astype(np.float32) * y_scale
    if sentinel_y is not None:
        sent_mask = (y < 0)
        y_mapped[sent_mask] = sentinel_y

    coords = np.empty((n, 3), dtype=np.float32)
    coords[:, 0] = x * X_SCALE
    coords[:, 1] = y_mapped
    coords[:, 2] = z.astype(np.float32) * Z_SCALE

    a_idx = np.clip(a_pow.astype(np.int32), 0, len(A_PALETTE) - 1)
    fallback_mask = (a_pow < 0)
    colors = A_PALETTE[a_idx]
    colors[fallback_mask] = GREY  # defensive, shouldn't trigger

    mesh = bpy.data.meshes.new(name)
    # Bulk vertex creation via foreach_set is fastest for millions of points
    mesh.vertices.add(n)
    mesh.vertices.foreach_set("co", coords.flatten())
    mesh.update()

    col_attr = mesh.color_attributes.new(name="Col", type="FLOAT_COLOR", domain="POINT")
    col_attr.data.foreach_set("color", colors.flatten())

    obj = bpy.data.objects.new(name, mesh)
    return obj


def setup_geo_nodes_points(obj, point_radius=POINT_RADIUS):
    """GeometryNodes: Mesh -> Points (with radius) -> Set Material -> Output."""
    mat = bpy.data.materials.new(name=f"{obj.name}_mat")
    mat.use_nodes = True
    mt = mat.node_tree
    for n in list(mt.nodes):
        mt.nodes.remove(n)
    attr = mt.nodes.new("ShaderNodeAttribute")
    attr.attribute_name = "Col"
    attr.location = (-400, 0)
    emit = mt.nodes.new("ShaderNodeEmission")
    emit.inputs["Strength"].default_value = 1.0
    emit.location = (-100, 0)
    out = mt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (200, 0)
    mt.links.new(attr.outputs["Color"], emit.inputs["Color"])
    mt.links.new(emit.outputs["Emission"], out.inputs["Surface"])

    mod = obj.modifiers.new(name="PointCloud", type="NODES")
    ng = bpy.data.node_groups.new(name=f"{obj.name}_geo", type="GeometryNodeTree")
    # Blender 4+ interface API
    ng.interface.new_socket(name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    ng.interface.new_socket(name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    nodes = ng.nodes
    links = ng.links
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
    # Wipe default scene contents without read_factory_settings (which fails on
    # broken addon unregister hooks in this install).
    for obj in list(bpy.data.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    for col in list(bpy.data.collections):
        bpy.data.collections.remove(col)
    scene = bpy.context.scene
    scene.render.engine = "BLENDER_EEVEE"
    scene.render.resolution_x = 1600
    scene.render.resolution_y = 1000
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    # Disable bloom-style effects in Eevee that wash out colors
    if hasattr(scene, "eevee"):
        for attr in ("use_bloom", "use_gtao", "use_ssr"):
            if hasattr(scene.eevee, attr):
                setattr(scene.eevee, attr, False)
    # Linear color management so emission colors render close to source RGB
    scene.view_settings.view_transform = "Standard"
    # Black background
    bg = scene.world or bpy.data.worlds.new("World")
    scene.world = bg
    bg.use_nodes = True
    bg_node = bg.node_tree.nodes.get("Background")
    if bg_node:
        bg_node.inputs["Color"].default_value = (0.02, 0.02, 0.03, 1.0)
        bg_node.inputs["Strength"].default_value = 1.0
    # No light needed (emission materials)


def make_camera(name, location, look_at, lens=35):
    cam_data = bpy.data.cameras.new(name)
    cam_data.lens = lens
    cam = bpy.data.objects.new(name, cam_data)
    cam.location = location
    direction = Vector(look_at) - Vector(location)
    cam.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.collection.objects.link(cam)
    return cam


def render_to(filepath, cam, hide_collections=()):
    bpy.context.scene.camera = cam
    # Toggle visibility
    for c in bpy.data.collections:
        c.hide_render = (c.name in hide_collections)
    bpy.context.scene.render.filepath = str(filepath)
    bpy.ops.render.render(write_still=True)


def main():
    setup_scene()

    # Load data (cache CSV -> NPZ first time)
    print("[load] Vis A (trajectory cloud)", flush=True)
    xa, ya, za, apa, sca = csv_to_arrays(VIZ_DIR / "vis_a.csv", VIZ_DIR / "vis_a.npz")
    print("[load] Vis B (stopping-time cloud)", flush=True)
    xb, yb, zb, apb, scb = csv_to_arrays(VIZ_DIR / "vis_b.csv", VIZ_DIR / "vis_b.npz")
    print("[load] Vis B' (divergence-aware cloud)", flush=True)
    xbd, ybd, zbd, apbd, scbd = csv_to_arrays(VIZ_DIR / "vis_b_div.csv", VIZ_DIR / "vis_b_div.npz")

    # Build three collections
    coll_a = bpy.data.collections.new("Vis_A_Trajectories")
    coll_b = bpy.data.collections.new("Vis_B_StoppingTime")
    coll_bd = bpy.data.collections.new("Vis_B_DivergenceAware")
    for c in (coll_a, coll_b, coll_bd):
        bpy.context.scene.collection.children.link(c)

    print(f"[mesh] building Vis A object ({len(xa):,} verts)", flush=True)
    obj_a = make_point_cloud_object("Vis_A", xa, ya, za, apa, Y_SCALE_A)
    coll_a.objects.link(obj_a)
    setup_geo_nodes_points(obj_a, point_radius=0.20)

    print(f"[mesh] building Vis B object ({len(xb):,} verts)", flush=True)
    obj_b = make_point_cloud_object("Vis_B", xb, yb, zb, apb, Y_SCALE_B, sentinel_y=SENTINEL_Y_B)
    coll_b.objects.link(obj_b)
    setup_geo_nodes_points(obj_b, point_radius=0.45)

    print(f"[mesh] building Vis B' object ({len(xbd):,} verts)", flush=True)
    obj_bd = make_point_cloud_object("Vis_B_div", xbd, ybd, zbd, apbd, Y_SCALE_BD)
    coll_bd.objects.link(obj_bd)
    setup_geo_nodes_points(obj_bd, point_radius=0.45)

    # Cameras / snapshots
    # Scene roughly spans x ~ 0..107, y ~ -20..107, z ~ 12..100. Center ~ (50, 40, 56).
    center = np.array([50, 40, 56])

    cam_iso = make_camera("cam_iso", (260, -200, 220), center)
    cam_front = make_camera("cam_front", (55, -240, 55), center)         # looking +Y, see X-Z
    cam_side = make_camera("cam_side", (280, 40, 55), center)            # looking -X, see Y-Z
    cam_top = make_camera("cam_top", (55, 40, 280), center)              # looking -Z, see X-Y
    cam_low = make_camera("cam_low", (260, -180, 40), (55, 50, 60))      # low iso

    # Save .blend BEFORE rendering (lets the user open it even if render fails)
    blend_path = VIZ_DIR / "qx1_3d.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    print(f"[save] {blend_path}", flush=True)

    # Render snapshots
    snapshots = [
        ("01_all_iso.png",       cam_iso,   ()),
        ("02_all_front.png",     cam_front, ()),
        ("03_all_side.png",      cam_side,  ()),
        ("04_all_top.png",       cam_top,   ()),
        ("05_visA_only_iso.png", cam_iso,   ("Vis_B_StoppingTime", "Vis_B_DivergenceAware")),
        ("06_visB_only_iso.png", cam_iso,   ("Vis_A_Trajectories", "Vis_B_DivergenceAware")),
        ("07_visBdiv_only_iso.png", cam_iso, ("Vis_A_Trajectories", "Vis_B_StoppingTime")),
        ("08_low_angle_all.png", cam_low,   ()),
    ]
    for fname, cam, hide in snapshots:
        out = SNAP_DIR / fname
        print(f"[render] {fname}  cam={cam.name}  hidden={hide}", flush=True)
        render_to(out, cam, hide_collections=hide)
        print(f"           -> {out}", flush=True)

    print("\n[done] all renders complete", flush=True)


if __name__ == "__main__":
    main()
