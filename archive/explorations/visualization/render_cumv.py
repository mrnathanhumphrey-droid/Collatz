"""
render_cumv.py - Render Vis 1 colored by cumulative-v residual.

Layout: same as descent_3d.blend Vis 1 (x=log_val, y=step, z=a_star_idx),
but instead of categorical a* coloring, color by:
  cumv_resid = cum_v - 2*n_odd_steps  (Cramer race deviation)

Negative resid (orbit "behind" heuristic, fewer halves than expected) -> cool.
Positive resid (orbit "ahead", more halves than expected) -> warm.
"""
import csv
from pathlib import Path

import bpy
import numpy as np
from mathutils import Vector

VIZ_DIR = Path("C:/Collatz/visualization")
DATA_DIR = VIZ_DIR / "viz_outputs"
SNAP_DIR = VIZ_DIR / "snapshots_cumv"
SNAP_DIR.mkdir(exist_ok=True)

X_SCALE = 4.0
Y_SCALE = 1.0
Z_SCALE = 20.0
POINT_RADIUS = 0.18

CUMV_RESID_MIN = -25.0
CUMV_RESID_MAX = 12.0


def csv_to_arrays(csv_path, cache_path):
    if cache_path.exists():
        d = np.load(cache_path)
        return d["x"], d["y"], d["z"], d["resid"]
    print(f"  parsing {csv_path.name} ...", flush=True)
    xs, ys, zs, rs = [], [], [], []
    with open(csv_path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            zs.append(float(row["z"]))
            rs.append(int(row["cumv_resid"]))
    arr_x = np.asarray(xs, dtype=np.float32)
    arr_y = np.asarray(ys, dtype=np.float32)
    arr_z = np.asarray(zs, dtype=np.float32)
    arr_r = np.asarray(rs, dtype=np.int32)
    np.savez(cache_path, x=arr_x, y=arr_y, z=arr_z, resid=arr_r)
    print(f"  cached -> {cache_path.name} ({len(arr_x):,} rows)", flush=True)
    return arr_x, arr_y, arr_z, arr_r


def diverging_palette(values, vmin, vmax):
    """Map values to diverging RGBA: blue (vmin) -> white (0) -> red (vmax)."""
    n = len(values)
    out = np.ones((n, 4), dtype=np.float32)
    v = np.clip(values.astype(np.float32), vmin, vmax)
    neg = v < 0
    pos = v >= 0
    # negative side: blue (0.1, 0.3, 1.0) -> white (1, 1, 1)
    t_neg = v[neg] / vmin  # vmin negative -> t in [0..1]
    out[neg, 0] = 1.0 - t_neg * 0.9
    out[neg, 1] = 1.0 - t_neg * 0.7
    out[neg, 2] = 1.0
    # positive side: white (1, 1, 1) -> red (1.0, 0.2, 0.1)
    t_pos = v[pos] / vmax
    out[pos, 0] = 1.0
    out[pos, 1] = 1.0 - t_pos * 0.8
    out[pos, 2] = 1.0 - t_pos * 0.9
    return out


def make_point_cloud(name, x, y, z, resid):
    n = len(x)
    coords = np.empty((n, 3), dtype=np.float32)
    coords[:, 0] = x * X_SCALE
    coords[:, 1] = y * Y_SCALE
    coords[:, 2] = z * Z_SCALE
    colors = diverging_palette(resid, CUMV_RESID_MIN, CUMV_RESID_MAX)

    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(n)
    mesh.vertices.foreach_set("co", coords.flatten())
    mesh.update()
    col_attr = mesh.color_attributes.new(name="Col", type="FLOAT_COLOR", domain="POINT")
    col_attr.data.foreach_set("color", colors.flatten())
    return bpy.data.objects.new(name, mesh)


def setup_geo_nodes_points(obj, point_radius=POINT_RADIUS):
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
        bgn.inputs["Color"].default_value = (0.05, 0.05, 0.07, 1.0)


def make_camera(name, location, look_at, lens=35):
    cd = bpy.data.cameras.new(name); cd.lens = lens
    cam = bpy.data.objects.new(name, cd)
    cam.location = location
    cam.rotation_euler = (Vector(look_at) - Vector(location)).to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.collection.objects.link(cam)
    return cam


def render_to(filepath, cam):
    bpy.context.scene.camera = cam
    bpy.context.scene.render.filepath = str(filepath)
    bpy.ops.render.render(write_still=True)


def main():
    setup_scene()
    print("[load] cumv data", flush=True)
    x, y, z, resid = csv_to_arrays(DATA_DIR / "cumv_a.csv", DATA_DIR / "cumv_a.npz")
    coll = bpy.data.collections.new("Vis1_CumV")
    bpy.context.scene.collection.children.link(coll)
    print(f"[mesh] cumv ({len(x):,} verts)", flush=True)
    obj = make_point_cloud("Vis1_CumV", x, y, z, resid)
    coll.objects.link(obj)
    setup_geo_nodes_points(obj)

    center = np.array([50, 50, 60])
    cam_iso = make_camera("cam_iso", (200, -150, 200), center)
    cam_front = make_camera("cam_front", (50, -200, 60), center)
    cam_side = make_camera("cam_side", (220, 50, 60), center)
    cam_top = make_camera("cam_top", (50, 50, 240), center)

    blend_path = VIZ_DIR / "cumv_3d.blend"
    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    print(f"[save] {blend_path}", flush=True)

    snaps = [
        ("cumv_iso.png", cam_iso),
        ("cumv_front.png", cam_front),
        ("cumv_side.png", cam_side),
        ("cumv_top.png", cam_top),
    ]
    for fname, cam in snaps:
        out = SNAP_DIR / fname
        print(f"[render] {fname}", flush=True)
        render_to(out, cam)

    print("\n[done]", flush=True)


if __name__ == "__main__":
    main()
