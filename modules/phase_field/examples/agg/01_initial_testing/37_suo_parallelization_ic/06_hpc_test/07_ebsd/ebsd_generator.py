import datetime

# ── Mesh parameters ────────────────────────────────────────────────────────────
nx, ny       = 100, 100
xmin, xmax   = 0.0, 1000.0
ymin, ymax   = 0.0, 1000.0
zmin, zmax   = 0.0, 0.0
dx           = (xmax - xmin) / nx   # 10.0
dy           = (ymax - ymin) / ny   # 10.0
dz           = 0.0

num_grains   = 15
n_cols_grain = 3   # partition x into 3 strips
n_rows_grain = 5   # partition y into 5 strips

# ── Grain assignment ───────────────────────────────────────────────────────────
def get_grain_id(cx, cy):
    col = min(int((cx - xmin) / (xmax - xmin) * n_cols_grain), n_cols_grain - 1)
    row = min(int((cy - ymin) / (ymax - ymin) * n_rows_grain), n_rows_grain - 1)
    return col * n_rows_grain + row + 1  # 1-indexed, range [1, 15]

# ── Fixed Euler angles ─────────────────────────────────────────────────────────
phi1, PHI, phi2 = 0.0, 0.0, 0.0
phase_id        = 1
symmetry        = 43

# ── Build data rows ────────────────────────────────────────────────────────────
rows = []
for j in range(ny):                # y outer loop
    cy = ymin + (j + 0.5) * dy
    for i in range(nx):            # x inner loop
        cx = xmin + (i + 0.5) * dx
        gid = get_grain_id(cx, cy)
        rows.append((phi1, PHI, phi2, cx, cy, zmin, gid, phase_id, symmetry))

# ── Write file ─────────────────────────────────────────────────────────────────
timestamp = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
filename  = "artificial_ebsd_15grains.txt"

with open(filename, "w") as f:
    f.write(f"# File written from Python Script\n")
    f.write(f"# DateTime: {timestamp}\n")
    f.write(f"# X_STEP: {dx:.6f}\n")
    f.write(f"# Y_STEP: {dy:.6f}\n")
    f.write(f"# Z_STEP: {dz:.6f}\n")
    f.write("#\n")
    f.write(f"# X_MIN: {xmin:.6f}\n")
    f.write(f"# Y_MIN: {ymin:.6f}\n")
    f.write(f"# Z_MIN: {zmin:.6f}\n")
    f.write("#\n")
    f.write(f"# X_MAX: {xmax:.6f}\n")
    f.write(f"# Y_MAX: {ymax:.6f}\n")
    f.write(f"# Z_MAX: {zmax:.6f}\n")
    f.write("#\n")
    f.write(f"# X_DIM: {nx}\n")
    f.write(f"# Y_DIM: {ny}\n")
    f.write(f"# Z_DIM: 0\n")
    f.write("#\n")
    f.write(f"# Phase_1: Primary\n")
    f.write(f"# Symmetry_1: {symmetry}\n")
    f.write(f"# Features_1: {num_grains}\n")
    f.write("#\n")
    f.write(f"# Num_Features: {num_grains} \n")
    f.write("#\n")
    f.write("# phi1 PHI phi2 x y z FeatureId PhaseId Symmetry\n")

    for r in rows:
        f.write(
            f"{r[0]:.6f} {r[1]:.6f} {r[2]:.6f} "
            f"{r[3]:.6f} {r[4]:.6f} {r[5]:.6f} "
            f"{r[6]} {r[7]} {r[8]}\n"
        )

print(f"Wrote {len(rows)} rows to '{filename}'")
print(f"Grain IDs present: {sorted(set(r[6] for r in rows))}")
