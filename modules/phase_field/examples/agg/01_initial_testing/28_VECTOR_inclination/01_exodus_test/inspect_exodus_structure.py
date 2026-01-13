#!/usr/bin/env python3
from __future__ import annotations

import argparse
from typing import Iterable, Tuple, Set, Dict, Any

try:
    from netCDF4 import Dataset
except ImportError as e:
    raise SystemExit(
        "Missing dependency netCDF4.\n"
        "Install with: mamba install -c conda-forge netcdf4\n"
        "or: pip install netCDF4\n"
    ) from e


def _indent(level: int) -> str:
    return "  " * level


def _safe_filters(var) -> Dict[str, Any] | None:
    # Only works for netCDF4 vars; classic returns AttributeError
    try:
        return var.filters()
    except Exception:
        return None


def _safe_chunking(var) -> Any:
    try:
        return var.chunking()
    except Exception:
        return None


def describe_group(g, level: int = 0, max_attr_values: int = 120) -> None:
    ind = _indent(level)
    print(f"{ind}Group: {g.path if hasattr(g, 'path') else '/'}")

    # Global/group attributes
    attrs = list(g.ncattrs())
    if attrs:
        print(f"{ind}  Attributes ({len(attrs)}):")
        for a in attrs:
            v = getattr(g, a)
            s = repr(v)
            if len(s) > max_attr_values:
                s = s[:max_attr_values] + " ... (truncated)"
            print(f"{ind}    - {a} = {s}")

    # Dimensions
    dims = g.dimensions
    print(f"{ind}  Dimensions ({len(dims)}):")
    for name, d in dims.items():
        # d.isunlimited() exists in netCDF4
        unlimited = ""
        try:
            if d.isunlimited():
                unlimited = " (UNLIMITED)"
        except Exception:
            pass
        print(f"{ind}    - {name}: {len(d)}{unlimited}")

    # Variables
    vars_ = g.variables
    print(f"{ind}  Variables ({len(vars_)}):")
    for name, v in vars_.items():
        dims_str = "(" + ", ".join(v.dimensions) + ")"
        print(f"{ind}    - {name}: dtype={v.dtype} shape={v.shape} dims={dims_str}")

        # Variable attributes (brief)
        vattrs = list(v.ncattrs())
        if vattrs:
            print(f"{ind}        attrs: {', '.join(vattrs)}")

        # NetCDF-4 niceties (chunking/compression)
        filters = _safe_filters(v)
        chunking = _safe_chunking(v)
        if filters or chunking not in (None, "contiguous"):
            if chunking not in (None, "contiguous"):
                print(f"{ind}        chunking: {chunking}")
            if filters:
                print(f"{ind}        filters: {filters}")

    # Recurse into subgroups (netCDF-4 only)
    subgroups = getattr(g, "groups", {})
    if subgroups:
        print(f"{ind}  Subgroups ({len(subgroups)}):")
        for sg_name, sg in subgroups.items():
            print(f"{ind}    - {sg_name}")
        for sg in subgroups.values():
            describe_group(sg, level=level + 1, max_attr_values=max_attr_values)


def collect_schema(ds) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Return (dimension_names, variable_paths, group_paths)
    Variable paths include group prefix, like '/group1/var'
    """
    dim_names: Set[str] = set()
    var_paths: Set[str] = set()
    group_paths: Set[str] = set()

    def walk(g):
        group_paths.add(getattr(g, "path", "/"))
        for dname in g.dimensions.keys():
            dim_names.add(f"{getattr(g, 'path', '/')}/{dname}")
        for vname in g.variables.keys():
            var_paths.add(f"{getattr(g, 'path', '/')}/{vname}")

        for sg in getattr(g, "groups", {}).values():
            walk(sg)

    walk(ds)
    return dim_names, var_paths, group_paths


def diff_schemas(a_name: str, a_ds, b_name: str, b_ds) -> None:
    a_dims, a_vars, a_grps = collect_schema(a_ds)
    b_dims, b_vars, b_grps = collect_schema(b_ds)

    def show_diff(label: str, only_a: Set[str], only_b: Set[str]) -> None:
        if not only_a and not only_b:
            print(f"\n{label}: identical")
            return
        print(f"\n{label}:")
        if only_a:
            print(f"  Only in {a_name} ({len(only_a)}):")
            for x in sorted(only_a):
                print(f"    - {x}")
        if only_b:
            print(f"  Only in {b_name} ({len(only_b)}):")
            for x in sorted(only_b):
                print(f"    - {x}")

    show_diff("Group paths", a_grps - b_grps, b_grps - a_grps)
    show_diff("Dimension names (scoped by group)", a_dims - b_dims, b_dims - a_dims)
    show_diff("Variable names (scoped by group)", a_vars - b_vars, b_vars - a_vars)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Inspect Exodus/NetCDF structure (classic and netcdf4/hdf5) and compare schemas."
    )
    ap.add_argument("--a", default="01_exodus_test_out.e", help="Non-hdf5 Exodus file")
    ap.add_argument("--b", default="01_exodus_test_h5out.e", help="HDF5-backed (netcdf4) Exodus file")
    ap.add_argument("--no-diff", action="store_true", help="Skip schema diff")
    args = ap.parse_args()

    print(f"\n=== OPEN A: {args.a} ===")
    with Dataset(args.a, "r") as dsa:
        describe_group(dsa)

        print(f"\n=== OPEN B: {args.b} ===")
        with Dataset(args.b, "r") as dsb:
            describe_group(dsb)

            if not args.no_diff:
                print("\n=== SCHEMA DIFF (metadata-level) ===")
                diff_schemas("A", dsa, "B", dsb)


if __name__ == "__main__":
    main()
