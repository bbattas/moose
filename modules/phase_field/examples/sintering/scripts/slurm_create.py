#!/usr/bin/env python3
"""
slurm_create.py
Generates (and optionally submits) SLURM batch scripts for MOOSE
phase-field jobs on HPG. Finds .i input files in the current directory
or one level of subdirectories and writes a corresponding .sh script.
"""

import os
import glob
import sys
import argparse
import logging
import subprocess
from pathlib import Path

# =============================================================================
# EDITABLE CONSTANTS
# =============================================================================

# Update this list when nodes are known to be problematic.
# Used by --exclude to add '#SBATCH --exclude=...' to the script.
BAD_NODES = "c0606a-s18,c0608a-s[6,22,24],c0609a-s10,c0610a-s10"

# RH9/EL9 module stack (only supported environment)
MODULES     = "ufrc mkl/2025.1.0 gcc/14.2.0 openmpi/5.0.10 python/3.12 cmake/3.30.5"
MPI_FLAG    = "--mpi=pmix_v5"
MPI_EXPORTS = "export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77"

# =============================================================================
# ARGUMENT PARSER
# =============================================================================

parser = argparse.ArgumentParser(
    description="Generate SLURM scripts for MOOSE phase-field jobs on HPG.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

# --- Behavior / Control ---
grp_ctrl = parser.add_argument_group("Behavior / Control")
grp_ctrl.add_argument("--verbose",      "-v", action="store_true",
    help="Enable verbose (INFO-level) output.")
grp_ctrl.add_argument("--subdirs",      "-s", action="store_true",
    help="Run in all direct subdirectories rather than only the current directory.")
grp_ctrl.add_argument("--overwrite",    "-o", action="store_true",
    help="Overwrite existing .sh scripts.")
grp_ctrl.add_argument("--run",          "-r", action="store_true",
    help="Submit all generated scripts via sbatch after writing.")
grp_ctrl.add_argument("--yes",          "-y", action="store_true",
    help="Skip the header preview confirmation prompt (non-interactive mode).")
grp_ctrl.add_argument("--preview-full",       action="store_true",
    help="Show the full script as the preview instead of just the header.")

# --- SLURM Header ---
grp_hdr = parser.add_argument_group("SLURM Header")
grp_hdr.add_argument("--nodes",         "-n", type=int,   default=1,
    help="Number of nodes.")
grp_hdr.add_argument("--tasks",               type=int,   default=30,
    help="Number of MPI tasks per node (ntasks-per-node).")
grp_hdr.add_argument("--cpus-per-task",       type=int,   default=1,
    help="Number of CPUs per task.")
grp_hdr.add_argument("--mem-per-cpu",         type=str,   default="7GB",
    help="Memory per CPU.")
grp_hdr.add_argument("--partition",           type=str,   default=None,
    help="SLURM partition. Omitted from script if not specified.")
grp_hdr.add_argument("--burst",         "-b", action="store_true",
    help="Submit as a burst job (adds --qos=michael.tonks-b).")
grp_hdr.add_argument("--time",          "-t", type=int,   default=72,
    help="Wall time in hours. Burst cap: 96h. Regular cap: 744h (31 days).")
grp_hdr.add_argument("--array",               type=str,   default=None,
    help="Relative path to a .txt file for SLURM array submission.")

# --- Job Configuration ---
grp_job = parser.add_argument_group("Extras")
grp_job.add_argument("--use-input-name",      action="store_true",
    help="Name the SLURM job after the .i filename instead of the directory name.")
grp_job.add_argument("--recover",             nargs="?",  const=True, default=None,
    metavar="CHECKPOINT_PATH",
    help="Add --recover to the MOOSE command. Optionally specify a checkpoint path "
         "(e.g. --recover 02_inc_checkpoint_cp/0052).")
grp_job.add_argument("--args",          "-c", type=str,   default=None,
    help="Extra arguments appended to the MOOSE command line.")
grp_job.add_argument("--no-email",            action="store_true",
    help="Disable email notifications on begin/end/fail.")
grp_job.add_argument("--exclude",             action="store_true",
    help=f"Add '#SBATCH --exclude={BAD_NODES}' to the script.")
grp_job.add_argument("--mpi",                 action="store_true",
    help="Add 'export OMPI_MCA_coll_hcoll_enable=0' as an OpenMPI error workaround.")

# # --- Environment ---
# grp_env = parser.add_argument_group("Environment")
# grp_env.add_argument("--mpi",                 action="store_true",
#     help="Add 'export OMPI_MCA_coll_hcoll_enable=0' as an OpenMPI error workaround.")

args = parser.parse_args()

# =============================================================================
# LOGGING SETUP
# =============================================================================

logging.basicConfig(
    level=logging.INFO if args.verbose else logging.WARNING,
    format="%(message)s",
)

# =============================================================================
# PATHING
# =============================================================================

script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
pf_opt = script_directory.rsplit("/", 3)[0] + "/phase_field-opt"

# =============================================================================
# FUNCTIONS
# =============================================================================

def check_for_slurm():
    """Return True if a .sh script already exists and --overwrite is not set."""
    if glob.glob("*.sh"):
        if args.overwrite:
            logging.info("  Overwriting existing .sh script.")
            return False
        else:
            logging.info("  Found existing .sh script, skipping.")
            return True
    return False


def check_for_input():
    """
    Return (True, inputName) if exactly one .i file is found.
    Warns and returns the first found if multiple exist.
    Returns (False, None) if none found.
    """
    matches = glob.glob("*.i")
    if not matches:
        logging.info("  No .i input file found.")
        return False, None
    if len(matches) > 1:
        logging.warning(
            f"  WARNING: Multiple .i files found {matches}. Using: {matches[0]}"
        )
    input_name = matches[0].rsplit(".", 1)[0]
    logging.info(f"  Input file: {matches[0]}")
    return True, input_name


def check_for_array():
    """
    Resolve the --array path to a .txt file and return (num_rows, full_path).
    Raises ValueError if zero or more than one match is found.
    """
    dir_lead = "../" if args.subdirs else ""
    p = Path(dir_lead + args.array).expanduser()

    if p.suffix.lower() == ".txt":
        matches = [p] if p.is_file() else []
    else:
        parent = p.parent if str(p.parent) != "" else Path(".")
        matches = list(parent.glob(p.name + "*.txt"))

    if len(matches) == 0:
        raise ValueError(f"  No .txt file found matching: {p}")
    if len(matches) > 1:
        raise ValueError(f"  More than one .txt file found matching: {p}\n"
                         + "\n".join(f"    {f}" for f in matches))

    full_path = os.path.abspath(matches[0])
    with open(full_path) as f:
        nrows = sum(1 for _ in f)
    logging.info(f"  Array file : {full_path}")
    logging.info(f"  Array rows : {nrows}")
    return nrows, full_path


def build_recover_flag():
    """Return the --recover string fragment based on args.recover."""
    if args.recover is None:
        return ""
    if args.recover is True:
        return " --recover"
    # A path string was provided
    return f" --recover {args.recover}"


def build_script(cwd, input_name):
    """
    Build and return the full SLURM script as a list of strings.
    This is the single source of truth for script content,
    used by both the preview and the file writer.
    """
    lines = []
    array_num  = None
    array_path = None

    if args.array is not None:
        array_num, array_path = check_for_array()

    # --- Header ---
    lines.append("#!/bin/bash")
    lines.append("")

    job_name = cwd.rsplit("/", 1)[1] if not args.use_input_name else input_name
    lines.append(f"#SBATCH --job-name={job_name}")
    lines.append(f"#SBATCH --nodes={args.nodes}")
    lines.append(f"#SBATCH --ntasks-per-node={args.tasks}")
    lines.append(f"#SBATCH --cpus-per-task={args.cpus_per_task}")
    lines.append(f"#SBATCH --mem-per-cpu={args.mem_per_cpu}")
    lines.append( "#SBATCH --distribution=cyclic:cyclic")

    if array_num is not None:
        lines.append(f"#SBATCH --array=1-{array_num}")

    lines.append("")

    if args.partition not in (None, "None"):
        lines.append(f"#SBATCH --partition={args.partition}")

    lines.append(f"#SBATCH --time={args.time}:00:00")

    if array_num is not None:
        lines.append("#SBATCH --output=moose_console_%A_%a.out")
    else:
        lines.append("#SBATCH --output=moose_console_%j.out")

    lines.append("#SBATCH --mail-user=bbattas@ufl.edu")
    if not args.no_email:
        lines.append("#SBATCH --mail-type=BEGIN,END,FAIL")

    lines.append("#SBATCH --account=michael.tonks")
    if args.burst:
        logging.info("    Burst allocation specified.")
        lines.append("#SBATCH --qos=michael.tonks-b")

    if args.exclude:
        lines.append(f"#SBATCH --exclude={BAD_NODES}")

    lines.append("#SBATCH --constraint=el9")

    # --- Environment setup ---
    lines.append("")
    lines.append("echo ${SLURM_JOB_NODELIST}")
    lines.append("")
    lines.append(f"MOOSE={pf_opt}")
    lines.append(f"OUTPUT={cwd}")
    if array_path is not None:
        lines.append(f"OV_FILE={array_path}")

    lines.append("")
    lines.append(MPI_EXPORTS)
    lines.append("module purge")
    lines.append(f"module load {MODULES}")

    if args.mpi:
        lines.append("")
        lines.append("export OMPI_MCA_coll_hcoll_enable=0")

    # --- Array task parsing block ---
    if array_path is not None:
        lines.extend([
            "",
            "LINE=\"$(awk 'NF>0 && $1 !~ /^#/' \"${OV_FILE}\" | sed -n \"${SLURM_ARRAY_TASK_ID}p\")\"",
            "",
            "if [[ -z \"${LINE}\" ]]; then",
            "  echo \"No overrides for task ${SLURM_ARRAY_TASK_ID}. Check the array size or ${OV_FILE}.\"",
            "  exit 1",
            "fi",
            "",
            "echo \"[$(date)] Task ${SLURM_ARRAY_TASK_ID} overrides: ${LINE}\"",
            "echo \"[$(date)] Running: srun ${MOOSE} -i INPUT.i ${LINE}\"",
        ])

    # --- Execution line ---
    recover_str = build_recover_flag()
    extra_args  = f" {args.args}" if args.args else ""
    lines.append("")
    lines.append("cd $OUTPUT")

    if array_path is not None:
        lines.append(
            f"srun {MPI_FLAG} $MOOSE -i $OUTPUT/{input_name}.i"
            + (f"{recover_str} ${{LINE}}" if recover_str else " ${LINE}")
            + extra_args
        )
    else:
        lines.append(
            f"srun {MPI_FLAG} $MOOSE -i $OUTPUT/{input_name}.i"
            + recover_str
            + extra_args
        )

    return lines


def show_preview(cwd, input_name):
    """
    Print the script preview (header only by default, full script with --preview-full).
    Prompts for confirmation unless --yes is set.
    """
    script_lines = build_script(cwd, input_name)

    if args.preview_full:
        preview_lines = script_lines
        label = "Full Script Preview:"
    else:
        # Header = everything up to and including the last #SBATCH line
        last_sbatch = 0
        for i, line in enumerate(script_lines):
            if line.startswith("#SBATCH"):
                last_sbatch = i
        preview_lines = script_lines[: last_sbatch + 1]
        label = "Header Preview:"

    logging.warning("")
    logging.warning(f"\033[1m{label}\033[0m")
    logging.warning("")
    for line in preview_lines:
        logging.warning(line)
    logging.warning("")

    if not args.yes:
        confirm = input("Continue with this script? [y]/n: ")
        if "n" in confirm.lower():
            logging.warning("Exiting.")
            sys.exit()


def write_script(cwd, input_name):
    """Build the script and write it to slurm_<inputName>.sh in the current directory."""
    script_lines = build_script(cwd, input_name)
    slurm_name = f"slurm_{input_name}.sh"
    logging.info(f"    Writing script: {slurm_name}")
    with open(slurm_name, mode="w") as f:
        f.write("\n".join(script_lines))
    return slurm_name


def run_slurm(slurm_name):
    """Submit a SLURM script by name via sbatch."""
    logging.info(f"  Submitting: {slurm_name}")
    subprocess.run(["sbatch", slurm_name])


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Validate executable path
    if not os.path.exists(pf_opt):
        logging.warning(f"WARNING: Cannot find phase_field-opt at: {pf_opt}")
        confirm = input("Would you like to continue anyway? y/[n]: ")
        if "y" not in confirm.lower():
            logging.warning("Exiting.")
            sys.exit()

    # Startup warnings
    if args.overwrite:
        logging.warning("WARNING: Overwriting existing .sh files.")
    if args.run:
        logging.warning("WARNING: Submitting all generated .sh files via sbatch.")
    if args.exclude:
        logging.warning(f"WARNING: Excluding nodes: {BAD_NODES}")

    # Show preview once before iterating (uses cwd as a stand-in for job name preview)
    show_preview(os.getcwd(), "[input_name]")

    # --- Directory iteration ---
    exclude_prefixes = (".", "00_")

    if args.subdirs:
        logging.info("Cycling through subdirectories.")
        working_dir = os.getcwd()
        dirs = sorted(
            d for d in os.listdir()
            if not d.startswith(exclude_prefixes) and os.path.isdir(d)
        )
        logging.info(f"Subdirectories found: {dirs}")

        for d in dirs:
            cwd = os.path.join(working_dir, d)
            os.chdir(cwd)
            logging.info(f"In directory: {cwd}")

            found, input_name = check_for_input()
            if found:
                if not check_for_slurm():
                    slurm_name = write_script(cwd, input_name)
                    if args.run:
                        run_slurm(slurm_name)
            logging.info("")

        os.chdir(working_dir)  # Return to original directory when done

    else:
        logging.info("Running in current directory only.")
        cwd = os.getcwd()
        logging.info(f"In directory: {cwd}")

        found, input_name = check_for_input()
        if found:
            if not check_for_slurm():
                slurm_name = write_script(cwd, input_name)
                if args.run:
                    run_slurm(slurm_name)

    logging.info("Done.")


if __name__ == "__main__":
    main()
