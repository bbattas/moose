#!/usr/bin/env python3
"""
slurm_create.py

Generates (and optionally submits) SLURM batch scripts for MOOSE
phase-field jobs. Finds .i input files in the current directory
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

# HPG RH9/EL9 module stack
MODULES = "ufrc mkl/2025.1.0 gcc/14.2.0 openmpi/5.0.10 python/3.12 cmake/3.30.5"
MPI_FLAG = "--mpi=pmix_v5"
MPI_EXPORTS = "export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77"

# INL configuration
INL_MODULES = "module load use.moose moose-dev-openmpi/2026.07.30"
INL_GENERAL_MAX_HOURS = 168
INL_SHORT_MAX_HOURS = 6

# =============================================================================
# ARGUMENT PARSER
# =============================================================================

parser = argparse.ArgumentParser(
    description="Generate SLURM scripts for MOOSE phase-field jobs.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

# --- Behavior / Control ---
grp_ctrl = parser.add_argument_group("Behavior / Control")
grp_ctrl.add_argument(
    "--verbose", "-v", action="store_true",
    help="Enable verbose (INFO-level) output.",
)
grp_ctrl.add_argument(
    "--subdirs", "-s", action="store_true",
    help="Run in all direct subdirectories rather than only the current directory.",
)
grp_ctrl.add_argument(
    "--overwrite", "-o", action="store_true",
    help="Overwrite existing .sh scripts.",
)
grp_ctrl.add_argument(
    "--run", "-r", action="store_true",
    help="Submit all generated scripts via sbatch after writing.",
)
grp_ctrl.add_argument(
    "--yes", "-y", action="store_true",
    help="Skip the header preview confirmation prompt (non-interactive mode).",
)
grp_ctrl.add_argument(
    "--preview-full", action="store_true",
    help="Show the full script as the preview instead of just the header.",
)
grp_ctrl.add_argument(
    "--inl", action="store_true",
    help="Use the INL SLURM, module, and execution configuration.",
)

# --- SLURM Header ---
grp_hdr = parser.add_argument_group("SLURM Header")
grp_hdr.add_argument(
    "--nodes", "-n", type=int, default=1,
    help="Number of nodes.",
)
grp_hdr.add_argument(
    "--tasks", type=int, default=None,
    help="Number of MPI tasks per node (default: 30 on HPG, 112 on INL).",
)
grp_hdr.add_argument(
    "--cpus-per-task", type=int, default=1,
    help="Number of CPUs per task.",
)
grp_hdr.add_argument(
    "--mem-per-cpu", type=str, default=None,
    help="Memory per CPU (default: 7GB on HPG, 2200MB on INL).",
)
grp_hdr.add_argument(
    "--full-mem", action="store_true",
    help="INL only: request all node memory using '#SBATCH --mem=0'.",
)
grp_hdr.add_argument(
    "--partition", type=str, default=None,
    help="SLURM partition. INL uses general by default or short with --burst.",
)
grp_hdr.add_argument(
    "--burst", "-b", action="store_true",
    help="Use burst resources. On INL, selects the short partition.",
)
grp_hdr.add_argument(
    "--time", "-t", type=int, default=None,
    help="Wall time in hours. INL defaults: 160 general, 5 short/burst.",
)
grp_hdr.add_argument(
    "--project", choices=("neams", "neup"), default="neams",
    help="INL WCKEY/project.",
)
grp_hdr.add_argument(
    "--array", type=str, default=None,
    help="Relative path to a .txt file for SLURM array submission.",
)

# --- Job Configuration ---
grp_job = parser.add_argument_group("Extras")
grp_job.add_argument(
    "--use-input-name", action="store_true",
    help="Name the SLURM job after the .i filename instead of the directory name.",
)
grp_job.add_argument(
    "--recover", nargs="?", const=True, default=None,
    metavar="CHECKPOINT_PATH",
    help="Add --recover to the MOOSE command. Optionally specify a checkpoint "
         "path (e.g. --recover 02_inc_checkpoint_cp/0052).",
)
grp_job.add_argument(
    "--args", "-c", type=str, default=None,
    help="Extra arguments appended to the MOOSE command line.",
)
grp_job.add_argument(
    "--no-email", action="store_true",
    help="Disable email notifications on begin/end/fail.",
)
grp_job.add_argument(
    "--exclude", action="store_true",
    help=f"Add '#SBATCH --exclude={BAD_NODES}' to the script.",
)
grp_job.add_argument(
    "--mpi", action="store_true",
    help="HPG only: add an OpenMPI collective-communication workaround.",
)

args = parser.parse_args()

# =============================================================================
# CONFIGURATION RESOLUTION / VALIDATION
# =============================================================================

def resolve_configuration():
    """Resolve cluster-specific settings and validate selected options."""
    if not args.inl:
        return {
            "tasks": args.tasks if args.tasks is not None else 30,
            "mem_per_cpu": args.mem_per_cpu if args.mem_per_cpu is not None else "7GB",
            "partition": args.partition,
            "time": args.time if args.time is not None else 72,
        }

    if args.full_mem and args.mem_per_cpu is not None:
        parser.error("--full-mem cannot be combined with --mem-per-cpu.")

    expected_partition = "short" if args.burst else "general"

    if args.partition not in (None, expected_partition):
        parser.error(
            f"INL {'burst' if args.burst else 'regular'} jobs must use "
            f"the '{expected_partition}' partition."
        )

    default_time = 5 if args.burst else 160
    maximum_time = INL_SHORT_MAX_HOURS if args.burst else INL_GENERAL_MAX_HOURS
    requested_time = args.time if args.time is not None else default_time

    if requested_time < 1:
        parser.error("--time must be at least one hour.")

    if requested_time > maximum_time:
        parser.error(
            f"INL {expected_partition} jobs are limited to "
            f"{maximum_time} hours."
        )

    return {
        "tasks": args.tasks if args.tasks is not None else 112,
        "mem_per_cpu": args.mem_per_cpu if args.mem_per_cpu is not None else "2200MB",
        "partition": expected_partition,
        "time": requested_time,
    }


config = resolve_configuration()

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

        logging.info("  Found existing .sh script, skipping.")
        return True

    return False


def check_for_input():
    """
    Return (True, input_name) if at least one .i file is found.

    Warn and return the first file if multiple .i files are present.
    Return (False, None) if no .i file is found.
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

    Raises ValueError if zero or more than one matching file is found.
    """
    dir_lead = "../" if args.subdirs else ""
    path = Path(dir_lead + args.array).expanduser()

    if path.suffix.lower() == ".txt":
        matches = [path] if path.is_file() else []
    else:
        parent = path.parent if str(path.parent) != "" else Path(".")
        matches = list(parent.glob(path.name + "*.txt"))

    if len(matches) == 0:
        raise ValueError(f"  No .txt file found matching: {path}")

    if len(matches) > 1:
        raise ValueError(
            f"  More than one .txt file found matching: {path}\n"
            + "\n".join(f"    {match}" for match in matches)
        )

    full_path = os.path.abspath(matches[0])

    with open(full_path) as array_file:
        nrows = sum(1 for _ in array_file)

    logging.info(f"  Array file : {full_path}")
    logging.info(f"  Array rows : {nrows}")
    return nrows, full_path


def build_recover_flag():
    """Return the --recover command-line fragment."""
    if args.recover is None:
        return ""

    if args.recover is True:
        return " --recover"

    return f" --recover {args.recover}"


def build_script(cwd, input_name):
    """
    Build and return the full SLURM script as a list of lines.

    This function is the single source of truth used by preview and writing.
    """
    lines = []
    array_num = None
    array_path = None

    if args.array is not None:
        array_num, array_path = check_for_array()

    # --- Header ---
    lines.append("#!/bin/bash")
    lines.append("")

    job_name = cwd.rsplit("/", 1)[1] if not args.use_input_name else input_name

    lines.append(f"#SBATCH --job-name={job_name}")
    lines.append(f"#SBATCH --nodes={args.nodes}")
    lines.append(f"#SBATCH --ntasks-per-node={config['tasks']}")
    lines.append(f"#SBATCH --cpus-per-task={args.cpus_per_task}")

    if args.inl and args.full_mem:
        lines.append("#SBATCH --mem=0")
    else:
        lines.append(f"#SBATCH --mem-per-cpu={config['mem_per_cpu']}")

    if not args.inl:
        lines.append("#SBATCH --distribution=cyclic:cyclic")

    if array_num is not None:
        lines.append(f"#SBATCH --array=1-{array_num}")

    lines.append("")

    if config["partition"] not in (None, "None"):
        lines.append(f"#SBATCH --partition={config['partition']}")

    lines.append(f"#SBATCH --time={config['time']}:00:00")

    if array_num is not None:
        lines.append("#SBATCH --output=moose_console_%A_%a.out")
    else:
        lines.append("#SBATCH --output=moose_console_%j.out")

    lines.append("#SBATCH --mail-user=bbattas@ufl.edu")

    if not args.no_email:
        lines.append("#SBATCH --mail-type=BEGIN,END,FAIL")

    if args.inl:
        lines.append(f"#SBATCH --wckey={args.project}")
    else:
        lines.append("#SBATCH --account=michael.tonks")

        if args.burst:
            logging.info("    Burst allocation specified.")
            lines.append("#SBATCH --qos=michael.tonks-b")

    if args.exclude:
        lines.append(f"#SBATCH --exclude={BAD_NODES}")

    if not args.inl:
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

    if args.inl:
        lines.append(INL_MODULES)
    else:
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
            "LINE=\"$(awk 'NF>0 && $1 !~ /^#/' \"${OV_FILE}\" | "
            "sed -n \"${SLURM_ARRAY_TASK_ID}p\")\"",
            "",
            "if [[ -z \"${LINE}\" ]]; then",
            "  echo \"No overrides for task ${SLURM_ARRAY_TASK_ID}. "
            "Check the array size or ${OV_FILE}.\"",
            "  exit 1",
            "fi",
            "",
            "echo \"[$(date)] Task ${SLURM_ARRAY_TASK_ID} overrides: ${LINE}\"",
        ])

    # --- Execution line ---
    recover_str = build_recover_flag()
    extra_args = f" {args.args}" if args.args else ""

    lines.append("")
    lines.append("cd $OUTPUT")

    if args.inl:
        total_tasks = args.nodes * config["tasks"]
        command = (
            f"mpiexec -n {total_tasks} moose-dev-exec $MOOSE "
            f"-i $OUTPUT/{input_name}.i"
        )
    else:
        command = f"srun {MPI_FLAG} $MOOSE -i $OUTPUT/{input_name}.i"

    if array_path is not None:
        command += recover_str if recover_str else ""
        command += " ${LINE}"
        command += extra_args
    else:
        command += recover_str
        command += extra_args

    lines.append(command)

    return lines


def show_preview(cwd, input_name):
    """
    Print the script preview.

    Show the header by default, or the complete script with --preview-full.
    Prompt for confirmation unless --yes is supplied.
    """
    script_lines = build_script(cwd, input_name)

    if args.preview_full:
        preview_lines = script_lines
        label = "Full Script Preview:"
    else:
        last_sbatch = 0

        for index, line in enumerate(script_lines):
            if line.startswith("#SBATCH"):
                last_sbatch = index

        preview_lines = script_lines[:last_sbatch + 1]
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
    """Build and write slurm_<input_name>.sh in the current directory."""
    script_lines = build_script(cwd, input_name)
    slurm_name = f"slurm_{input_name}.sh"

    logging.info(f"    Writing script: {slurm_name}")

    with open(slurm_name, mode="w") as script_file:
        script_file.write("\n".join(script_lines) + "\n")

    return slurm_name


def run_slurm(slurm_name):
    """Submit a SLURM script by name via sbatch."""
    logging.info(f"  Submitting: {slurm_name}")
    subprocess.run(["sbatch", slurm_name], check=False)


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Validate application executable path.
    if not os.path.exists(pf_opt):
        logging.warning(f"WARNING: Cannot find phase_field-opt at: {pf_opt}")
        confirm = input("Would you like to continue anyway? y/[n]: ")

        if "y" not in confirm.lower():
            logging.warning("Exiting.")
            sys.exit()

    if args.overwrite:
        logging.warning("WARNING: Overwriting existing .sh files.")

    if args.run:
        logging.warning("WARNING: Submitting all generated .sh files via sbatch.")

    if args.exclude:
        logging.warning(f"WARNING: Excluding nodes: {BAD_NODES}")

    # Show one preview before processing directories.
    show_preview(os.getcwd(), "[input_name]")

    exclude_prefixes = (".", "00_")

    if args.subdirs:
        logging.info("Cycling through subdirectories.")
        working_dir = os.getcwd()

        directories = sorted(
            directory for directory in os.listdir()
            if not directory.startswith(exclude_prefixes)
            and os.path.isdir(directory)
        )

        logging.info(f"Subdirectories found: {directories}")

        for directory in directories:
            cwd = os.path.join(working_dir, directory)
            os.chdir(cwd)

            logging.info(f"In directory: {cwd}")

            found, input_name = check_for_input()

            if found and not check_for_slurm():
                slurm_name = write_script(cwd, input_name)

                if args.run:
                    run_slurm(slurm_name)

            logging.info("")

        os.chdir(working_dir)

    else:
        logging.info("Running in current directory only.")
        cwd = os.getcwd()

        logging.info(f"In directory: {cwd}")

        found, input_name = check_for_input()

        if found and not check_for_slurm():
            slurm_name = write_script(cwd, input_name)

            if args.run:
                run_slurm(slurm_name)

    logging.info("Done.")


if __name__ == "__main__":
    main()
