#!/usr/bin/env python3

# python3 mpi_memory.py --label test_04 -- mpiexec -n 4 /Users/bbattas/projects/moose/modules/phase_field/phase_field-opt -i 04_4cpu_main.i

import argparse
import csv
import os
import subprocess
import sys
import time

import psutil


parser = argparse.ArgumentParser(description='Memory tracker, example: python3 mpi_memory.py --label test_04 -- mpiexec -n 4 /Users/bbattas/projects/moose/modules/phase_field/phase_field-opt -i 04_4cpu_main.i')
parser.add_argument("--label", required=True)
parser.add_argument("--interval", type=float, default=0.01)
parser.add_argument("--executable", default="phase_field-opt")
parser.add_argument("command", nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.command and args.command[0] == "--":
    args.command = args.command[1:]

if not args.command:
    parser.error("Provide the mpiexec command after --")

csv_name = f"{args.label}_memory.csv"
launcher = subprocess.Popen(args.command)

known_processes = {}
peak_by_pid = {}
peak_total = 0
start = time.monotonic()

with open(csv_name, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow([
        "elapsed_seconds",
        "pid",
        "rss_mib",
        "simultaneous_total_rss_mib",
    ])

    while True:
        # Discover MPI processes launched under mpiexec.
        try:
            descendants = psutil.Process(launcher.pid).children(recursive=True)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            descendants = []

        for process in descendants:
            try:
                command = process.cmdline()
                if command and os.path.basename(command[0]) == args.executable:
                    known_processes[process.pid] = process
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

        live = []

        for pid, process in list(known_processes.items()):
            try:
                if process.is_running() and process.status() != psutil.STATUS_ZOMBIE:
                    rss = process.memory_info().rss
                    live.append((pid, rss))
                    peak_by_pid[pid] = max(peak_by_pid.get(pid, 0), rss)
                else:
                    del known_processes[pid]
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                known_processes.pop(pid, None)

        elapsed = time.monotonic() - start
        total_rss = sum(rss for _, rss in live)
        peak_total = max(peak_total, total_rss)

        for pid, rss in live:
            writer.writerow([
                f"{elapsed:.4f}",
                pid,
                f"{rss / 1024**2:.3f}",
                f"{total_rss / 1024**2:.3f}",
            ])

        csv_file.flush()

        if launcher.poll() is not None and not live:
            break

        time.sleep(args.interval)

return_code = launcher.wait()

print(f"\nMemory results for {args.label}")
print(f"CSV file: {csv_name}")
print(f"Peak simultaneous total RSS: {peak_total / 1024**2:.2f} MiB")

for pid, peak in sorted(peak_by_pid.items()):
    print(f"PID {pid}: peak RSS {peak / 1024**2:.2f} MiB")

sys.exit(return_code)
