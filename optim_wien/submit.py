"""
Job Submission Support — .machines files, SLURM/PBS templates.

References:
  - WIEN2k User's Guide, Section 3.5 (parallel execution)
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
"""

import os
from dataclasses import dataclass, field


@dataclass
class MachinesConfig:
    granules: list = field(default_factory=lambda: [1])
    processes_per_granule: list = field(default_factory=lambda: [1])
    total_processes: int = 1


def generate_machines(basename: str, nproc: int = 1, output_dir: str = ".") -> str:
    """Generate a case.machines file for WIEN2k parallel execution.

    Format per WIEN2k User's Guide, Section 3.5:
      granularity:processes:host1:host2:...

    Args:
        basename: WIEN2k case name
        nproc: Number of parallel processes
        output_dir: Output directory

    Returns:
        Path to generated .machines file
    """
    granules = max(1, nproc)
    lines = [f"{granules}:{nproc}"]
    for i in range(nproc):
        lines.append(f"1:{i+1}:localhost")
    lines.append("")

    path = os.path.join(output_dir, f"{basename}.machines")
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def generate_machines_advanced(
    basename: str,
    nproc_k: int = 2,
    nproc_lapw0: int = 1,
    nproc_lapw1: int = 4,
    nproc_lapw2: int = 4,
    output_dir: str = ".",
) -> str:
    """Generate an advanced .machines file with per-program granularity.

    WIEN2k allows different parallelization for each program (kgen, lapw0,
    lapw1, lapw2). This fine-grained control is documented in the WIEN2k
    User's Guide for optimal performance on heterogeneous clusters.

    Args:
        basename: WIEN2k case name
        nproc_k: Processes for kgen
        nproc_lapw0: Processes for lapw0
        nproc_lapw1: Processes for lapw1
        nproc_lapw2: Processes for lapw2
        output_dir: Output directory

    Returns:
        Path to generated .machines file
    """
    lines = [
        "############################################################",
        "# WIEN2k .machines file — auto-generated",
        "# Format: granularity:processes:host1:host2:...",
        "# Reference: WIEN2k User's Guide, Section 3.5",
        "############################################################",
        "",
        f"kpar:granularity=1:processes={nproc_k}",
        f"lapw0:granularity=1:processes={nproc_lapw0}",
        f"lapw1:granularity=1:processes={nproc_lapw1}",
        f"lapw2:granularity=1:processes={nproc_lapw2}",
        "",
        ":default",
    ]
    for i in range(max(nproc_k, nproc_lapw0, nproc_lapw1, nproc_lapw2)):
        lines.append(f"1:{i+1}:localhost")
    lines.append("")

    path = os.path.join(output_dir, f"{basename}.machines")
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def generate_slurm_script(
    basename: str,
    job_name: str = "",
    nproc: int = 4,
    walltime: str = "24:00:00",
    memory: str = "8G",
    partition: str = "normal",
    work_dir: str = ".",
    rkmax: float = 7.0,
    numk: int = 1000,
    ecut: int = 6,
    parallel: bool = True,
    email: str = "",
    output_dir: str = ".",
) -> str:
    """Generate a SLURM submission script for WIEN2k.

    Produces a ready-to-use batch script that runs init_lapw and
    run_lapw with the optimized parameters.

    References:
      - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
      - SLURM documentation: https://slurm.schedmd.com/

    Args:
        basename: WIEN2k case name
        job_name: SLURM job name (default: basename)
        nproc: Number of MPI processes
        walltime: Max wall time (HH:MM:SS)
        memory: Memory per node
        partition: SLURM partition
        work_dir: Working directory for the calculation
        rkmax: RKMAX value
        numk: Total k-points
        ecut: Ecut absolute value (Ry)
        parallel: Use -p flag for run_lapw
        email: Email for SLURM notifications
        output_dir: Output directory for script

    Returns:
        Path to generated .sh script
    """
    jn = job_name or basename
    par_flag = " -p" if parallel else ""
    email_line = f"#SBATCH --mail-user={email}\n#SBATCH --mail-type=ALL" if email else ""

    script = f"""#!/bin/bash
#SBATCH --job-name={jn}
#SBATCH --nodes=1
#SBATCH --ntasks={nproc}
#SBATCH --time={walltime}
#SBATCH --mem={memory}
#SBATCH --partition={partition}
#SBATCH --output={jn}_%j.out
{email_line}

# Load WIEN2k environment (adjust path as needed)
# source /path/to/wien2k/SOURCEME.sh

# Set up WIEN2k parallel
export WIEN_PARALLEL=mpi
export SCRATCH=$SLURM_TMPDIR

# Copy input files to scratch
cp -r {work_dir}/* $SCRATCH/
cd $SCRATCH

# Initialize WIEN2k
init_lapw -b -rkmax {rkmax} -numk {numk} -ecut {ecut}

# Run self-consistent field calculation
run_lapw -ec 0.0001 -cc 0.001 -i 80{par_flag}

# Copy results back
cp -r $SCRATCH/* {work_dir}/

echo "WIEN2k calculation completed."
"""
    path = os.path.join(output_dir, f"submit_{basename}_slurm.sh")
    with open(path, "w") as f:
        f.write(script)
    os.chmod(path, 0o755)
    return path


def generate_pbs_script(
    basename: str,
    job_name: str = "",
    nproc: int = 4,
    walltime: str = "24:00:00",
    memory: str = "8gb",
    queue: str = "workq",
    work_dir: str = ".",
    rkmax: float = 7.0,
    numk: int = 1000,
    ecut: int = 6,
    parallel: bool = True,
    email: str = "",
    output_dir: str = ".",
) -> str:
    """Generate a PBS submission script for WIEN2k.

    Produces a ready-to-use batch script compatible with Torque/PBS
    schedulers commonly found on HPC clusters running WIEN2k.

    References:
      - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
      - PBS documentation: https://www.altair.com/pbs-works/

    Args:
        basename: WIEN2k case name
        job_name: PBS job name
        nproc: Number of MPI processes
        walltime: Max wall time (HH:MM:SS)
        memory: Memory per node
        queue: PBS queue name
        work_dir: Working directory
        rkmax: RKMAX value
        numk: Total k-points
        ecut: Ecut absolute value (Ry)
        parallel: Use -p flag for run_lapw
        email: Email for PBS notifications
        output_dir: Output directory for script

    Returns:
        Path to generated .pbs script
    """
    jn = job_name or basename
    par_flag = " -p" if parallel else ""
    email_line = f"#PBS -M {email}\n#PBS -m abe" if email else ""

    script = f"""#!/bin/bash
#PBS -N {jn}
#PBS -l nodes=1:ppn={nproc}
#PBS -l walltime={walltime}
#PBS -l mem={memory}
#PBS -q {queue}
#PBS -j oe
#PBS -o {jn}_$PBS_JOBID.log
{email_line}

# Load WIEN2k environment (adjust path as needed)
# source /path/to/wien2k/SOURCEME.sh

# Set up WIEN2k parallel
export WIEN_PARALLEL=mpi

cd $PBS_O_WORKDIR

# Initialize WIEN2k
init_lapw -b -rkmax {rkmax} -numk {numk} -ecut {ecut}

# Run self-consistent field calculation
run_lapw -ec 0.0001 -cc 0.001 -i 80{par_flag}

echo "WIEN2k calculation completed."
"""
    path = os.path.join(output_dir, f"submit_{basename}_pbs.pbs")
    with open(path, "w") as f:
        f.write(script)
    os.chmod(path, 0o755)
    return path
