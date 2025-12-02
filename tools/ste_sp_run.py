#!/usr/bin/env python3
import subprocess
import glob
import os
from ase import io
import argparse

def auto_detect_prefix_suffix():
    g_files = glob.glob("*_g_*.inp")
    if not g_files:
        raise FileNotFoundError("The * g *.inp file was not found in the current directory.")
    fname = os.path.splitext(g_files[0])[0]
    parts = fname.split('_g_')
    if len(parts) != 2:
        raise ValueError(f"The file name format is incorrect (*_g_*): {g_files[0]}")
    prefix = parts[0]
    suffix = parts[1]
    return prefix, suffix

def ste_interp(prefix, suffix):
    nimgs = 7
    total_imgs = 13
    initial = io.read(f'{prefix}_g_{suffix}.inp', format='cp2k-restart')
    final = io.read(f'{prefix}_ste_{suffix}.inp', format='cp2k-restart')

    images = [initial]
    images += [initial.copy() for i in range(nimgs - 2)]
    images += [final]
    images += [final.copy() for i in range(total_imgs - nimgs)]

    pos0 = initial.get_positions()
    pos1 = final.get_positions()
    dist = (pos1 - pos0) / (nimgs - 1)

    for i in range(total_imgs):
        images[i].set_positions(pos0 + i * dist)
        
        with open(f'{prefix}_g_{suffix}.inp', 'r') as template_file:
            template_lines = template_file.readlines()
        
        new_lines = []
        in_coord_block = False
        for line in template_lines:
            if '&COORD' in line:
                new_lines.append(line)
                for atom in images[i]:
                    new_lines.append(f'{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n')
                in_coord_block = True
            elif in_coord_block:
                if '&END COORD' in line:
                    new_lines.append(line)
                    in_coord_block = False
            elif f'PROJECT {prefix}_g_{suffix}' in line:
                new_lines.append(line.replace(f'{prefix}_g_{suffix}', f'{prefix}_{i+1}_{suffix}'))
            elif f'WFN_RESTART_FILE_NAME {prefix}_g_{suffix}-RESTART.wfn' in line:
                if i == 0:
                    new_lines.append(line)
                else:
                    new_lines.append(line.replace(f'{prefix}_g_{suffix}-RESTART.wfn', f'{prefix}_{i}_{suffix}-RESTART.wfn'))
            else:
                new_lines.append(line)
        
        with open(f'{prefix}_{i+1}_{suffix}.inp', 'w') as new_file:
            new_file.writelines(new_lines)

def get_user_name():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-name', '-n', dest='user_name', help='User name for job name', type=str)
    args, _ = parser.parse_known_args()
    user_name = args.user_name
    if not user_name or not user_name.strip():
        try:
            user_name = input("user name please: ").strip()
            while not user_name:
                user_name = input("user name can not be empty, input again: ").strip()
        except EOFError:
            raise SystemExit("user name not be provided, program exited.")
    return user_name.strip().replace(' ', '_')
def main():
    prefix, suffix = None, None
    prefix, suffix = auto_detect_prefix_suffix()
    ste_interp(prefix, suffix)

    exe_name = "cp2k.psmp" if suffix == "gs" else "cp2kSTE.psmp"
    user_name = get_user_name()
    sbatch_script = f"""#!/bin/bash
#SBATCH -J {user_name}-{prefix}_{suffix}
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH -o {prefix}_{suffix}.log
#SBATCH -e {prefix}_{suffix}.err

export cp2kroot=/home/think/app/cp2k-2024.1
source $cp2kroot/tools/toolchain/install/setup
export PATH=$cp2kroot/exe/local:$PATH
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=2


for i in $(seq 1 13)
do
    input_file="{prefix}_${{i}}_{suffix}.inp"
    output_file="{prefix}_${{i}}_{suffix}.out"
    echo "task running $i/13"
    start_time=$(date +%s)
    mpirun -np 64 -map-by ppr:4:L3cache:pe=2 {exe_name} $input_file > $output_file 2>&1
    status=$?
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    if [ $status -eq 0 ]; then
        echo "task $i finished, take $elapsed s"
    else
        echo "task $i failed, exit code $status"
    fi
    echo ""
done

echo "All tasks completed."
"""

    script_path = f"{prefix}_{suffix}.sh"
    with open(script_path, "w") as script_file:
        script_file.write(sbatch_script)

    subprocess.run(["sbatch", script_path])

if __name__ == "__main__":
    print("""A script for interpolation and single-point calculation of the STE state.
{prefix}_g/ste_{suffix}.inp files are required.
Author: Ning.Ding  Email: ningdingac@163.com"""
)
    main()