import glob
import re
import csv
import os
import numpy as np
from ase import io

def get_prefix_suffix():
	candidate_pairs = [('gs_sp', 'ex_sp'), ('gs-sp', 'ex-sp')]
	pattern = re.compile(r"^(.+?)_\d+_(.+?)\.out$")
	for gs_folder, ex_folder in candidate_pairs:
		if os.path.isdir(gs_folder) and os.path.isdir(ex_folder):
			for f in glob.glob(os.path.join(gs_folder, "*.out")):
				match = pattern.match(os.path.basename(f))
				if match:
					prefix = match.group(1)
					return prefix, gs_folder, ex_folder
	for gs_folder in ['gs_sp', 'gs-sp']:
		if os.path.isdir(gs_folder):
			for f in glob.glob(os.path.join(gs_folder, "*.out")):
				match = pattern.match(os.path.basename(f))
				if match:
					prefix = match.group(1)
					for ex_folder in ['ex_sp', 'ex-sp']:
						if os.path.isdir(ex_folder):
							return prefix, gs_folder, ex_folder
					return prefix, gs_folder, None
	raise ValueError("The .out files following naming rules not found in expected subfolders")

def calc_dq_values(prefix, gs_folder):
	dq_values = []
	ref_path = os.path.join(gs_folder, f'{prefix}_1_gs.inp')
	try:
		gs = io.read(ref_path, format='cp2k-restart')
	except Exception:
		return ["MISSING"] * 13
	dq_values.append("0")
	for i in range(2, 14):
		path_i = os.path.join(gs_folder, f'{prefix}_{i}_gs.inp')
		try:
			es = io.read(path_i, format='cp2k-restart')
			dr = es.get_positions() - gs.get_positions()
			m = es.get_masses()
			dq = np.sqrt(np.dot(m, np.square(dr)).sum())
			dq_values.append(f"{dq:.6f}")
		except Exception:
			dq_values.append("MISSING")
	return dq_values

def extract_energies(prefix, gs_folder, ex_folder):
	pattern = re.compile(
		r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[(?:a\.u\.|hartree)\]\s*:?\s*(-?\d+\.\d+)"
	)
	energies = []

	hartree_ref = None
	ref_out = os.path.join(gs_folder, f"{prefix}_1_gs.out")
	try:
		with open(ref_out, 'r', encoding='utf-8') as f:
			ref_text = f.read()
		ref_match = pattern.search(ref_text)
		if ref_match:
			hartree_ref = float(ref_match.group(1))
	except Exception:
		hartree_ref = None

	for i in range(1, 14):
		gs_file = os.path.join(gs_folder, f"{prefix}_{i}_gs.out")
		try:
			with open(gs_file, 'r', encoding='utf-8') as f:
				text_gs = f.read()
			match_gs = pattern.search(text_gs)
			if match_gs and hartree_ref is not None:
				gs_val = f"{(float(match_gs.group(1)) - hartree_ref) * 27.2114:.4f}"
			else:
				gs_val = "" if match_gs is None else ""
		except FileNotFoundError:
			gs_val = "MISSING"

		if ex_folder:
			ex_file = os.path.join(ex_folder, f"{prefix}_{i}_ex.out")
			try:
				with open(ex_file, 'r', encoding='utf-8') as f:
					text_ex = f.read()
				match_ex = pattern.search(text_ex)
				if match_ex and hartree_ref is not None:
					ex_val = f"{(float(match_ex.group(1)) - hartree_ref) * 27.2114:.4f}"
				else:
					ex_val = "" if match_ex is None else ""
			except FileNotFoundError:
				ex_val = "MISSING"
		else:
			ex_val = "MISSING"

		energies.append([gs_val, ex_val])
	return energies

def main():
	prefix, gs_folder, ex_folder = get_prefix_suffix()
	dq_values = calc_dq_values(prefix, gs_folder)
	energies = extract_energies(prefix, gs_folder, ex_folder)
	rows = []
	for dq, en in zip(dq_values, energies):
		rows.append([dq] + en)
	csv_file = f"{prefix}_CCD.csv"
	with open(csv_file, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(["dq", "gs", "ex"])
		writer.writerows(rows)
	print(f"Results updated in {csv_file}")


if __name__ == "__main__":
	main()