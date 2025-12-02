import numpy as np
from ase import io


gs_file = input("请输入 gs 文件路径: ")
es_file = input("请输入 es 文件路径: ")

gs = io.read(gs_file, format='cp2k-restart')
es = io.read(es_file, format='cp2k-restart')

dr = es.get_positions() - gs.get_positions()
m = es.get_masses()
dq_values = np.sqrt(np.dot(m, np.square(dr)).sum())

print(f"dq = {dq_values:.6f} Å AMU^1/2")