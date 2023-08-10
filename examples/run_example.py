import subprocess
import os

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
subprocess.call(f'mkdir {CURRENT_DIR}/../log', shell=True)

GRAPH_DIR="/ssd0/graphs/link"
graphs = ["Epinions1", "Slashdot", "DBLP", "Youtube"]
for g in graphs:
    out_file = f"{CURRENT_DIR}/../log/{g}_lb.out"
    in_file = f"{GRAPH_DIR}/{g}_sym.bin"
    cmd = f"{CURRENT_DIR}/2hop_cover {in_file} | tee {out_file}"