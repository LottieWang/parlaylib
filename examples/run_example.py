import subprocess
import os

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(f"{CURRENT_DIR}/../log", exist_ok=True)  

GRAPH_DIR="~/data"
graphs = ["Epinions1", "Slashdot", "DBLP"]

# def collect_data(file, keywords):


if __name__=="__main__":
    numKbitParallel=0
    out_file = f"{CURRENT_DIR}/../log/stats_{numKbitParallel}.out"
    for g in graphs:
        in_file = f"{GRAPH_DIR}/{g}_sym.bin"
        cmd = f"{CURRENT_DIR}/../build/examples/BFS_labeling {in_file} | tee -a {out_file}"
        subprocess.call(cmd, shell=True)
    