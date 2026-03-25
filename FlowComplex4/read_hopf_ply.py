import numpy as np
from plyfile import PlyData

def load_hopf_vertices(hpath):
    hopf = []
    with open(hpath, 'rb') as f:
        plydata = PlyData.read(f)
        for i in list(range(len(plydata.elements[0]))):
            hopf.append(plydata.elements[0].data[i].tolist()[:4])
    return np.array(hopf)

if __name__ == "__main__":
    filename = "happy_vrip_res4"
    hopf = load_hopf_vertices(filename + ".ply")

    f = open(filename + ".obj", "w")
    for i in range(len(hopf)):
        f.write("v %f %f %f %f\n" % (hopf[i][0], hopf[i][1], hopf[i][2], 0))

    f.close()