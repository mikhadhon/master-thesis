from sys import argv
import polyscope as ps
import numpy as np

def main():
    ps.init()

    obj_file = open("/Users/khoa/CLionProjects/master-thesis/" + argv[1])
    vertices = []
    simplices = []
    dim = 0
    simplex_dim = 0

    for line in obj_file:
        data = line.split(" ")

        if data[0] == "v":
            if dim == 0:
                dim = len(data) - 1

            vertex_coords = []
            for i in range(dim):
                vertex_coords.append(data[i + 1])

            vertices.append(vertex_coords)

        elif data[0] == "f":
            if simplex_dim == 0:
                simplex_dim = len(data) - 1

            simplex_indices = []
            for i in range(simplex_dim):
                simplex_indices.append(int(data[i + 1]) - 1)

            simplices.append(simplex_indices)

    np_vertices = np.array(vertices)
    np_simplices = np.array(simplices)

    if simplex_dim == 2:
        ps.register_curve_network("Reconstruction", np_vertices, np_simplices)

    elif simplex_dim == 3:
        ps.register_surface_mesh("Reconstruction", np_vertices, np_simplices)

    ps.show()

if __name__ == "__main__":
    main()
