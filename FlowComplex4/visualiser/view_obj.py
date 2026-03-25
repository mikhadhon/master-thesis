from sys import argv
import polyscope as ps
import numpy as np

def stereo_projection(object_4d):
    l = 1.0
    projection_matrix = np.array([
        [l, 0.0, 0.0, 0.0],
        [0.0, l, 0.0, 0.0],
        [0.0, 0.0, l, 0.0]
    ])

    projected = []
    for o in object_4d:
        projected.append((projection_matrix / (l - o[3])) @ o)
    return np.array(projected)

def main():
    ps.init()

    obj_file = open("/Users/khoa/CLionProjects/master-thesis/FlowComplex4/" + argv[1])
    vertices = []
    simplices = []
    dim = 4
    simplex_dim = 3

    for line in obj_file:
        data = line.split(" ")

        if data[0] == "v":
            vertex_coords = []
            for i in range(dim):
                vertex_coords.append(float(data[i + 1]))

            vertices.append(vertex_coords)

        elif data[0] == "f":
            simplex_indices = []
            for i in range(simplex_dim):
                simplex_indices.append(int(data[i + 1]) - 1)

            simplices.append(simplex_indices)

    np_vertices = stereo_projection(np.array(vertices))
    np_simplices = np.array(simplices)

    ps.register_surface_mesh("vertices", np_vertices, np_simplices)

    ps.show()

if __name__ == "__main__":
    main()
