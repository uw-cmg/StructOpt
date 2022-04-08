import numpy
def get_cluster_volume(cluster):
    max_pos = numpy.maximum.reduce(cluster.get_positions())
    min_pos = numpy.minimum.reduce(cluster.get_positions())
    diff_pos = [max_pos[i]-min_pos[i] for i in range(3)]
    vol = diff_pos[0]*diff_pos[1]*diff_pos[2]
    return vol
