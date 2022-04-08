import numpy as np
from structopt.cluster import Cluster
from structopt.common.individual.fingerprinters.all_close_atom_positions import all_close_atom_positions


def test_all_close_atom_positions():
    a1 = Cluster(id=0, generator_parameters={
            "sphere": {"atomlist": [["Au", 100]],
                "cell": [100.0, 100.0, 100.0],
                "fill_factor": 0.7
                }
            }
    )
    a2 = a1.copy()
    assert all_close_atom_positions(a1, a2)

    shuffle = np.arange(0, len(a2.positions), 1)
    np.random.shuffle(shuffle)
    a2.set_positions(a2.positions[shuffle])
    assert all_close_atom_positions(a1, a2)

    a2.set_positions(a2.positions + np.random.normal(0, 0.01, len(a2.positions)*3).reshape((len(a2.positions),3)))
    assert all_close_atom_positions(a1, a2, rtol=0.001, atol=0.1)


if __name__ == "__main__":
    test_all_close_atom_positions()

