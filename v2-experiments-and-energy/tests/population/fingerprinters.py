from structopt.common.population import Population
from structopt.tools.dictionaryobject import DictionaryObject
import structopt


parameters = structopt.setup(DictionaryObject({
    "structure_type": "cluster",
    "generators": {
	"sphere": {"number_of_individuals": 2,
		   "kwargs": {"atomlist": [["Au", 55]],
			      "cell": [20, 20, 20]}}
    },
    "fingerprinters": {
        "all_close_atom_positions": {"probability": 1.0}
    },
}))


def test_all_close_atom_postions():
    parameters.fingerprinters = {
        "keep_best": True,
        "all_close_atom_positions": {"probability": 1.0, "kwargs": {}}
    }
    pop = Population(parameters=parameters)
    pop.initial_number_of_individuals = 1  # Need to override this because this is the value of nkeep that gets passed to the fingerprinter

    pop[1].set_positions(pop[0].get_positions())

    pop.apply_fingerprinters()
    assert len(pop) == 1



def test_diversify_module():
    parameters.fingerprinters = {
        "keep_best": True,
        "diversify_module": {"probability": 1.0,
            "kwargs": {"module": "LAMMPS"}
        }
    }

    pop = Population(parameters=parameters)
    pop.initial_number_of_individuals = 1  # Need to override this because this is the value of nkeep that gets passed to the fingerprinter

    pop[0].LAMMPS = 1.234
    pop[1].LAMMPS = 1.234

    pop.apply_fingerprinters()
    assert len(pop) == 1



if __name__ == "__main__":
    test_all_close_atom_postions()
    test_diversify_module()

