from StructOpt.Optimizer_new import Optimizer

optimizer = Optimizer(relaxation_method="LAMMPS", modules=["FEMSIM", "Random"], weights=[1.0, 2.0])

for i in range(10):
    optimizer.step()
    print('Step {}:'.format(optimizer.step_number))
    for ind in optimizer.individuals:
        print('  Fitness for individual {} = {}'.format(ind, ind.fitness))

