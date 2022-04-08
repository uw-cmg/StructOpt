import math
from mpi4py import MPI
from importlib import import_module

def parallel_mpi4py(individ,module,relax):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()

        if rank==0:
            ntimes=int(math.ceil(float(len(individ))/float(comm.Get_size())))
            nadd=int(ntimes*comm.Get_size()-len(individ))
            maplist=[[] for n in range(ntimes)]
            strt=0
            for i in range(len(maplist)):
                maplist[i]=[indi for indi in individ[strt:comm.Get_size()+strt]]
                strt+=comm.Get_size()
            for i in range(nadd):
                maplist[len(maplist)-1].append(None)
        else:
            ntimes=None
        ntimes = comm.bcast(ntimes,root=0)
        outs=[]
        for i in range(ntimes):
            if rank==0:
                one=maplist[i]
            else:
                one=None
            ind = comm.scatter(one,root=0)
            
            if ind == None:
                rank = MPI.COMM_WORLD.Get_rank()
                stro = 'Evaluated none individual on {0}\n'.format(rank)
                if relax:
                    out = (0, ind, stro)
                else:
                    out = (0, stro)
            else:
                mod = import_module('StructOpt.fitness.%s_eval'%module)
                klass = getattr(mod, '%s_eval'%module)
                out = klass().evaluate_indiv(ind, rank, relax)

            outt = comm.gather(out,root=0)
            if rank == 0:
                outs.extend(outt)
        return outs
