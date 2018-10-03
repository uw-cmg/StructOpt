import logging
import numpy
def compose_structure(Optimizer,individ):
    # Establish individual structure for evaluation.  Piece together regions when necessary.
    logger = logging.getLogger(Optimizer.loggername)
    indiv=individ[0]
    if 'EE' in Optimizer.debug:
        debug = True
    else:
        debug = False
    
    if Optimizer.structure=='Defect':
        indi=indiv.copy()
        bulk=individ.bulki
        nat=indi.get_number_of_atoms()
        if debug: 
            logger.info('Extending defect structure to include bulk len(r1+r2)={0} len(bulk)={1}'.format(nat,len(bulk)))
        csize=bulk.get_cell()                                                                                                         
        totalsol=Atoms(cell=csize, pbc=True)
        totalsol.extend(indi)
        totalsol.extend(bulk)
        for sym,c,m,u in Optimizer.atomlist:
            nc=len([atm for atm in totalsol if atm.symbol==sym])
            STR+='Defect configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n'
    elif Optimizer.structure=='Surface':
        totalsol=Atoms()
        totalsol.extend(indiv)
        nat=indiv.get_number_of_atoms()
        totalsol.extend(individ.bulki)
        if debug:
            logger.info('Extending surface structure to include bulk len(r1+r2)={0} len(bulk)={1}'.format(nat,len(individ.bulki)))
        for sym,c,m,u in Optimizer.atomlist:
            nc=len([atm for atm in totalsol if atm.symbol==sym])
            STR+='Surface-Bulk configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n'
        cell=numpy.maximum.reduce(indiv.get_cell())
        totalsol.set_cell([cell[0],cell[1],500])
        totalsol.set_pbc([True,True,False])
    elif Optimizer.structure=='Cluster':
        totalsol = indiv.copy()
        nat = len(totalsol)
        if debug:
            logger.info('Extending cluster with {0} atoms to center of evaluation box of size {1}'.format(nat,Optimizer.large_box_size))
        origcell = indiv.get_cell()
        totalsol.set_cell([Optimizer.large_box_size,Optimizer.large_box_size,Optimizer.large_box_size])
        totalsol.translate([Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0,Optimizer.large_box_size/2.0])
    elif Optimizer.structure=='Crystal':
        totalsol = indiv.copy()
        nat = len(totalsol)
    else:
        print 'WARNING: In EvalEnergy. Optimizer.structure not recognized'
        logger.warning('Optimizer.structure not recognized')

    return totalsol
