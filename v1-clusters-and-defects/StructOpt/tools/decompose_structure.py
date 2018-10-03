import numpy
from ase import Atom, Atoms
from StructOpt.tools.find_defects import find_defects
from StructOpt.tools.check_cell_type import check_cell_type
from StructOpt.tools.get_cluster_volume import get_cluster_volume
def decompose_structure(Optimizer, totalsol, individ):
   # Separate structures into distinct pieces
    if 'EE' in Optimizer.debug:                        
        debug = True                            
    else:                                                
        debug = False
    
    if Optimizer.structure=='Defect':
        if Optimizer.fixed_region==True or Optimizer.finddefects==False:
            if debug:
                logger.info('Identifying atoms in defect structure based on ID')
            individ[0]=totalsol[0:nat]
            bul=totalsol[(nat):len(totalsol)]
            individ[0].set_cell(csize)
        else:
            if debug:
                logger.info('Applying find defects scheme to identify R1 and R2 for Defect')
            if 'FD' in Optimizer.debug:
                outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=Optimizer.debugfile)
            else:
                outt=find_defects(totalsol,Optimizer.solidbulk,Optimizer.sf,atomlistcheck=Optimizer.atomlist,trackvacs=Optimizer.trackvacs,trackswaps=Optimizer.trackswaps,debug=False)
            individ[0]=outt[0]
            bul=outt[1]
            individ.vacancies = outt[2]
            individ.swaps = outt[3]
            STR += outt[4]
        indiv=individ[0]
    elif Optimizer.structure=='Surface':
        if debug:
            logger.info('Finding surface top layer')
        top,bul=find_top_layer(totalsol,Optimizer.surftopthick)
        indiv=top.copy()
        individ[0]=top.copy()
        bul = Atoms()
    elif Optimizer.structure=='Crystal':
        if debug:
            logger.info('Checking crystal cell type')
        celltype = check_cell_type(totalsol)
        STR+='Cell structure = {0}\n'.format(celltype)
        bul = Atoms()
        individ[0] = totalsol.copy()
    elif Optimizer.structure=='Cluster':
        individ.volume = get_cluster_volume(totalsol)
        bul = Atoms()
        if debug:
            logger.info('Translating cluster back to smaller box size location')
        totalsol.translate([-Optimizer.large_box_size/2.0,-Optimizer.large_box_size/2.0,-Optimizer.large_box_size/2.0])
        totalsol.set_cell(individ[0].get_cell())
        individ[0] = totalsol.copy()

    return individ, bul

