"""A radial distribution function observer."""
from __future__ import print_function
 
import numpy
import asap3
from asap3.Internal.Subject import Subject
from asap3.Internal.ListOfElements import ListOfElements
from asap3.mpi import world
try:
    import cPickle as pickle
except ImportError:
    import pickle
import sys
import numbers
 
class RadialDistributionFunction(Subject):
    def __init__(self, atoms, rMax, nBins, groups=None, interval=1, average=1,
                 autoclear=False, verbose=False):
        """Create a RadialDistributionFunction observer.
 
        Arguments:
        
        atoms: The atoms being observed.
        
        rMax: The maximum distance in the RDF.
        
        nBins: The number of bins in the histograms.
        
        groups (optional): A non-negative integer per atom, used to
        classify atoms into groups, the RDF is calculated for each
        group.
 
        interval (optional, default=1): How often should the RDF be
        calculated.  DEPRECATED: Use interval argument when attaching
        to dynamics instead.
 
        average (optional, default=1): How many times should the RDF
        be calculated and averaged before any observer of this object
        is notified and/or the RDF is saved to a file.
 
        autoclear (optional, default=False): Should the RDF be cleared
        after the RDF has been processed by observers and/or written
        to a file?  The default is to continue to accumulate data.
        """
        Subject.__init__(self)
        self.atoms = atoms
        self.natoms = atoms.get_number_of_atoms()
        self.globalrdf = None
        self.rMax = rMax * 1.0  # Guard against integer divisions.
        self.nBins = nBins
        self.dr = self.rMax / nBins
        self.interval = interval
        self.average = average
        self.autoclear = autoclear
        self.verbose = verbose
        self.autosave = False
        self.n1 = 0  # Counter associated with interval
        self.n2 = 0  # Counter associated with average
        self.clearnext = False
        self.countRDF = 0
        self.listOfElements = ListOfElements(atoms)
        
        if groups is None:
            self.groups = numpy.zeros(len(atoms), numpy.int32)
            self.ngroups = 1
        else:
            if groups.shape != (len(atoms),):
                raise ValueError("groups must be an integer per atom")
            if min(groups) < 0 or max(groups) >= len(atoms):
                raise ValueError("groups array is unreasonable")
            self.groups = groups.astype(numpy.int32)
            self.ngroups = max(self.groups) + 1
 
    def update(self, atoms=None):
        """Calculate the RDF of the atoms.
 
        Make an RDF calculation now (or if interval=n was specified,
        every n time this method is called.
 
        If average = 1, then any observers of this object is notified.
 
        If average > 1, the calculated RDFs are accumulated, and once
        a sufficient number of RDFs have been accumulated, the
        observers are notified and the averaged RDFs are ready.
        """
        if atoms is not None:
            self.atoms = atoms
 
        self.n1 += 1
        if self.n1 >= self.interval:
            if self.clearnext:
                self.clear()
                self.clearnext = False
            self.n1 = 0
            self.do_rdf()
            self.n2 += 1
            if self.n2 >= self.average:
                self.clearnext = self.autoclear
                self.call_observers()  # I have data ready
                if self.autosave:
                    self.save()
                self.n2 = 0
 
    def do_rdf(self):
        """Do the actual RDF calculation.  Do not call directly."""
        if self.verbose:
            print("Calculating RDF", file=sys.stderr)
        assert self.natoms == self.atoms.get_number_of_atoms()
        
        # Check if we need to replace the atoms with a copy.
        atoms = self.atoms
        if world.size > 1:
            try:
                calc = atoms.get_calculator()
            except AttributeError:
                calc = None
            if calc is not None:
                cutoff = calc.get_cutoff()
                if cutoff < self.rMax:
                    if self.verbose:
                        print("Cutoff too long for calculator - copying atoms", file=sys.stderr)
                    atoms = asap3.MakeParallelAtoms(atoms, atoms.nCells)
            
        grdf, rdfs, cnts = asap3._asap.RawRDF(atoms, self.rMax, self.nBins,
                                              self.groups, self.ngroups,
                                              self.listOfElements)
        # Global communication in parallel simulations
        if world.size > 1:
            if self.verbose:
                print("Collecting data", file=sys.stderr)
            world.sum(grdf)
            for i in range(len(rdfs)):
                keys = rdfs[i].keys()
                for key in sorted(keys):
                    world.sum(rdfs[i][key])
            for cnt in cnts:
                keys = cnt.keys()
                for key in sorted(keys):
                    cnt[key] = world.sum(cnt[key])
                
        # Sanity check on counts
        c = 0
        for cnt in cnts: # Loop over groups
            for t in cnt.keys(): # Loop over elements
                c += cnt[t]
        assert c == self.natoms
        # Now store the data
        if self.globalrdf is None:
            self.globalrdf = grdf
            self.rdfs = rdfs
            self.atomcounts = cnts
            # Ensure that all elements are represented
            for i in range(len(rdfs)): # Loop over groups
                for e in self.listOfElements: # Loop over elements
                    if e not in self.atomcounts[i]:
                        self.atomcounts[i][e] = 0
        else:
            self.globalrdf += grdf
            for i in range(len(rdfs)): # Loop over groups
                for t in rdfs[i].keys(): # Loop over pairs
                    self.rdfs[i][t] += rdfs[i][t]
                for t in cnts[i].keys(): # Loop over elements
                    self.atomcounts[i][t] += cnts[i][t]
        self.countRDF += 1
        self.volume = atoms.get_volume()
        if self.verbose:
            print("RDF done!", file=sys.stderr)
 
    def clear(self):
        """Clear the accumulated RDFs.  Called automatically."""
        if self.verbose:
            print("Clearing RDF", file=sys.stderr)
        self.globalrdf = self.rdfs = self.atomcounts = None
        self.countRDF = 0
        
    def get_rdf(self, groups=None, elements=None, c='volume'):
        """Get an RDF.
 
        Arguments (both optional):
        
        groups: Only get the RDF for atoms in these groups. Can either be an
        integer or a list integer refering to the groups (default: all groups)
        
        elements: Get the partial RDF for one element given as an atomic number
        or two elements givenas a tuple of two atomic numbers (a, b). The
        returned RDF either tells how many neighbors an atom has or how many b
        neighbors an a atom has respectively.
 
        If get_rdf is called on a newly created
        RadialDistributionFunction object, and this object has been
        created without specifying interval and average parameters, it
        is assumed that the user is not using the object as an
        observer, but wants an immediate calculation of the RDF.  In
        that case, calling get_rdf triggers the calculation of the
        RDFs.  In all other cases previously stored RDFs are returned.
        """
        if self.globalrdf is None and self.interval == 1 and self.average == 1:
            self.update()
            
        if groups is None and elements is None:
            # Return global RDF
            return self.normalize(self.globalrdf,
                                  self.countRDF * self.natoms,
                                  type=c)
 
        # Sum over selected groups
        if groups is None:
            groups = range(len(self.rdfs))
        else:
            if not isinstance(groups, list):
                groups = [groups]
 
        rdfs = self.rdfs[groups[0]].copy()  # Shallow copy
        atomcounts = self.atomcounts[groups[0]].copy()  # Shallow copy
        for i in groups[1:]: # Loop over groups
            if i >= len(self.rdfs):
                continue
            for t in self.rdfs[i].keys(): # Loop over pairs
                rdfs[t] += self.rdfs[i][t]
            for t in self.atomcounts[i].keys(): # Loop over elements
                atomcounts[t] += self.atomcounts[i][t]
 
        # Sum over selected element pairs
        elementpairs = elements
        if elementpairs is None:
            elements = self.listOfElements
            elementpairs = rdfs.keys()
        else:
            if isinstance(elementpairs, tuple):
                elements = [elementpairs[0]]
                elementpairs = [elementpairs]
            elif isinstance(elementpairs, numbers.Integral):
                elements = [elementpairs]
                elementpairs = []
                for e1 in self.listOfElements:
                    elementpairs.append((elements[0], e1))
            else:
                raise ValueError('Elements must be either an interger or ' +
                                 'a tuple of two integers.')
 
        rdf = None
        for pair in elementpairs:
            if rdf is None:
                rdf = numpy.array(rdfs[pair])
            else:
                rdf += rdfs[pair]
 
        # The atomcounts should be summed over the selected element(s)!
        atomcount = None
        for e in elements:
            if atomcount is None:
                atomcount = atomcounts[e]
            else:
                atomcount += atomcounts[e]
        if self.verbose:
            print('Number of selected atoms:', atomcount, file=sys.stderr)
 
        return self.normalize(rdf, atomcount, normalize)
 
    def normalize(self, rdf, ncount, type):
        """Normalize the raw RDF returned by the C++ module."""
        if type == 'volume':
            factor = (4 * numpy.pi / 3.0) * ncount * self.natoms / self.volume
            r = numpy.arange(self.nBins) * self.dr
            r3low = r * r * r
            r += self.dr
            r3high = r * r * r
            normalization = factor * (r3high - r3low)
            return rdf / normalization
        elif type == 'atoms':
            normalization = 1.0 * ncount
            return rdf / normalization
        elif type == 'reduced':
            factor = (4 * numpy.pi / 3.0) * ncount * self.natoms / self.volume
            r = numpy.arange(self.nBins) * self.dr
            r3low = r * r * r
            r += self.dr
            r3high = r * r * r
            normalization = factor * (r3high - r3low)            
            return 4 * numpy.pi * r * (rdf/normalization - 1) * self.natoms / self.volume
        else:
            normalization = 1.0
            return rdf / normalization
    
    # Saving and loading RDF's
    def output_file(self, prefix):
        """Give the file prefix for saving RDFs, and turn on saving."""
        self.autosave = True
        self.savefilecounter = 0
        if "%" in prefix:
            # Assume the user knows what (s)he is doing
            self.filenames = prefix
        else:
            self.filenames = prefix + "%04d.rdf"
 
    def save(self, filename=None):
        if self.verbose:
            print("Saving RDF", file=sys.stderr)
        if filename == None:
            filename = (self.filenames % (self.savefilecounter,))
            self.savefilecounter += 1
        data = {"globalrdf": self.globalrdf,
                "rdfs" : self.rdfs,
                "atomcounts": self.atomcounts,
                "countRDF": self.countRDF,
                "listOfElements": self.listOfElements,
                "rMax": self.rMax,
                "dr": self.dr,
                "nBins": self.nBins,
                "natoms": self.natoms,
                "volume": self.volume}
        f = open(filename, "wb")
        pickle.dump(data, f, 2) # Protocol 2 (efficient bin from Python 2.3)
        f.close()
