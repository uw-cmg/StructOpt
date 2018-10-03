import numpy
import os
import sys
import copy
import random
from tempfile import mkdtemp
import math
import cmath
#matplotlib.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm  as cm
from PIL import Image, ImageChops
from ase.io import write
from ase import Atom, Atoms
import scipy.interpolate
import scipy.special
import scipy.misc
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from scipy import optimize
import time
import json
import logging
try:
    from mpi4py import MPI
except:
    pass

class STEM_eval(object):
    def __init__(self):
        self.args = self.read_inputs()
        for k,v in self.args.items():
            setattr(self,k,v)

    def read_inputs(self):
        args = json.load(open('stem_inp.json'))
        with open(args['psf']) as fp:
            psf_data = fp.read() 
        args['psf'] = numpy.reshape(numpy.array(psf_data.split(),dtype=float),[args['pixels']]*2)
        return args

    def evaluate_fitness(self, Optimizer, individ):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()

        if rank==0:
            ntimes=int(math.ceil(1.*len(individ)/comm.Get_size()))
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
                out = (None, stro)
            else:
                out = self.evaluate_indiv(Optimizer, ind, rank)

            outt = comm.gather(out,root=0)
            if rank == 0:
                outs.extend(outt)
        return outs

    def evaluate_indiv(self, Optimizer, individ, rank):

        logger = logging.getLogger(Optimizer.loggername)

        logger.info('Received individual HI = {0} for STEM evaluation'.format(
            individ.history_index))
        atms = individ[0]
        expfun = self.calculate_simp_function(self.stemref2image(str(self.stem_ref)))
        #Get the image from ConvStem
        try:
           Grid_sim2exp = self.grid_sim2exp
        except:
           Grid_sim2exp = 1
        if self.pixelshift == True: 
            points = [-0.5993457/2,0,0.5993457/2]
            chisq = []
            for x in points:
                for y in points:
                    pixelshift = [x,y,0]
                    simfun=self.get_image(self.psf,atms,self.slice_size,self.pixels,pixelshift)
                    simfun_resize = simfun
                    chisq.append(self.compare_functions(expfun,simfun_resize))
            return min(chisq)       

        else:
            simfun=self.get_image(self.psf,atms,self.slice_size,self.pixels)
            if Grid_sim2exp > 1: 
                 #simfun_resize = numpy.zeros([150,150],dtype=float)       
                simfun_resize = numpy.zeros([self.pixels,self.pixels],dtype=float)
                for x in range(len(simfun)):
                    for y in range(len(simfun[0])):
                        simfun_resize[x/Grid_sim2exp][y/Grid_sim2exp]+=simfun[x][y]
                simfun_resize /= 0.5993457**2
                simfun_resize *= 0.12**2
            else:
                simfun_resize = simfun

            chisq=self.compare_functions(expfun,simfun_resize)  
        logger.info('M:finish chi2 evaluation, chi2 = {0} @ rank ={1}'.format(chisq,rank))
        signal = 'Evaluated individual {0} on {1}\n'.format(individ.index,rank)
        return chisq, signal

    def calculate_simp_function(self,image):
        return image
    
    def calculate_comp_function(self,imagefile):
        img=Image.open(imagefile).convert('LA')
        data=numpy.zeros((len(a),len(a[0])))
        for i in range(len(a)):
            for j in range(len(a[0])):
                data[i][j]=a[i][j]
        
        #Calculate the atom column positions
        positions=get_atom_pos(self,data)
        x=positions[0]
        y=positions[1]
        
        #Apply rotation/translation
        mvs=[]
        for i in range(len(x)):
            mvs.append(data[y[i],x[i]])
        top2=sorted(enumerate(mvs), key=lambda a:a[1], reverse=True)[:2]
        #Set new origin at max value
        #Translation vector becomes:
        trans=[x[top2[0][0]],y[top2[0][0]]]
        #Find rotation angle alpha
        slope=float(y[top2[0][0]]-y[top2[1][0]])/float(x[top2[0][0]]-x[top2[1][0]])
        alpha=math.atan(abs(slope))
        #Adjust for quadrant
        if x[top2[1][0]] < x[top2[0][0]]:
            if y[top2[1][0]] > y[top2[0][0]]:
                #quadrant 2
                alpha=alpha+math.pi/2.0
            else:
                #quadrant 3
                alpha=alpha+math.pi
        else:
            if y[top2[1][0]] < y[top2[0][0]]:
                #quadrant 4
                alpha=alpha+math.pi*3.0/2.0
        xs=numpy.arange(len(data))
        xys=[]
        z=[]
        for i in range(len(xs)):
            for j in range(len(xs)):
                ara=numpy.array([xs[i], xs[j], 1])
                xys.append(ara)
                z.append(data[j][i])
        transfmat=numpy.array([[math.cos(alpha), -1.0*math.sin(alpha), -1.0*x[top2[0][0]]],[math.sin(alpha),math.cos(alpha),-1.0*y[top2[0][0]]],[0.0,0.0,1.0]])
        nx=[]
        ny=[]
        for one in xys:
            ara=numpy.dot(transfmat,one)
            nx.append(ara[0])
            ny.append(ara[1])
        
        gfuntot=scipy.interpolate.RectBivariateSpline(x,y,z)
        
        return gfuntot
    
    def compare_functions(self, expfun, simfun):
        """Function compares two matrices and calculates chisq.  Matrices must be same size"""
        try:
          Grid_sim2exp = self.grid_sim2exp
        except:
          Grid_sim2exp = 1
        #if self.args['Pixelshift'] == True:
        chisq = 0.0
        for i in range(len(expfun)):
              for j in range(len(expfun[0])):
                 if expfun[i][j] > -1000:
                    chisq +=(simfun[i][j]-expfun[i][j])**2
        #else:
        #   chisq = abs(sum(sum((simfun-expfun)**2)))
        chisq = chisq/len(simfun)/len(simfun[0])
        return chisq
    
    def stemref2image(self,stemref):
        if stemref.split('/')[-1].endswith('xyz'): # xyz coordinates mean a phantom
            from StructOpt import inp_out
            atoms_ref=inp_out.read_xyz(stemref,0)
            return self.get_image(self.psf,atoms_ref,self.slice_size,self.pixels)
        else: # experimental image
            fileobj = open(stemref.split('/')[-1], 'r')
            lines = fileobj.readlines()
            nl = len(lines)
            Ims = numpy.empty([nl,nl],dtype=float)
            for x in range(0,len(lines)):
                Ims[x] = lines[x].split() 
            
            nk = self.pixels
            Icfit = numpy.zeros([nk,nk],dtype=float)
            if nl < nk:
                raise RuntimeError("Required pixels should be smaller than provided experimental image pixel!")
            for x in range((nl-nk)/2,(nl-nk)/2+nk): # cropping the image if nl != nk, nl > nk
                for y in range((nl-nk)/2,(nl-nk)/2+nk):
                    Icfit[x-(nl-nk)/2][y-(nl-nk)/2] = self.IcImsfit(Ims[x][y],self.fitting_coeff)
                    if Icfit[x-(nl-nk)/2][y-(nl-nk)/2] < 0:
                        Icfit[x-(nl-nk)/2][y-(nl-nk)/2] = -1000
            return Icfit

    def IcImsfit(self,Ims,coeff): 
        # Fitting of convolutional method vs multislice method (this is a polynomial fit)`
        Ic = 0.0;
        for i in range(len(coeff)):
            Ic = Ic + coeff[i]*Ims**i
        return Ic

    
    def get_image(self, psf, atms, rmax, nx, pixelshift=[0,0,0], scalefactor=1.0):
        """Function to get image based point spread function and atoms
        rmax=Size of slice in Angstoms
        nx=number of pixels in slice"""

        try:
            rank = MPI.COMM_WORLD.Get_rank()
        except:
            rank = 0
        
        psf2d = numpy.fft.fft2(psf)

        pot = self.stempot(rmax,rmax,len(psf),len(psf[0]),atms,pixelshift,scalefactor)

        potft = numpy.fft.fft2(pot)

        potm = potft*psf2d

        zcon_im = numpy.fft.ifft2(potm,axes=(0,1)).real

        return zcon_im
    
    def get_atom_pos(self, data):
		"""Function to identify the location of the atom columns
		Inputs are:
			data = 2D matrix with imagefile value
			self.args:
				neighborhood_size = size of square stamp to search array
				threshold = difference for change"""
		
	
		if 'neighborhood_size' in self.args:
			neighborhood_size = self.args['neighborhood_size']
		else:
			neighborhood_size = 30
		if 'threshold' in self.args:
			threshold = self.args['threshold']
		else:
			threshold = 30
		
		#Use filters to calculate peaks
		data_max = filters.maximum_filter(data, neighborhood_size)
		maxima = (data == data_max)
		data_min = filters.minimum_filter(data, neighborhood_size)
		diff = ((data_max - data_min) > threshold)
		maxima[diff == 0] = 0

		labeled, num_objects = ndimage.label(maxima)
		slices = ndimage.find_objects(labeled)
		x, y = [], []
		for dy,dx in slices:
			x_center = (dx.start + dx.stop - 1)/2
			x.append(x_center)
			y_center = (dy.start + dy.stop - 1)/2    
			y.append(y_center)
		

		posiitons=[x,y]
		
		return positions


        
    def get_probe_function(self, args):
        """Function to get the probe function based on input args
        kev=Electron energy in keV
        ap = Objective aperture semiangle in mrad
        Cc = Chromatic aberration coefficient in mm
        dE = Delta E in eV
        ds=Source size in Angstroms
        rmax=Size of slice in Angstroms
        Cs=Spherical aberration Cs in mm
        df=Defocus in Angstroms
        nx=number of pixels in slice"""
        
        keV = args['electron_energy']
        ap = args['aperture_semiangle']
        try:
            Cc = args['chromatic_aberration_coefficient']
        except:
            Cc = 0.0
        try:
            dE = args['delta_E']
        except:
            dE = 0.0
        #Cs=args['Spherical aberration']
        #df=args['Defocus']
        ds=args['source_size']
        rmax=args['slice_size']
        nx=args['pixels']
       
        try:
            aber=args['aber']
        except:
            aber=[[0,0] for i in range(12)]
         
        nk = math.floor(40*(0.001*ap / self.wavlen(keV))*rmax)
        if math.fmod(nk, 2):
            nk += 1
        if nx > nk:
          nk = nx
        #Get probe function
        if Cc == 0.0 and dE == 0.0:
            psf2D = self.STEMPSF2DCoh(aber, keV, ap, nk, rmax)
            #Redimension/D/C probe2DCoh #Converts wave to double precision complex wave
            #wave/c psf2D = $"probe2DCoh"
        else:
            psf2D = self.STEMPSF2DIncoh(aber, keV, Cc, dE, ds, ap, nk, rmax)
            #Redimension/D/C probe2DIncoh #Converts wave to double precision complex wave
            #wave/c psf2D = $"probe2DIncoh"
        #Shift minimum to zero
        #mn = numpy.minimum.reduce(numpy.minimum.reduce(psf2D))
        #for i in range(len(psf2D)):
        #    for j in range(len(psf2D[0])):
        #        psf2D[i][j] += abs(mn)
        
        return psf2D

    def STEMPSF2DCoh(self, aber, keV, ap, nk,rmax):
        kmax = 0.001*ap / self.wavlen(keV)
        phase = self.ChiPhase2D(aber, keV, ap, nk/10)
        xp=numpy.linspace(-kmax*2.0,(kmax*2.0),nk/10)
        yp=numpy.linspace(-kmax*2.0,(kmax*2.0),nk/10)
        phasef = scipy.interpolate.RectBivariateSpline(xp,yp,phase)
        #print "M:kmax",kmax
        #print phasef
        
        kbound = nk/4.0/(rmax/2.0)
        # need to pad the phase with zeros here.
        probe2DCoh = numpy.zeros((nk,nk),dtype=complex)
        x=numpy.linspace(-kbound,kbound,nk)
        y=numpy.linspace(-kbound,kbound,nk)
        
        #SetScale/I x -20*kmax, 40*kmax, "", probe2DCoh
        #SetScale/I y -20*kmax, 40*kmax, "", probe2DCoh
        for i in range(len(x)):
            for j in range(len(y)):
                if x[i]**2 + y[j]**2 < kmax**2:
                    probe2DCoh[i][j] = complex(1, phasef(x[i],y[j]))
                else:
                    probe2DCoh[i][j] = complex(0, 0)
        #probe2DCoh = ( (x^2 + y^2 < kmax^2) ?  p2rect(cmplx(1, phase(x)(y)))  : cmplx(0, 0))

        p2dcf = numpy.fft.fft2(probe2DCoh)
        for i in range(len(x)):
            for j in range(len(y)):
                probe2DCoh[i][j] = complex((abs(p2dcf[i][j]))**2, 0)
        probe2DCoh = numpy.real(probe2DCoh)    
        #print "M:probe2DCoh" 
        #print probe2DCoh[0]
        #print "M:sum:probe2DCoh",sum(sum(probe2DCoh)) 
        probe2DCoh /= sum(sum(probe2DCoh))
        return probe2DCoh

    def STEMPSF2DIncoh(self, aber, keV, Cc, dE, ds, ap, nk, rmax):
        kmax = 0.001*ap / self.wavlen(keV)
        defocus_distribution = self.ChromaticDefocusDistribution(Cc, dE, keV, ap)
        start_df = aber[0][0]
        df_range = 2.5*(Cc*1e7*dE /(1e3*keV)*( (1+keV/511.0)/(1+keV/1022.0) ))
        x = numpy.linspace(-df_range,df_range, num=len(defocus_distribution))
        for i in range(len(defocus_distribution)):
            #aber_step = numpy.copy(aber)
            aber[0][0] = start_df + x[i]
            #print aber_step
            probe2DCoh = self.STEMPSF2DCoh(aber, keV, ap, nk, rmax)
            if i==0:
               # probe2DIncoh = copy.deepcopy(probe2DCoh)
                probe2DIncoh = numpy.copy(probe2DCoh)
                probe2DIncoh *= defocus_distribution[0]
            else:
                for j in range(len(probe2DIncoh)):
                    for k in range(len(probe2DIncoh[0])):
                        probe2DIncoh[j][k] += defocus_distribution[i]*probe2DCoh[j][k]
            #print i
            #print probe2DIncoh

        probe2DIncoh /= sum(defocus_distribution)
        #print "probe2DIncoh",probe2DIncoh[0],probe2DIncoh[100]
        #print "ave",sum(defocus_distribution),probe2DIncoh
        aber[0][0] = start_df
        
        if ds != 0:
            fds = ds / (2*(2*math.log(2))**0.5)    # real-space standard deviation for Gaussian with FWHM ds
            #fds = 1/(2*math.pi*fds)    # FT standard deviation
            #Redimension/C probe2DIncoh
            #wave/C probe2DSS = $"probe2DIncoh"
            probe2DSS = numpy.fft.fft2(probe2DIncoh)
            #g = Gauss(0.0, 0.0, fds, fds)#, len(probe2DSS), len(probe2DSS[0]))
            #xs = numpy.linspace(-2.19234, 2.19234, num=len(probe2DSS))  #256 pts
            #xs = numpy.linspace(-4.76855, 4.76855, num=len(probe2DSS))   #560 pts
            xs = numpy.linspace(-rmax/2, rmax/2, num=len(probe2DSS))   #560 pts
            #xs = numpy.linspace(-31.9809, 31.9809, endpoint=False, num=len(probe2DSS))   #3750 pts
            #leftend = -60*kmax/(len(probe2DSS)-1)*((len(probe2DSS)-1)/2+1)  # nk =258 left 257->129 points
            #rightend = 60*kmax/(len(probe2DSS)-1)*((len(probe2DSS)-1)/2)    # nk =258 right 257->128 points
            #xs = numpy.linspace(leftend, rightend, num=len(probe2DSS))
            #print 'len',len(probe2DSS),xs[0],xs[1],xs[len(probe2DSS)-1]
            #print 'Gauss',xs[124],Gauss2D(xs[124],0.0,fds,xs[124],0.0,fds)
            Gaussfunction = numpy.zeros((nk,nk),dtype=complex)
            for j in range(len(probe2DSS)):
                    for k in range(len(probe2DSS[0])):
                        Gaussfunction[j][k]=complex(self.Gauss2D(xs[j],0.0,fds,xs[k],0.0,fds),0.0)
                       # if Gaussfunction[j][k] > 0.1:
                      #     print j,k, Gaussfunction[j][k]
            #probe2DSS *= complex(Gauss(x, 0.0, fds, y, 0.0, fds), 0.0)
            Gaussft = numpy.fft.fft2(Gaussfunction)
            probe2DSS *= Gaussft
            probe2Dreal = numpy.fft.ifft2(probe2DSS)
            probe2DSS = numpy.real(probe2Dreal)
            #print "probe2DSS",probe2DSS[0],probe2DSS[100]
            probe2DSS = numpy.fft.fftshift(probe2DSS)
            #probe2DSS = cmplx(sqrt(magsqr(probe2DSS)), 0)
            probe2DSS /= sum(sum(probe2DSS))
            probe2DIncoh = numpy.copy(probe2DSS)
            #SetScale/p x dimoffset(probe2dcoh, 0), dimdelta(probe2dcoh, 0), "", probe2dIncoh
            #SetScale/P y dimoffset(probe2dcoh, 1), dimdelta(probe2dcoh, 1), "", probe2dIncoh
             
        return probe2DIncoh

    def ChiPhase2D(self, aber_in, keV, ap, nk):
        aber = self.SwitchToAngstroms(aber_in)
        wl = self.wavlen(keV)
        kmax = (0.001*ap / wl)  # maximum k through aperture
        x=numpy.linspace(-kmax*2.0,(kmax*2.0),nk)
        y=numpy.linspace(-kmax*2.0,(kmax*2.0),nk)
        astack = [numpy.zeros((len(x),len(y))) for one in range(12)]
        # Evaluate the phase shifts from the various aberrations
        for i in range(len(x)):
            for j in range(len(y)):
                #C1, defocus
                astack[0][i][j] = (1.0/2.0)*wl*aber[0][0]*(x[i]**2.0 + y[j]**2.0)
                #A1, 2-fold astigmatism
                astack[1][i][j] = (1.0/2.0)*wl*aber[1][0]*(x[i]**2.0 - y[j]**2.0)
                #A2, 3-fold astigmatism
                astack[2][i][j] = (1.0/3.0)*wl**2.0*aber[2][0]*(x[i]**3.0 - 3.0*x[i]*y[j]**2.0)
                #B2, axial coma
                astack[3][i][j] = wl**2.0*aber[3][0]*(x[i]**3.0 + x[i]*y[j]**2.0)
                #C3, primary spherical aberration
                astack[4][i][j] = (1.0/4.0)*wl**3.0*aber[4][0]*(x[i]**4.0 + 2.0*x[i]**2.0*y[j]**2.0 + y[j]**4.0)
                #A3, 4-fold astigmatism
                astack[5][i][j] = (1.0/4.0)*wl**3.0*aber[5][0]*(x[i]**4.0 - 6.0*x[i]**2.0*y[j]**2.0 + y[j]**4.0)
                #S3, star aberration
                astack[6][i][j] = wl**3.0*aber[6][0]*(x[i]**4.0 - y[j]**4.0)
                #A4, 5-fold astigmatism
                astack[7][i][j] = (1.0/5.0)*wl**4.0*aber[7][0]*(x[i]**5.0 - 10.0*x[i]**3.0*y[j]**2.0 + 5.0*x[i]*y[j]**4.0)
                #D4, 3-lobe aberration
                astack[8][i][j] = wl**4.0*aber[8][0]*(x[i]**5.0 - 2.0*x[i]**3.0*y[j]**2.0 - 3.0*x[i]*y[j]**4.0)
                #B4, axial coma
                astack[9][i][j] = wl**4.0*aber[9][0]*(x[i]**5.0 + 2.0*x[i]**3.0*y[j]**2.0 + x[i]*y[j]**4.0)
                #C5, 5th order spherical aberration
                astack[10][i][j] = (1.0/6.0)*wl**5.0*aber[10][0]*(x[i]**6.0 + 3.0*x[i]**4.0*y[j]**2.0 + 3.0*x[i]**2.0*y[j]**4.0 + y[j]**6.0)
                #A5, 5th order spherical aberration
                astack[11][i][j] = (1.0/6.0)*wl**5.0*aber[11][0]*(x[i]**6.0 - 15.0*x[i]**4.0*y[j]**2.0 + 15.0*x[i]**2.0*y[j]**4.0 - y[j]**6.0)
        
        #Set minimum to zero
        #nnastack = [numpy.zeros((len(x),len(y))) for one in range(12)] 
        # for i in range(len(x)):
    #         for j in range(len(y)):
    #             for k in range(len(astack)):
    #                 if astack[k][i][j] < 0.0:
    #                     astack[k][i][j]=0.0
        #Shift minimum to zero
    #    for k in range(len(astack)):
    #        mn = numpy.minimum.reduce(numpy.minimum.reduce(astack[k]))
    #        for i in range(len(x)):
    #            for j in range(len(y)):
    #                astack[k][i][j]+= abs(mn)
    #    print 'mn',mn
        # rotate the phase shifts of the non-centrosymmetric aberrations
        for i in range(12):
            if aber[i][1] != 0.0:
                astack[i] = scipy.ndimage.rotate(astack[i],aber[i][1],reshape=False)
                #ImageTransform/P=(i) getplane astack
                #wave aphase = $"M_ImagePlane"
                #ImageRotate/E=0/O/A=(aber[i][1]) aphase
                #SetScale/P x -DimSize(aphase, 0)*DimDelta(astack, 0) / 2, DimDelta(astack, 0), "", aphase
                #SetScale/P y -DimSize(aphase, 1)*DimDelta(astack, 1) / 2, DimDelta(astack, 0), "", aphase
                #astack[][][i] = aphase(x)(y)
           # print "M:astack",i,aber[i][0],aber[i][1]
           # print astack[i][0]

        # sum all the aberration contributions
        fsum = numpy.zeros((len(x),len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                fsum[i][j] = sum([one[i][j] for one in astack])*2*math.pi
        #MatrixOp/O phase = 2*Pi*sumbeams(astack)
        #SetScale/I x -2*kmax, 2*kmax, "", phase
        #SetScale/I y -2*kmax, 2*kmax, "", phase
        
        return fsum

    def ChromaticDefocusDistribution(self, Cc, dE, keV, ap):
        Cc *= 1e7  # mm to Angstroms
        df_phase_max = 2*math.pi / 50.0  # maximum phase step at the aperture edge due to chromatic aberration
        kmax = 0.001*ap/self.wavlen(keV)
        # defocus range and form from Reimer
        H = (Cc*dE /(1e3*keV)*( (1+keV/511.0)/(1+keV/1022.0) ))
        N = (2*(math.log(2)**0.5) / ((math.pi**0.5)*H))
        df_range = 2.5*H
        ndf = math.ceil(df_range * self.wavlen(keV) * kmax**2 / df_phase_max)
        if ndf < 31:
            ndf = 31
        
        #ndf = (ndf < 31 ? 31 : ndf)
        if not math.fmod(ndf,2):
            ndf+=1
        
        #ndf = (!mod(ndf, 2) ? ndf+1 : ndf)
        defocus_distribution = numpy.linspace(-df_range, df_range, num = ndf)
        #Make/O/N=(ndf) defocus_distribution
        #SetScale/I x -df_range, df_range, "", defocus_distribution
        for x in range(len(defocus_distribution)):
            defocus_distribution[x] = N*math.exp( -math.log(2) * (2*defocus_distribution[x] / H)**2 )
        return defocus_distribution

    def stempot(self,xmax,ymax,nx,ny,atms,pixelshift,scalefactor):
        """function V = stempot(xmax,ymax,nx,ny,potfile)
        %  STEMPOT Generate a Projected Potential
        %  inputs xmax, ymax are the size of the slice in angstroms
        %  nx,ny are the number of pixels in the x and y directions """
        #zed=2 for rutherford scattering of the nucleus, less for screening
        zed = 1.7

        ix = numpy.arange(1.0,nx)
        iy = numpy.arange(1.0,ny)
        dx = xmax/nx
        dy = ymax/ny
        rx = numpy.arange(0,xmax-dx,dx)
        ry = numpy.arange(0,ymax-dy,dy)

        Zatom = atms.get_atomic_numbers()
        #translate atoms such that the center of mass is in the center of the computational cell
        com = atms.get_center_of_mass()
       # com = [ 44.40963074 , 44.65497562 , 44.90406073] #for AuNP
       # com = numpy.array(com)
        #print 'com',com -0.149836425, 0.29967285, 0
        #com += [0.41205016875, 0.6742639125, 0] #for rotated line profile 
      #  com += [-0.149836425, 0.29967285, 0]  #for AuNP
        #com += pixelshift
        #print 'com+pixelshift',com
        cop = xmax/2.0
        trans = [cop-i for i in com]
        atms.translate(trans)
        positions=atms.get_positions()
        ax=[]
        ay=[]
        az=[]
        for o,t,h in positions:
            ax.append(o)
            ay.append(t)
            az.append(h)
        ax = numpy.array(ax)
        ay = numpy.array(ay)
        az = numpy.array(az)
        amax = len(Zatom)

        #find boundaries of slice
        axmin = min(ax)
        axmax = max(ax)
        aymin = min(ay)
        aymax = max(ay)
        
        V= numpy.zeros((nx,ny))

        #map x and y coords of the atoms to the nearest grid points
        #A fraction of the atom must be assigned to the closest gridpoints
        #to avoid sum and difference frequencies appearing in the image
        #grid point to the left of the atom
        ix = numpy.array([math.floor(axi/dx) for axi in ax])
        #apply periodic boundary conditions
        iax = numpy.array([math.fmod(iaxi,nx) for iaxi in ix])
        ibx = numpy.array([math.fmod(iaxi+1,nx) for iaxi in ix])
        #fraction of atom at iax
        fax = numpy.array([1-math.fmod((axi/dx),1 ) for axi in ax])
        #grid point above the atom
        iy = numpy.array([math.floor(ayi/dy) for ayi in ay])
        #apply periodic boundary conditions
        iay = numpy.array([math.fmod(iayi,ny) for iayi in iy])
        iby = numpy.array([math.fmod(iayi+1,ny) for iayi in iy])
        #fraction of atom at iay 
        fay = numpy.array([1-math.fmod((ayi/dy),1 ) for ayi in ay])
        #Add each atom to the potential grid
        V1 = numpy.array([fax[i] * fay[i] * (Zatom[i]**zed) for i in range(len(fax))])
        V2 = numpy.array([(1-fax[i]) * fay[i] * (Zatom[i]**zed) for i in range(len(fax))])
        V3 = numpy.array([fax[i] * (1-fay[i]) * (Zatom[i]**zed) for i in range(len(fax))])
        V4 = numpy.array([(1-fax[i]) * (1-fay[i]) * (Zatom[i]**zed) for i in range(len(fax))])
    #     V1 = numpy.array([fax[i] * fay[i] * scalefactor for i in range(len(fax))])
    #     V2 = numpy.array([(1-fax[i]) * fay[i] * scalefactor for i in range(len(fax))])
    #     V3 = numpy.array([fax[i] * (1-fay[i]) * scalefactor for i in range(len(fax))])
    #     V4 = numpy.array([(1-fax[i]) * (1-fay[i]) * scalefactor for i in range(len(fax))])

        for j in range(amax):
           V[iax[j],iay[j]] += V1[j]
           V[ibx[j],iay[j]] += V2[j]
           V[iax[j],iby[j]] += V3[j]
           V[ibx[j],iby[j]] += V4[j]
        rev_trans = [-1.0*i for i in trans]
        atms.translate(rev_trans)
        return V

    def wavlen(self, keV):
        #calculate electron wavelength in Angstroms given the energy in keV.
        wav = 12.3986 / ( (2*511.0 + keV) * keV)**0.5
        return wav

    def SwitchToAngstroms(self, aber):
        '''Function to switch aberrations wave into Angstroms instead of natural units'''
        naber = numpy.copy(aber)
        #C1, A1, A3, B2, all start in nm
        for i in range(0,4):
            naber[i][0] = 10*aber[i][0]
        #C3, A3, S3, A4, D4, B4 start in um
        for i in range(4,10):
            naber[i][0] = 1e4*aber[i][0]
        #C5, A5 in mm
        for i in range(10,12):
            naber[i][0] = 1e7*aber[i][0]
        return naber

    def Gauss2D(self, x, center_x, width_x, y, center_y, width_y, height=1.0):
        """Returns a 2D gaussian """
        g = math.exp(-0.5*((center_x-x)/width_x)**2)/(width_x*(2.0*math.pi)**0.5)
        g *= math.exp(-0.5*((center_y-y)/width_y)**2)/(width_y*(2.0*math.pi)**0.5)
        g *= height
        return g

    def Gauss(self, center_x, center_y, width_x, width_y, height=1.0):
        """Returns a gaussian function with the given args"""
        width_x = float(width_x)
        width_y = float(width_y)
        return lambda x,y: height*math.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

    def moments(self, data):
    	"""Returns (height, x, y, width_x, width_y)
    		the gaussian args of a 2D distribution by calculating its moments """
    	total = data.sum()
    	X, Y = indices(data.shape)
    	x = (X*data).sum()/total
    	y = (Y*data).sum()/total
    	col = data[:, int(y)]
    	width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    	row = data[int(x), :]
    	width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    	height = data.max()
    	return x, y, width_x, width_y, height

    def fitgaussian(self, data):
    	"""Returns (height, x, y, width_x, width_y)
    		the gaussian args of a 2D distribution found by a fit"""
    	params = self.moments(data)
    	errorfunction = lambda p: ravel(self.Gauss(*p)(*indices(data.shape)) - data)
    	p, success = optimize.leastsq(errorfunction, params)
    	return p

    def find_stem_coeff(self, Optimizer, indiv):
        from StructOpt.tools.eval_energy import eval_energy
        outs = eval_energy([Optimizer,indiv])
        indiv.energy = outs[0]
        stro=outs[3]
        if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
            indiv.bulki = outs[1]
        chisq = Optimizer.stemcalc.run(indiv[0])
        aflag=True
        alpha = 1.0
        while True:
            value=alpha*chisq
            div=abs(indiv.energy)/value
            if div <1:
                alpha*=0.1
            elif div >10:
                alpha*=10
            else:
                break
        indiv.fitness=indiv.energy/indiv[0].get_number_of_atoms()+alpha*chisq
        return alpha, indiv
