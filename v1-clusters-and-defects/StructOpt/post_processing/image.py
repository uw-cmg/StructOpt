import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import ticker
from StructOpt.tools.StemCalc import ConvStem
from StructOpt import Optimizer
from StructOpt import inp_out
import numpy
import ast
fp = open('structoptinput.txt','r')
data = fp.readlines()
for i in range(len(data)):
    data[i] = data[i].strip()
data=' '.join(data)
data='{'+data+'}'
parameters = ast.literal_eval(data)
Au=inp_out.read_xyz('STEM_ref.xyz',0)

fileobj = open(parameters['psf'],'r')
lines = fileobj.readlines()
nk = len(lines)
parameters['pixels'] = nk
parameters['psf'] = numpy.empty([nk,nk],dtype=float)
for x in range(0,nk):   
    parameters['psf'][x] = lines[x].split()
A = ConvStem(parameters=parameters,calc_exp=False)
imAu = A.get_image(parameters['psf'], Au, parameters['slice_size'], parameters['pixels'])
IPlot = imAu.T
plt.imshow(IPlot,origin="lower",cmap = plt.get_cmap('hot'))
cb=plt.colorbar()
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()
for t in cb.ax.get_yticklabels():
     t.set_fontsize(24)
plt.axis('off')
plt.savefig('image.eps')

