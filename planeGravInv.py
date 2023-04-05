#!/usr/bin/python3
__copyright__ = "Copyright (c) 2021 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Andrea Codd"

import importlib, sys, os
sys.path.insert(0, os.getcwd())
import argparse

from esys.escript import *
from esys.finley import ReadGmsh, ReadMesh
import esys.escript.unitsSI as U
import numpy as np
from scipy.interpolate import griddata
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions
from esys.escript.pdetools import PCG
from esys.downunder import *
from esys.weipa import *
from scipy.io import netcdf_file

def grepValuesByMaskPrint(xi, data, mask, name):
    X=data.getX()
    x=[]
    y=[]
    z=[]
    values=[]
   
    for i in range(mask.getNumberOfDataPoints()):
        if mask.getTupleForDataPoint(i)[0] > 0:
            x1 = X.getTupleForDataPoint(i)[0]
            y1 = X.getTupleForDataPoint(i)[1]
            z1 = X.getTupleForDataPoint(i)[2]
            v1 = data.getTupleForDataPoint(i)[0]
            x.append(x1)
            y.append(y1)
            z.append(z1)
            values.append(v1)
     
    if len(xi) == 2:
        r=griddata((np.array(x), np.array(y)), np.array(values), tuple(xi), method='linear',  fill_value=np.nan, rescale=False)  
        xc=np.array(xi[0].reshape(-1,1))
        yc=np.array(xi[1].reshape(-1,1))
        data = np.array(r.reshape(-1,1))     
        np.savetxt(name, np.hstack([xc, yc, data]), delimiter =',', fmt='%1.4e' )
    else:
        r=griddata((np.array(x), np.array(y), np.array(z)), np.array(values), tuple(xi), method='linear', fill_value=np.nan, rescale=False) 
        xc=np.array(xi[0].reshape(-1,1))
        yc=np.array(xi[1].reshape(-1,1))
        zc=np.array(xi[2].reshape(-1,1))        
        data = np.array(r.reshape(-1,1))     
        np.savetxt(name, np.hstack([xc, yc, zc, data]), delimiter =',', fmt='%1.4e' )    
    return

class ACGravity(object):
    def __init__(self, domain, w_e, gz_e, rho_0, mu, m0, dataGrid, densGrid,
         atol=1.0, rtol=1.0, iter_max=100, pde_tol=1e-8, output_name='solutions', 
         csv_name = 'csvs', verboseLevel="low"):

        self.domain = domain
        self.w_e = w_e
        self.gz_e = gz_e
        self.gg   = integrate(gz_e*gz_e)                                 
        self.mu   = mu*self.gg
        self.m0   = m0       
        self.atol = atol
        self.rtol = rtol
        self.iter_max = iter_max
        self.pdetol = pdetol
        self.verboseLevel = verboseLevel
        self.rho_0 = rho_0                             
        self.beta = 4.0*np.pi*U.Gravitational_Constant*self.rho_0

        self.mfs = []
        self.smooths = []
        self.output_name = output_name

        z=self.domain.getX()[2]
        top = whereZero(z-sup(z))
        self.pdeu = self.setupPDE()
        self.pdeu.setValue(q=top)
        
        # output stuff - grids and masks
        self.densGrid = densGrid
        self.dataGrid = dataGrid
        self.mskBaseRF = Scalar(0, ReducedFunction(self.domain))
        self.mskBaseRF.setTaggedValue("Base", 1)
        self.mskBaseRF.expand()
        mskDataRF = Scalar(0, ReducedFunction(self.domain))
        mskDataRF.setTaggedValue("DataArea", 1)
        mskDataRF.expand()
        self.mskDataF=interpolate(mskDataRF,Function(domain)) 
        self.outputname = output_name  
        self.csv_name = csv_name

    def setupPDE(self):
        pde=LinearSinglePDE(self.domain, isComplex=False)
        pde.setValue(A=kronecker(3))  
        pde.setSymmetryOn()
        options=pde.getSolverOptions()
        options.setTolerance(self.pdetol)
        options.setSolverMethod(SolverOptions.PCG)
        options.setPackage(SolverOptions.TRILINOS)
        options.setPreconditioner(SolverOptions.AMG)
        options.setTrilinosParameter("max levels", 10)  
        options.setTrilinosParameter("problem: symmetric", True)
        options.setTrilinosParameter("reuse: type", "full")
        return pde

    def RHS(self): 
        self.pdeu.setValue(X = self.gz_e*kronecker(3)[2], Y = Scalar(0., Solution(self.domain)))
        u = self.pdeu.getSolution()
        return ArithmeticTuple ((self.beta/self.mu)*u, Vector(0., Solution(self.domain)))         

    def Aprod(self,m):
        self.pdeu.setValue(Y = -self.beta*m, X = Vector(0., Solution(self.domain)))   
        u = self.pdeu.getSolution()                            
        self.pdeu.setValue(Y = Scalar(0., Solution(self.domain)), X = self.w_e*grad(u)[2]*kronecker(3)[2])   
        u = self.pdeu.getSolution()
        return ArithmeticTuple(-(self.beta/self.mu)*u, grad(m))     
        
    def Msolve(self,r):
        self.pdeu.setValue(X=r[1], Y=r[0])
        return self.pdeu.getSolution()       

    def bilinearform(self, m, r):
        return integrate(inner(grad(m), r[1])+ m*r[0])

    def myPCG(self, x,r,itermax,rtol):
       piter=0
       rhat=self.Msolve(r)
       d = rhat
       rhat_dot_r = self.bilinearform(rhat, r)
       if rhat_dot_r<0: print("negative norm.")
       norm_r0=np.sqrt(rhat_dot_r)
       atol2=self.rtol*norm_r0
       if atol2<=0:
          print("Non-positive tolarance.")
       print(("PCG: initial residual norm = %e (absolute tolerance = %e)"%(norm_r0, atol2)))
       mfold=0.5
       mf=0.5
       self.smooths.append(0.)
       self.mfs.append(0.5)
       donep05=False
       donep01=False
       donep008=False
       donep005=False
       while not np.sqrt(rhat_dot_r) <= atol2:
           piter+=1
           if piter  >= iter_max: 
               print("maximum number of %s steps reached."%iter_max)
               break 
           q=self.Aprod(d)
           alpha = rhat_dot_r / self.bilinearform(d, q)
           x += alpha * d
           r += q * (-alpha)      
           rhat = self.Msolve(r)
           rhat_dot_r_new = self.bilinearform(rhat, r)
           beta = rhat_dot_r_new / rhat_dot_r
           rhat += beta * d
           d = rhat

           rhat_dot_r = rhat_dot_r_new
           if rhat_dot_r<0: print("negative norm.")
           self.pdeu.setValue(Y = -self.beta*x, X = Vector(0., Solution(self.domain)))
           u = self.pdeu.getSolution()       
           u_z = grad(u)[2]
           mf=0.5*integrate((self.w_e*u_z+self.gz_e)**2)/self.gg
           if mf<0.05 and not donep05:
               print("0.05",piter)
               saveSilo(self.output_name+"p05", gravity=-u_z*self.w_e, m=x*self.rho_0)
               saveVTK(self.output_name+"p05", gravity=-u_z*self.w_e, m=x*self.rho_0)
               grepValuesByMaskPrint(self.dataGrid, -u_z, self.mskDataF,self.csv_name+"_computedgravity_p05.csv")
               grepValuesByMaskPrint(self.densGrid, x*self.rho_0, self.mskBaseRF, self.csv_name+"_density_p05.csv")               
               donep05=True
           if mf<0.01 and not donep01:
               print("0.01",piter)
               saveSilo(self.output_name+"p01", gravity=-u_z*self.w_e, m=x*self.rho_0)
               saveVTK(self.output_name+"p01", gravity=-u_z*self.w_e, m=x*self.rho_0)
               grepValuesByMaskPrint(self.dataGrid, -u_z, self.mskDataF,self.csv_name+"_computedgravity_p01.csv")
               grepValuesByMaskPrint(self.densGrid, x*self.rho_0, self.mskBaseRF, self.csv_name+"_density_p01.csv")               
               donep01=True
           if mf<0.008 and not donep008:
               print("0.008",piter)
               saveSilo(self.output_name+"p008", gravity=-u_z*self.w_e, m=x*self.rho_0)
               saveVTK(self.output_name+"p008", gravity=-u_z*self.w_e, m=x*self.rho_0)
               grepValuesByMaskPrint(self.dataGrid, -u_z, self.mskDataF,self.csv_name+"_computedgravity_p008.csv")
               grepValuesByMaskPrint(self.densGrid, x*self.rho_0, self.mskBaseRF, self.csv_name+"_density_p008.csv")               
               donep008=True
           if mf<0.005 and not donep005:
               print("0.005",piter)
               saveSilo(self.output_name+"p005", gravity=-u_z*self.w_e, m=x*self.rho_0)
               saveVTK(self.output_name+"p005", gravity=-u_z*self.w_e, m=x*self.rho_0)
               grepValuesByMaskPrint(self.dataGrid, -u_z, self.mskDataF,self.csv_name+"_computedgravity_p005.csv")
               grepValuesByMaskPrint(self.densGrid, x*self.rho_0, self.mskBaseRF, self.csv_name+"_density_p005.csv")       
               donep005=True
           self.mfs.append(mf)
           smooth=integrate(inner(grad(x),grad(x)))
           self.smooths.append(smooth)
           print('mf ', mf,' smooth ',smooth)
           print(("PCG: iteration step %s: residual norm = %e"%(piter, np.sqrt(rhat_dot_r))))
       print(("PCG: tolerance reached after %s steps."%piter))
       allsmooths=np.array(self.smooths)
       allmfs=np.array(self.mfs)
       np.savetxt(self.csv_name+"smooths.csv",allsmooths, delimiter=",")   #
       np.savetxt(self.csv_name+"mfs.csv",allmfs,delimiter=",")
       return x,

    def solve(self):
        r = self.RHS()
        if self.verboseLevel=="low":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=False)
        elif self.verboseLevel=="medium":
            m = PCG(r, self.Aprod, self.m0, self.Msolve, self.bilinearform,
                 atol=self.atol, rtol=self.rtol, iter_max=self.iter_max, 
                 initial_guess=True, verbose=True)
        elif self.verboseLevel == "high":
            m = self.myPCG(self.m0, r,self.iter_max,self.rtol)
        self.pdeu.setValue(Y = -self.beta*m[0], X = Vector(0., Solution(self.domain)))
        u = self.pdeu.getSolution()       
        u_z = grad(u)[2]
        print("    Gravity data",   self.gz_e)
        print("Computed gravity", -self.w_e*u_z)
        mf0 = 0.5*integrate((self.gz_e)**2)/self.gg 
        mf = 0.5*integrate((self.w_e*u_z+self.gz_e)**2)/self.gg
        print('Initial misfit ', mf0)
        print('  Final misfit ', mf)
        print("    difference ", mf0-mf)
        return m[0], -u_z

class DataReader(object):
    EARTH_R=6378.0088*1000.
    toobig=1.e99
    def __init__(self, filename):
        f=netcdf_file(filename, 'r', mmap=False)
        self.setReferencePoint()
        latitude=None
        longitude=None
        for latvar in "latitude", "lat", "y_range", "y":
            try:
                latitude=f.variables[latvar]
                break
            except KeyError:
                pass
        for lonvar in "longitude", "lon", "x_range", "x":
            try:
                longitude=f.variables[lonvar]
                break
            except KeyError:
                pass
        assert latitude is not None, "no latitude dimension found."
        assert longitude is not None, "no longitude dimension found."
        
        LAT_SOUTH=min(latitude)
        LAT_NORTH=max(latitude)
        print("Latitude range = %s, %s"%(LAT_NORTH, LAT_SOUTH))
        DataCenterLat=(LAT_SOUTH+LAT_NORTH)/2
        GRID_STEP_X=longitude[1]-longitude[0]
        GRID_STEP_Y=latitude[1]-latitude[0]

        # data region:
        DataSpacingX=GRID_STEP_X
        DataSpacingY=GRID_STEP_Y
        self.DataSpacing=(DataSpacingX, DataSpacingY)
        self.table=None
        for datavar in f.variables.keys():
            if datavar not in [latvar, lonvar]:
                DATA=f.variables[datavar]
                print("Gravity data variable '%s' read from %s."%(datavar,data_file))
                # convert NaN (not-a-number) values to -1000 so the plotting works
                DATA=np.where(np.isnan(DATA[:]), -1000, DATA[:])

                # make sure data is in native format
                if DATA.dtype.byteorder != '=':

                    DATA=np.array(DATA, dtype=np.float32)
                self.table=DATA
        assert self.table is not None, "no gravity data found."
        print("data range : %s - %s"%(self.table.min(), self.table.max()))

    def setReferencePoint(self, Xref=(0., 0.)):
        self.Xref=Xref
        
    def extractData(self, domain):
        X=Function(domain).getX()
        result=interpolateTable(self.table, X[:2], self.Xref, self.DataSpacing, self.toobig)
        print('extract data')
        print(result)
        return result

########################################################################
### Input files and variables from file 
parser = argparse.ArgumentParser(description='Gravity inversion for plane data in netcdf format.', epilog="version 01/2021 by a.codd@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')
args = parser.parse_args()
config = importlib.import_module(args.config)
print("Configuration "+args.config+".py imported.")

mu       = config.mu
rho_0    = config.rho_0
atol     = config.atol  
rtol     = config.rtol
pdetol   = config.pdetol
iter_max = config.iter_max 
data_file= config.data_file
data_scale = config.data_scale

# build domain 
filename, file_extension = os.path.splitext(config.mesh_name)
if file_extension == ".msh":
    dom=ReadGmsh(config.mesh_name, numDim = 3)
else:
    dom=ReadMesh(config.mesh_name, numDim = 3)
print("Mesh read from "+config.mesh_name)

# define data area
w_e = Scalar(0,Function(dom))        
w_e.setTaggedValue("DataArea",1)
w_e.expand()

# define ground
rho_e = Scalar(0,Function(dom))    
rho_e.setTaggedValue("Base",1)
rho_e.setTaggedValue("PaddingBase",1)
rho_e.expand()
rho_e = rho_0*rho_e

# read in data and make g
dr=DataReader(config.data_file)
gz_e = dr.extractData(dom)*wherePositive(w_e)*data_scale#/1.e6   #### convert from micrometres/s^2 to m/s^2
del dr

# grids for output computed gravity and density
# from mesh
dataRefX = config.DataRefX
dataRefY = config.DataRefY
DataSpacingX = config.DataSpacingX
DataSpacingY = config.DataSpacingY
spZ = config.DataMeshSizeVertical
coreDepth = config.CoreThickness
VerticalMeshCB = config.MeshSizeCoreFactor*spZ
numX = config.DataNumX
numY = config.DataNumY

spanX = dataRefX + (numX-1)*DataSpacingX
spanY = dataRefY + (numY-1)*DataSpacingY

print(spanX,spanY)


print("coreDepth",coreDepth)
print("spZ/2",spZ/2)
print("VerticalMeshCB",VerticalMeshCB)

totZD = coreDepth - spZ/2. - VerticalMeshCB
print("totZD",totZD)
numZ = 1 + int(np.floor(totZD/spZ))
print("numZ",numZ)
totZD = -(numZ-1)*spZ-spZ/2.
print("totZD",totZD)

# data grid
DXg = np.linspace(dataRefX + 0.5*DataSpacingX, spanX-0.5*DataSpacingX, numX-1, endpoint=True)
DYg = np.linspace(dataRefY + 0.5*DataSpacingY, spanY-0.5*DataSpacingY, numY-1, endpoint=True) 
dataGrid = np.meshgrid(DXg, DYg)

# density grid
DX = np.linspace(dataRefX + 1.5*DataSpacingX, spanX-1.5*DataSpacingX, numX-3, endpoint=True)
DY = np.linspace(dataRefY + 1.5*DataSpacingY, spanY-1.5*DataSpacingY, numY-3, endpoint=True)
DZbase = np.linspace(totZD, -spZ/2, numZ, endpoint=True)
densGrid = np.meshgrid(DX, DY, DZbase)

# initial guess
m0 = Scalar(0., ContinuousFunction(dom)) 

grav = ACGravity(dom, w_e, gz_e, rho_e, mu, m0, dataGrid, densGrid,
                   atol, rtol, iter_max, pdetol, config.output_name, config.csv_name, config.VerboseLevel)
m, gz =grav.solve()

wgz = w_e*gz
RMS=np.sqrt(integrate(wherePositive(w_e)*(gz-gz_e)**2)/integrate(wherePositive(w_e)))
print("RMS",RMS)
gravdiff=wgz-gz_e
densitydiff=m*rho_e
saveSilo(config.output_name+"_final", gravity=gz, wgz=wgz, gravidiff=gravdiff, m=m, data=gz_e, densitydiff=densitydiff)
saveVTK(config.output_name+"final", gravity=gz, wgz=wgz, gravidiff=gravdiff, m=m, data=gz_e, densitydiff=densitydiff)
mskBaseRF = Scalar(0, ReducedFunction(dom))
mskBaseRF.setTaggedValue("Base", 1)
mskBaseRF.expand()
mskDataRF = Scalar(0, ReducedFunction(dom))
mskDataRF.setTaggedValue("DataArea", 1)
mskDataRF.expand()
mskDataF=interpolate(mskDataRF,Function(dom)) 
grepValuesByMaskPrint(dataGrid, gz, mskDataF, config.csv_name+"_computedgravity_final.csv")
grepValuesByMaskPrint(densGrid, densitydiff, mskBaseRF, config.csv_name+"_density_final.csv")   

print('results silo saved to '+config.output_name+'_final.silo')
print("finished")




    
