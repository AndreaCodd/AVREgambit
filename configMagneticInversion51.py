#
# gambit
#
# Configuration file for gravity inversion for use by planeMagneticInv.py 
# mesh has been made with mkGeoWithData2D.py
#
# Inversion constants:
#
# scale between misfit and regularization
mu       = 1.e-16
#
# used ti scale computed suscetibility
k_0      = 1.
#
# IPCG tolerance *|r| <= atol+rtol*|r0|*  (energy norm)
# absolute tolerance for IPCG interations
atol     = 0.        
#
# relative tolerance for IPCG iterations   
rtol     = 1.e-3
#
# tolerance for solving PDEs 
# make sure this is not more than the square of rtol     
pdetol   = 1.e-10    
#
# maximum number of IPCG iterations
iter_max = 5000 
#
# data scale.
data_scale = 1.
#
# Background magnetic field
# magnitude
Magmag =  58014.0333    #nT
#
# Inclination relative to the horizontal
Incl   = -64.68784390   # positive is “down”
#
#  and declination clockwise from true north
Decl   =  0.86171459   #
# gambit configuration file for magnetic inversion 
# for a real data set and mesh file generated 
#
# File names
# mesh file name.  This needs to be in msh or fly format.  
# Magnetic inversion uses same mesh as gravity inversion.
mesh_name = "mesh_51x85.msh"
#
# data file name in netcdf format.  See readme.md for more details.
data_file = "data/Magnetic_51x85.nc"
#
# output file name for .csv output and silo output
output_name = "silos/M_51x85_k0_{0:1.3e}_mu_{1:1.3e}".format(k_0,mu)
csv_name = "csv/M_51x85_k0_{0:1.3e}_mu_{1:1.3e}".format(k_0,mu)
geofile ='mesh_51x85.geo'
meshfile='mesh_51x85.msh'



# It is assumed that the XY-data arrays are flat and parallel to the surface at a given height and
# corresponding data do not change vertically across a thin layer. 
#
# .... these are the horizontal coordinates of the lower-left (south-west) end of the data array in the mesh [m]:
DataRefX=0.0
DataRefY=0.0

# ... this is the height of the grav data above ground [m] (can be zero)
DataHeightAboveGround = 0

# ... number of data points in east-west (X) and north-south (Y) direction:
DataNumX = 51
DataNumY = 85

# .... this spacing of the data array [m]:
DataSpacingX = 3200.
DataSpacingY = 3200.

# Note: the resolution specified here should roughly match the resolution of the actual data as input data are interpolated to the resolution in the mesh

# ... this is the "thickness" of the data array = the thickness of the vertical layer. 
DataMeshSizeVertical = 3200

# ... this is the thickness of region below the data area. In essence it defines the depth of the inversion
CoreThickness = 50000

# ... these are factors by which the DataMeshSizeVertical is raised in the air layer and in the core. 

MeshSizeAirFactor = 10
MeshSizeCoreFactor = 5

#
#
# Level for the verbosity of the output, "low", "medium" or "high".
# low: 
#   screen outputs:
#      data range, 
#      summaries of magnetic data and final magnetism
#      initial, final and difference misfits
#   file output:
#      silo of final solution
# medium: low outputs + 
#   screen outputs:
#      residual norm from the IPCG iterations
# high: medium outputs + 
#   screen outputs:
#      misfit and smoothing value at each iteration step
#   file outputs:
#      csv files for misfit and smoothing at each IPCG iteration
#      silos at misfit values of 0.05, 0.01, 0.008 and 0.005. (Initial misfit is 0.5.)
#VerboseLevel = "low"
#VerboseLevel = "medium"
VerboseLevel = "high"
