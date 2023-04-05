#!/usr/bin/python3
#import esys.escript.unitsSI as U
import numpy as np
import importlib, sys, os
sys.path.insert(0, os.getcwd())
import json
import argparse
parser = argparse.ArgumentParser(description='this generates gmsh geo file for gravity or magnetic application', epilog="version 11/2020 by l.gross@uq.edu.au")
parser.add_argument(dest='config', metavar='CONFIG', type=str, help='configuration file.')

print("** Generates a geo file for gravity and/or magnetic data file generation **")
args = parser.parse_args()

config = importlib.import_module(args.config)
print("configuration "+args.config+".py imported.")
if config.DataHeightAboveGround >0:
    GEOTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "geoTemplates/AboveGroundTemplate.geo")
else:
    GEOTEMPLATE=os.path.join(os.path.dirname(os.path.abspath(__file__)), "geoTemplates/OnGroundTemplate.geo")

DataSpacingX=config.DataSpacingX
DataSpacingY=config.DataSpacingY

mapping= { 
"DataRefX" : config.DataRefX, 
"DataRefY" : config.DataRefY,
"DataHeightAboveGround" : config.DataHeightAboveGround,
"DataSpacingX" : DataSpacingX,
"DataSpacingY" : DataSpacingY,
"DataNumX" : config.DataNumX,
"DataNumY" : config.DataNumY,
"DataMeshSizeVertical" : config.DataMeshSizeVertical,
"CoreThickness" : config.CoreThickness,
"AirLayerThickness" : config.AirLayerThickness,
"PaddingX" : config.PaddingX,
"PaddingY" : config.PaddingY,
"PaddingZ" : config.PaddingZ,
"PaddingAir" : config.PaddingAir,
"MeshSizeAirFactor" : config.MeshSizeAirFactor,
"MeshSizeCoreFactor" : config.MeshSizeCoreFactor,
"MeshSizePaddingFactor" : config.MeshSizePaddingFactor
}
print("Setting: ")
print(json.dumps(mapping, sort_keys=True, indent=4)[2:-2])

if config.geofile:
    GEOFN=config.geofile
else:
    GEOFN=config.project+".geo"
text=open(GEOTEMPLATE,'r').read().format(**mapping)
open(GEOFN, "w").write(text)
print("GMSH geofile has been written to ",GEOFN)
print("to generate mesh run:")
print("     gmsh -3 -format msh2 -o %s %s"%(config.meshfile, GEOFN))


