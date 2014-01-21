
###_* load libraries

import os
import sys

libpath = os.path.expanduser('~/git/projects/dlplib/src')
if libpath not in sys.path:
    sys.path.append(libpath)
    
from fieldutil import MolecGraph
from gaffutil import GAFFparms
from acutil import convert_filetype, run_antechamber, read_ac
from operator import itemgetter
import pandas as pd

###_* define variables

###_ . file name
molecfile = 'CID_176' # this is adipic acid

###_ . environment variables
os.environ['ACHOME'] = os.path.join(os.environ['HOME'],'Programs/antechamber-1.27')
os.environ['AMBERHOME'] = os.path.join(os.environ['HOME'],'Programs/amber11')

###_* run antechamber

convert_filetype(molecfile+'.sdf',molecfile+'.pdb')
run_antechamber(os.path.join(os.environ['ACHOME'],'exe/'),molecfile)

###_* read antechamber file

molec = read_ac(molecfile+'.ac')

###_* find angles and dihedrals from vertices and edges

mgraph = MolecGraph(atoms=molec['atoms'].index,
                    bonds=map(tuple,molec['bonds'][['origin','target']].values),
                    atomtypes=molec['atoms']['type'])
mgraph.find_connections()
mgraph.match_labels()

###_* read and match GAFF parameters

gp = GAFFparms()
gp.read('/Users/stakahama/Programs/antechamber-1.27/dat/leap/parm/gaff.dat')
gp.get_matches(mgraph)
gp.extract_parms()

###_* export tables corresponding to molecule

directory = molecfile#os.path.splitext(molecfile)[0]
if not os.path.exists(directory):
    os.mkdir(directory)

for x in ['atoms','bonds','angles','dihedrals']:
    gp.extracted[x].to_csv(os.path.join(directory,x+'.dat'),
                           sep='\t',index=False)

###_* calculate Lennard Jones interactions and export table

water = ['ow','hw']
allatoms = gp.index['atoms'] + water

gp.lennjopairs(allatoms)
gp.ljparms.to_csv('lennjo.dat',sep='\t',index=False)
