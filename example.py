
###_* load libraries

import os
import sys

libpath = os.path.expanduser('~/git/projects/dlplib/src')
if libpath not in sys.path:
    sys.path.append(libpath)
    
from fieldsutil import MolecGraph
from gaffutil import GAFFparms
from acutil import read_ac
from operator import itemgetter
import pandas as pd

###_* define variables

acfile = 'CID_176.ac' # this is adipic acid

###_* read file

molec = read_ac(acfile)

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

directory = os.path.splitext(acfile)[0]
if not os.path.exists(directory):
    os.mkdir(directory)

for x in ['atoms','bonds','angles','dihedrals']:
    gp.extracted[x].to_csv(os.path.join(directory,x+'.dat'),
                           sep='\t',index=False)

