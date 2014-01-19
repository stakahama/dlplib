############################################################
# Tools for creating field files
# ~fieldsutil.py~
# created Jan. 2014
# Satoshi Takahama (satoshi.takahama@epfl.ch)
#
# included: 
#   class MolecGraph
############################################################

###_* load libraries
from collections import OrderedDict
from operator import itemgetter

###_* define class MolecGraphs
class MolecGraph:
    ## for user: find_connections(); match_labels()
    
    def __init__(self,atoms,bonds,atomtypes):
        ## atoms (vertices): integers (any sequence obj)
        ## bonds (edges): pairs of integers (list of tuples)
        ## atomtypes: names of atoms (pandas Series) 
        self.index = OrderedDict([('atoms',atoms), ('bonds',bonds)])
        self.labels = OrderedDict([('atoms',atomtypes)])

    @staticmethod        
    def find_angles(V,E):
        ## returns list of tuples
        ## from Takahama and Russell (2011) Appendix
        return [(i,j,k) for i in V for j in V for k in V
                if i < k
                and ((i,j) in E or (j,i) in E) 
                and ((j,k) in E or (k,j) in E)]

    @staticmethod
    def find_dihedrals(V,E):
        ## from Takahama and Russell (2011) Appendix        
        return [(i,j,k,l) for i in V for j in V for k in V for l in V
                if i!=k and j!=l and j < k
                and ((i,j) in E or (j,i) in E)
                and ((j,k) in E or (k,j) in E)
                and ((k,l) in E or (l,k) in E)]

    def find_connections(self):
        ## returns list of tuples        
        ## assigns to self.index
        V = self.index['atoms']
        E = self.index['bonds']
        self.index['angles'] = self.find_angles(V,E)
        self.index['dihedrals'] = self.find_dihedrals(V,E)

    def lookup(self,attr):
        ## returns list of tuples        
        atomtypes = self.labels['atoms']
        out = []
        for x in self.index[attr]:
            out.append(tuple(atomtypes.ix[list(x)].tolist()))
        return out

    def match_labels(self):
        ## assigns to self.labels
        for x in ['bonds','angles','dihedrals']:
            self.labels[x] = self.lookup(x)
        

###_* example usage

# if __name__ = '__main__':
#     molec = acutils.read_ac(acfile)
#     mgraph = MolecGraph(atoms=molec['atoms'].index,
#                         bonds=map(tuple,molec['bonds'][['origin','target']].values),
#                         atomtypes=molec['atoms']['type'])
#     mgraph.find_connections()
#     mgraph.match_labels()

