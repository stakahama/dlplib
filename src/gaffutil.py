############################################################
# Tools for working with GAFF parameter tables
# ~acutil.py~
# created Jan. 2014
# Satoshi Takahama (satoshi.takahama@epfl.ch)
#
# included:
#   class GAFFparms
############################################################

###_* load libraries
import re
from collections import OrderedDict
from operator import itemgetter
import pandas as pd

###_* define class GAFFparms
class GAFFparms:

    def __init__(self):
        pass

###_ . read gaff.dat
    def read(self,filename):
        # main function; calls self.make_tables()
        main = []
        main.append([])
        with open(filename) as f:
            ix = 0
            for line in f:
                if line == '\n':
                    main.append([])
                    ix += 1
                    continue
                main[ix].append(line.strip())
        self.main = main
        self.make_tables()

    def _asdframe(self,table,skip=0,patt=None,names=None):
        # called by self.make_tables()
        text = table[skip]
        n = -1
        if patt:
            lastm = [m for m in re.finditer(patt, text, re.IGNORECASE)][n]
            lastpos = lastm.start()+1
        else:
            lastpos = min(map(len,table[skip:]))+1
        index = []
        fields = []
        for line in table[skip:]:
            fd = line[:lastpos].replace(' -','-').split()
            ix = tuple(map(str.strip,fd[0].split('-'))) if '-' in fd[0] else fd[0]
            index.append(ix)
            fields.append(fd[1:])
        return pd.DataFrame(fields,index=index,columns=names)

    def make_tables(self):
        # calls self._asdframe()
        # called by self.read()
        main = self.main
        defpatt = '[ ]{2,}'
        tables = OrderedDict()
        tables['atoms'] = self._asdframe(main[0],1,defpatt,('mass','unknown'))
        tables['bonds'] = self._asdframe(main[1],1,defpatt,('k','r'))
        tables['angles'] = self._asdframe(main[2],0,defpatt,('k','theta'))
        tables['imprdihedrals'] = self._asdframe(main[3],0,None,('int','A','gamma','m'))
        tables['propdihedrals'] = self._asdframe(main[4],0,None,('A','gamma','m'))
        tables['lennjo'] = self._asdframe(main[7],1,defpatt,('Rstar','epsilon'))
        self.parms = tables

###_ . match labels from MolecGraph

    # match each type. test its reverse and also X-substituted cases (for dihedrals)
    @staticmethod
    def test_reverse(index):
        # searches for matches;
        # both straight and its reverse case in key index
        # used by match_*
        def fn(x):
            if x in index:
                match = x
            else:
                rev = tuple(reversed(x))
                match = rev if rev in index else None
            return match
        return fn

    def match_atoms(self):
        return self.keys.labels['atoms'].tolist()

    def match_bonds(self):
        index = self.parms['bonds'].index        
        keys = self.keys.labels['bonds']
        return map(self.test_reverse(index),keys)

    def match_angles(self):
        index = self.parms['angles'].index        
        keys = self.keys.labels['angles']
        return map(self.test_reverse(index),keys)

    def match_imprdihedrals(self):
        index = self.parms['imprdihedrals'].index        
        keys = self.keys.labels['dihedrals']
        fn = self.test_reverse(index)
        out = []
        for x in keys:
            match = fn(x)
            if match:
                out.append(match)
                continue
            match = fn(('X',)+x[1:-1]+('X',))
            if match:
                out.append(match)
                continue
            out.append(None)
        return out

    def match_propdihedrals(self):
        index = self.parms['propdihedrals'].index        
        keys = self.keys.labels['dihedrals']
        fn = self.test_reverse(index)
        out = []    
        for x in keys:
            match = fn(x)
            if match:
                out.append(match)
                continue
            match = fn(('X',)+x[1:])
            if match:
                out.append(match)
                continue
            match = fn(('X','X')+x[2:])
            if match:
                out.append(match)
                continue
            match = fn(x[:-1]+('X',))
            if match:
                out.append(match)
                continue
            match = fn(x[:-2]+('X','X'))
            if match:
                out.append(match)
                continue
            match = fn(('X',)+x[1:-1]+('X',))
            if match:
                out.append(match)
                continue
            out.append(None)
        return out

    def set_keys(self,keys):
        # keys is the object created by MolecGraph class
        self.keys = keys

    def get_matches(self,keys):
        # keys is the object created by MolecGraph class        
        self.keys = keys
        index = OrderedDict()
        for x in ['atoms','bonds','angles','propdihedrals','imprdihedrals']:
            index[x] = getattr(self,'match_'+x)()
        self.index = index

    def find_missing_matches(self,attr,seq=None):
        # diagnostic function
        # attr: character string
        # seq: character string ('atoms','bonds','angles',etc.)
        # if seq is not supplied, just returns indices
        index = [i for i,x in enumerate(getattr(self,attr)) if x is None]
        if len(index) > 0:
            if seq:
                value = itemgetter(*index)(getattr(self.keys,seq))
            else:
                value = index
        else:
            value = None
        return value

###_ . extract parameters from saved tables
    def merge(self,parms,keys,index):
        if type(parms.index[0]) is tuple:
            n = len(parms.index[0])
            replace = lambda n: lambda x: x if x else (pd.np.nan,)*n
            newkeys = map(replace(n),keys)
            atypes = map('atomtype_{:d}'.format,range(n))
            atoms = map('atom_{:d}'.format,range(n))
        else:
            newkeys = keys
            atypes = ['atomtype']
            atoms = ['atom']
        molec = pd.DataFrame(index,index,atoms).\
                join(pd.DataFrame(newkeys,index,atypes))
        gaff = pd.DataFrame(parms.index.tolist(),parms.index,atypes).\
               join(parms)
        # out = pd.merge(molec.reset_index(),gaff,how='left',on=atypes).set_index('index')
        out = pd.merge(molec,gaff,how='left',on=atypes)
        
        return out.drop(atypes,axis=1)

    def extract_parms(self):
        extracted = OrderedDict()
        for x in ['atoms','bonds','angles','propdihedrals','imprdihedrals']:
            x_ = x if 'dihedrals' not in x else 'dihedrals'
            extracted[x] = self.merge(self.parms[x],self.index[x],self.keys.index[x_])
        missing = extracted['propdihedrals'].ix[:,-1].isnull()
        atoms = map('atom_{:d}'.format,range(4))
        if pd.np.any(missing):
            extracted['dihedrals'] = pd.merge(extracted['propdihedrals'].ix[missing,atoms],
                                              extracted['imprdihedrals'],on=atoms,
                                              how='outer')
        else:
            extracted['dihedrals'] = extracted['propdihedrals']
        self.extracted = extracted

###_* example usage
# if __name__ == '__main__':
#     gp = GAFFparms()
#     gp.read('/Users/stakahama/Programs/antechamber-1.27/dat/leap/parm/gaff.dat')
#     gp.get_matches(mgraph)
#     gp.extract_parms()


###_* scratch
############################################################
# alternative GAFFParms.read()

# class GAFFParms:

#     def __init__(self):
#         pass

#     def read(self,filename):

#         def split_string(string,ix):
#             ix = list(ix)+[len(string)]
#             st = 0
#             substr = []
#             for en in ix:
#                 substr.append(string[st:en].strip())
#                 st = en
#             return substr

#         def read_block(f,fw=None,n=None):
#             out = []    
#             if n is None:
#                 for line in f:
#                     line = line.strip()
#                     if len(line)==0:
#                         break
#                     out.append(split_string(line,fw))
#             else:
#                 for i,line in enumerate(f):
#                     line = line.strip()
#                     if len(line)>0:
#                         out.append(line)
#                     if i >= n-1:
#                         break
#             return out

#         def make_dframe(contents,names):
#             fn = lambda x: '-' in x and tuple(map(str.strip,x.split('-'))) or x
#             index = []
#             table = []
#             for x in contents:
#                 index.append(fn(x[0]))
#                 table.append(x[1:][:len(names)])
#             dframe = pd.DataFrame(table,index=index,columns=names)
#             return dframe.convert_objects(convert_numeric=True)
        
#         def fw(n):
#             return 2*n+(n-1)
        
#         with open(filename) as f:
#             misc = []
#             misc += read_block(f,n=1)
#             atoms = read_block(f,[fw(1),8,22])
#             misc += read_block(f,n=1)
#             bonds = read_block(f,[fw(2),12,22,39,42,52])
#             angles = [re.split('[ ]{2,}',x) for x in read_block(f,n=2)]
#             angles += read_block(f,[fw(3),16,28,39,44,53,62])
#             imprdihedrals = read_block(f,[fw(4),15,24,38,54])
#             propdihedrals = read_block(f,[fw(4),15,24,38,54])
#             misc += read_block(f,n=5)
#             ljpars = read_block(f,[fw(1)+2,20,28])
            
#         self.parms = OrderedDict([('atoms',make_dframe(atoms,('mass','unknown'))),
#                                   ('bonds',make_dframe(bonds,('k','r'))),
#                                   ('angles',make_dframe(angles,('k','theta'))),
#                                   ('imprdihedrals',make_dframe(imprdihedrals,('int','A','gamma','m'))),
#                                   ('propdihedrals',make_dframe(propdihedrals,('A','gamma','m'))),
#                                   ('ljpars',make_dframe(ljpars,('Rstar','epsilon'))),
#                                   ('misc',misc)])
