############################################################
# Tools for interacting with antechamber
# ~acutil.py~
# created Jan. 2014
# Satoshi Takahama (satoshi.takahama@epfl.ch)
#
# included: 
#   def run_antechamber
#   def convert_filetype
#   def read_ac
############################################################

###_* load libraries
import os
import sys
import subprocess
import pandas as pd
from collections import OrderedDict

###_* define functions
def run_antechamber(execpath,filename,printonly=False):
    ## expecpath: path name
    ## filename : name of file
    ## printonly: show only commands
    ## see http://ambermd.org/antechamber/ac.html
    ## for future: call convert_filetype depending on 
    ##   what is accepted by antechamber

    ## local variables
    name, inp = os.path.splitext(filename)
    inp = inp[1:]
    extensions = ['ac','prepi','frcmod']
    files = dict(map(lambda x: (x,'%s.%s'% (name,x)),extensions))
    files.update({'ext':inp,'orig':filename})
    cmds = ['antechamber -fi %(ext)s -fo ac -i %(orig)s -o %(ac)s' % files,
            'am1bcc -i %(ac)s -f ac -o %(ac)s -j 4' % files,
            'atomtype -i %(ac)s -o %(ac)s -f ac -p gaff' % files,
            'prepgen -i %(ac)s -o %(prepi)s -f car ' % files,
            'parmchk -i %(prepi)s -o %(frcmod)s -f prepi' % files]
    ## print only
    if printonly:
        for line in cmds:
            print execpath+line
        return
    ## execute
    for line in cmds:
        print execpath+line
        try:
            subprocess.call(execpath+line,shell=True)
        except:
            pass

def convert_filetype(origin,target):
    ## origin, target: full filenames of 'from' and 'to'
    ## calls openbabel
    subprocess.call('babel -i %s %s -o %s %s' %
                    (os.path.splitext(origin)[0][1:],origin,
                     os.path.splitext(target)[0][1:],target),
                    shell=True)

def read_ac(filename):
    ## filename: character string
    def set_dtypes(x):
        ## x: 'atoms' or 'bonds'
        ## lexically scoped: table, index, columns, dtypes
        dtype = dtypes[x]
        dframe = pd.DataFrame(table[x],index=index[x],columns=columns[x])
        for i,c in enumerate(dframe.columns):
            dframe[c] = dframe[c].astype(dtype[i])
        return dframe
    ## define local variables
    columns = {'atoms':['atom','subname','subid','x','y','z','charge','type'],
               'bonds':['origin','target','type','origin_name','target_name']}
    dtypes = {'atoms':['object']*2+['int']+['float']*4+['object'],
              'bonds':['int']*3+['object']*2}
    header = []
    index = {'atoms':[],'bonds':[]}
    table = {'atoms':[],'bonds':[]}
    ## read file
    with open(filename) as f:
        ## read header
        for i,line in enumerate(f): 
            header.append(line.strip().split())
            if i > 0:
                break
        ## read table
        for line in f: 
            if line[:4]=='ATOM':
                rest = line.strip().split()[1:]
                index['atoms'].append(int(rest[0]))
                table['atoms'].append(rest[1:])
            elif line[:4]=='BOND':
                rest = line.strip().split()[1:]
                index['bonds'].append(int(rest[0]))
                table['bonds'].append(rest[1:])
            else:
                pass
    ## return ordered dictionary of list, dataframe, dataframe
    return OrderedDict([('header',header),
                        ('atoms',set_dtypes('atoms')),
                        ('bonds',set_dtypes('bonds'))])

###_* example usage
# if __name__ == '__main__':
#     os.environ['ACHOME'] = os.path.join(os.environ['HOME'],'Programs/antechamber-1.27')
#     os.environ['AMBERHOME'] = os.path.join(os.environ['HOME'],'Programs/amber11')
#     molecule = sys.argv[-1]
#     run_antechamber(os.path.join(os.environ['ACHOME'],'exe/'),molecule,True)
