import pandas, numpy, scipy
import random
import statsmodels.formula.api as sm
import pylab
import os
import sys
import operator
import optparse
from scipy import stats


def get_pdb_coors(file):
    tmp=dict()
    fhandle=open(file)
    for line in fhandle.readlines():
        type=line.split()[0]
        if type=='ATOM': # parse normal PDB line
            atom=line.split()[2]
            resname=line.split()[3]
            x=float(line.split()[6])
            y=float(line.split()[7])
            z=float(line.split()[8])
            try:
                int(line.split()[4])
                pdbnum=int(line.split()[4]) # a check for columns
            except ValueError:
                print "PDB contains chain column"
                pdbnum=int(line.split()[5]) # a check for columns
            if pdbnum not in tmp.keys():
                tmp[pdbnum]=[]
            tmp[pdbnum].append(numpy.array([x,y,z]))
    return tmp
                

def get_distance(a, b):
    dist=numpy.linalg.norm(a-b)
    return dist

def get_mol2_coors(file):
    tmp=[]
    read=False
    for line in file.readlines():
        if 'ATOM' in line:
            read=True
            continue
        if 'BOND' in line:
            break
        if read==True:
            x=float(line.split()[2])
            y=float(line.split()[3])
            z=float(line.split()[4])
            tmp.append(numpy.array([x,y,z]))
    return tmp


def main(protfile, liglist):
    coors=dict()
    coors['prot']=get_pdb_coors(protfile)
    coors['lig']=dict()
    list_handle=open(liglist)
    for ligfile in list_handle.readlines():
        ligfile=ligfile.rstrip('\n')
        name=os.path.basename(ligfile).split('.mol2')[0]
        lighandle=open(ligfile)
        coors['lig'][name]=get_mol2_coors(lighandle)
    pocket=[]
    restrain=[]
    for res in coors['prot'].keys():
        for a in coors['prot'][res]:
            for ligand in coors['lig'].keys():
                for b in coors['lig'][ligand]:
                    dist=get_distance(a,b)
                    if dist <= 5.0:
                        if res not in pocket:
                            print "adding pocket %s" % res
                            pocket.append(res)
                        else:
                            continue
                    else:
                        if res not in restrain:
                            restrain.append(res)
                        else:
                            continue
    print pocket 
    new_restrain=[]
    for i in restrain:
        if i not in pocket:
            new_restrain.append(i)
    numpy.savetxt('pocket_consensus.txt', pocket, fmt='%i')
    numpy.savetxt('restrain_consensus.txt', new_restrain, fmt='%i')
    
def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-p','--protfile',dest='protfile',type='string')
    parser.add_option('-l','--liglist',dest='liglist',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(protfile=options.protfile, liglist=options.liglist)

