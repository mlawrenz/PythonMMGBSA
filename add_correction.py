import numpy, pylab
import os
import sys
import operator
import optparse
from scipy import stats


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''Add precalculated global correcion to data''', formatter_class=ArgumentDefaultsHelpFormatter)

#always need these user inputs
group = parser.add_argument_group('Necessary user input',)
group.add_argument('-data','--data',dest='data',help='data to be corrected')
group.add_argument('-g','--globalcorr',dest='globalcorr',help='correction list')



def get_ref(refdata):
    ref=dict()
    refhandle=open(refdata)
    for line in refhandle.readlines():
        if 'name' in line:
            continue
        else:
            name=line.split()[0]
            k=float(line.split()[1])
            ref[name]=k
    sorted_ref=sorted(ref.iteritems(), key=operator.itemgetter(1))
    return sorted_ref

def check_data(data1, data2, names=False):
    try:
        int(data1[0])
    except ValueError:
        print "removing title %s" % data1[0]
        data1=data1[1:]
    if names==True:
        for i in data2:
            if i not in data1:
                print "missing %s" % i
                print "WARNING: different ligand names, ensure same ligands in files"
            elif len(data1) != len(data2):
                print "WARNING: different # of ligands, ensure same ligands in files"
            else:
                pass
    return data1


def read_file(file, column):
    if column==0:
        print "ASSUMING LIGAND NAME IN COLUMN 0"
    fhandle=open(file)
    array=[]
    for line in fhandle.readlines():
        if column==0:
            array.append(line.split()[column].split('.')[0])
        else:
            try:
                array.append(float(line.split()[column]))
            except ValueError:
                pass
    return numpy.array(array)
    

def main(args):
    data=args.data
    sorted_ref=get_ref(data)
    ref_values=numpy.array([item[1] for item in sorted_ref])
    ligand_names=[item[0] for item in sorted_ref]

    data=args.globalcorr
    correction=dict()
    correction['names']=read_file(data, 0)
    correction['names']=check_data(correction['names'], ligand_names, names=True)
    correction['values']=read_file(data, 1)
    correction['values']=check_data(correction['values'], ligand_names, names=True)
    correction['sorted']=numpy.zeros((len(data)))
    refdata=numpy.zeros(len(data))

    count=0
    for item in sorted_ref:
        name=item[0]
        location=numpy.where(correction['names']==name)[0]
        if len(location) > 1:
            print "HAVE 2 ENTRIES FOR %s" % name
            sys.exit()
        if location.size:
            correction['sorted'][count]=correction['values'][location]
            refdata[count]=item[1]
            count+=1
        else:
            print "no correction for: ", item
            correction['sorted'][count]=0
            count+=1
    location=numpy.where(data==0)[0]
    if location.size:
        print "MISMATCH IN DATA"
        sys.exit()
    print ligand_names
    print ref_values
    print correction['sorted']
    ohandle=open('corr-%s' % args.data, 'w')
    for (lig, val, cor) in zip(ligand_names, ref_values, correction['sorted']):
        total=val+cor
        ohandle.write('%s\t%s\n' % (lig, total))
    ohandle.close()


if __name__=="__main__":	
    args = parser.parse_args()
    main(args)
