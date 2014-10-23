import numpy, pylab
import os
import sys
import operator
import optparse
from scipy import stats


def get_ref(refdata):
    ref=dict()
    refhandle=open(refdata)
    for line in refhandle.readlines():
        if 'root' in line:
            continue
        else:
            name=line.split()[0]
            k=float(line.split()[1])
            ref[name]=round(0.6*numpy.log(k*10**(-9)), 2)
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
    

def main(refdata, adata, bdata, cdata=None, ddata=None, output=False):
    sorted_ref=get_ref(refdata)
    ref_values=numpy.array([item[1] for item in sorted_ref])
    ligand_names=[item[0] for item in sorted_ref]
    calc=dict()
    if bdata==None:
        namelist=[adata, ]
    elif cdata==None:
        namelist=[adata, bdata]
    elif ddata==None:
        namelist=[adata, bdata, cdata]
    else:
        namelist=[adata, bdata, cdata, ddata]
    for (n, data) in enumerate(namelist):
        print "gathering and checking data from %s" % data
        calc[n]=dict()
        calc[n]['names']=read_file(data, 0)
        calc[n]['names']=check_data(calc[n]['names'], ligand_names, names=True)
        calc[n]['values']=read_file(data, 1)
        calc[n]['values']=check_data(calc[n]['values'], ref_values)
        calc[n]['sorted']=numpy.zeros((len(calc[n]['values'])))
    colors=['r', 'b', 'k', 'g', 'm']
    pylab.figure()
    print "ligand names: ", ligand_names
    print "reference values are: ", ref_values
    if output==True:
        ohandle=open('statistics.output', 'w')
    for n in sorted(calc.keys()):
        dataname=namelist[n]
        data=numpy.zeros((len(calc[n]['values'])))
        refdata=numpy.zeros((len(calc[n]['values'])))
        count=0
        for item in sorted_ref:
            name=item[0]
            location=numpy.where(calc[n]['names']==name)[0]
            if location.size:
                data[count]=calc[n]['values'][location]
                refdata[count]=item[1]
                count+=1
        print "--------------"
        print "%s values are: " % dataname, data
        print "--------------"
        slope, intercept, r_val, p_val, std_err=stats.linregress(data, refdata)
        print "R^2=%s" % round(r_val**2, 2)
        print "p=%s" % round(p_val, 4)
        print "coeff error =%s" % round(std_err, 2)
        line=slope*numpy.array(data)+intercept
        residuals=[(i-j)**2 for (i,j) in zip(refdata, line)]
        newerr=numpy.sqrt(sum(residuals)/(len(residuals)-2))
        print "residual std. err: %s" % round(newerr, 2)
        if output==True:
            ohandle.write('%s\t%s\t%s\t%s\n' % (namelist[n], round(r_val**2, 2), round(p_val, 4), round(newerr, 2)))
        pylab.scatter(data, refdata, c=colors[n], label='%s R^2=%s, N=%s' % (os.path.basename(dataname), round(r_val**2, 2), len(data)))
        pylab.plot(data, line, '-', c=colors[n]) 
        pylab.hold(True)
        lg=pylab.legend(loc=2)
        lg.draw_frame(False)
        pylab.xlabel('calc dG (MMGBSA)')
        pylab.ylabel('exp dG (from IC50)')
        n+=1
    if output==True:
        ohandle.close()
        pylab.savefig('plots.png', dpi=300)
    pylab.show()


def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-r','--refdata',dest='refdata',type='string')
    parser.add_option('-a','--adata',dest='adata',type='string')
    parser.add_option('-b','--bdata',dest='bdata',type='string')
    parser.add_option('-c','--cdata',dest='cdata',type='string')
    parser.add_option('-d','--ddata',dest='ddata',type='string')
    parser.add_option('-o', action="store_true", dest="output", help="using -o will save statistics to file statistics.output")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    if options.output==True:
        main(refdata=options.refdata, adata=options.adata, bdata=options.bdata, cdata=options.cdata, ddata=options.ddata, output=True)
    else:
        main(refdata=options.refdata, adata=options.adata, bdata=options.bdata, cdata=options.cdata, ddata=options.ddata)

