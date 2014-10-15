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

def check_data(data1, data2):
    if len(data1) != len(data2):
        data1=data1[1:]
        if len(data1) != len(data2):
            print "problem with parsing data"
            sys.exit()
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
    

def main(refdata, adata, bdata, cdata=None):
    sorted_ref=get_ref(refdata)
    ref_values=numpy.array([item[1] for item in sorted_ref])
    ligand_names=[item[0] for item in sorted_ref]
    calc=dict()
    if cdata==None:
        namelist=[adata, bdata]
    else:
        namelist=[adata, bdata, cdata]
    for (n, data) in enumerate(namelist):
        calc[n]=dict()
        calc[n]['names']=read_file(data, 0)
        calc[n]['names']=check_data(calc[n]['names'], ref_values)
        calc[n]['values']=read_file(data, 1)
        calc[n]['values']=check_data(calc[n]['values'], ref_values)
        calc[n]['sorted']=numpy.zeros((len(calc[n]['values'])))
        for (ref_location, item) in enumerate(sorted_ref):
            name=item[0]
            location=numpy.where(calc[n]['names']==name)[0]
            calc[n]['sorted'][ref_location]=calc[n]['values'][location]
    colors=['r', 'b', 'k', 'g']
    pylab.figure()
    print "ligand names: ", ligand_names
    print "reference values are: ", ref_values
    for n in sorted(calc.keys()):
        name=namelist[n]
        data=calc[n]['sorted']
        print "--------------"
        print "%s values are: " % name, data
        print "--------------"
        slope, intercept, r_val, p_val, std_err=stats.linregress(numpy.array(ref_values), numpy.array(data))
        print "R^2=%s" % round(r_val**2, 2)
        print "p=%s" % round(p_val**2, 2)
        print "sigma=%s" % round(std_err, 2)
        line=slope*numpy.array(ref_values)+intercept
        #R=pylab.corrcoef(numpy.array(ref_values), numpy.array(data))
        pylab.scatter(numpy.array(ref_values), numpy.array(data), c=colors[n],
label='%s R^2=%s' % (os.path.dirname(name), round(r_val**2, 2)))
        pylab.plot(ref_values, line, '-', c=colors[n]) 
        pylab.hold(True)
        lg=pylab.legend()
        lg.draw_frame(False)
        pylab.ylabel('calc dG (MMGBSA)')
        pylab.xlabel('exp dG (from IC50)')
        n+=1
    pylab.show()


def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-r','--refdata',dest='refdata',type='string')
    parser.add_option('-a','--adata',dest='adata',type='string')
    parser.add_option('-b','--bdata',dest='bdata',type='string')
    parser.add_option('-c','--cdata',dest='cdata',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(refdata=options.refdata, adata=options.adata, bdata=options.bdata, cdata=options.cdata)

