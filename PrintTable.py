import optparse
import os
import glob
import numpy

def print_table(dir, errorwrite=False):
    all_errors=['MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB', 'E_surf', 'E_lig']
    all_values=dict()
    all_errors=dict()
    files=glob.glob('%s/*-cpx-*FINAL*' % dir)
    allcomponents=['']
    for file in files:
        base=os.path.basename(file)
        ligand=base.split('-cpx-')[0]
        if ligand not in all_values.keys():
            all_values[ligand]=dict()
            all_errors[ligand]=dict()
        collect=False
        print "on %s" % file
        fhandle=open(file)
        for line in fhandle.readlines():
            if 'Differences' in line:
                collect=True
            if collect==True:
                if 'DELTA TOTAL' in line:
                    all_values[ligand]['MMGB']=float(line.split()[2])
                    all_errors[ligand]['MMGB']=float(line.split()[3])
                if 'VDWAALS' in line:
                    all_values[ligand]['vdW']=float(line.split()[1])
                    all_errors[ligand]['vdW']=float(line.split()[2])
                if 'EEL' in line and '1-4 EEL' not in line:
                    all_values[ligand]['eel_inter']=float(line.split()[1])
                    all_errors[ligand]['eel_inter']=float(line.split()[2])
                if 'EGB' in line:
                    all_values[ligand]['EGB']=float(line.split()[1])
                    all_errors[ligand]['EGB']=float(line.split()[2])
                if 'ESURF' in line:
                    all_values[ligand]['E_surf']=float(line.split()[1])
                    all_errors[ligand]['E_surf']=float(line.split()[2])
        all_values[ligand]['eel/EGB']=all_values[ligand]['eel_inter']+all_values[ligand]['EGB']
        all_errors[ligand]['eel/EGB']=numpy.sqrt(all_errors[ligand]['eel_inter']**2 + all_errors[ligand]['EGB']**2)
    states=['ligcpx', 'ligsolv'] # cpx - solv = strain
    components=dict()
    for ligandstate in states:
        files=glob.glob('%s/*-%s-*FINAL*' % (dir, ligandstate))
        values=[]
        errors=[]
        components[ligandstate]=dict()
        for file in files:
            base=os.path.basename(file)
            ligand=base.split('-%s-' % ligandstate)[0]
            components[ligandstate][ligand]=dict()
            fhandle=open(file)
            for line in fhandle.readlines():
                if 'TOTAL' in line:
                    components[ligandstate][ligand]['value']=float(line.split()[1])
                    components[ligandstate][ligand]['err']=float(line.split()[2])
                if ligandstate=='ligsolv' and 'G gas' in line:
                    all_values[ligand]['E_lig']=line.split()[2]
                    all_errors[ligand]['E_lig']=line.split()[3]
    for ligand in components['ligcpx'].keys():
        all_values[ligand]['strain']=components['ligcpx'][ligand]['value']-components['ligsolv'][ligand]['value']
        all_errors[ligand]['strain']=numpy.sqrt(components['ligcpx'][ligand]['err']**2+components['ligsolv'][ligand]['err']**2)
    ligands=all_values.keys()
    sorted_ligands=sorted(ligands, key=lambda x: all_values[x]['MMGB'])
    ohandle=open('%s/sorted_results.tbl' % dir, 'w')
    formatkeyorder=['name  MMGB+str', 'MMGB  strain  vdW  eel_inter  eel/EGB  EGB  E_surf  E_lig'] 
    entry=''.join(['%s\t' % x for x in formatkeyorder])
    entry=''.join([ entry, '\n'])
    ohandle.write(entry)
    print entry
    keyorder=['MMGB+str', 'MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB' , 'E_surf',  'E_lig'] 
    for ligand in sorted_ligands:
        all_values[ligand]['MMGB+str']=all_values[ligand]['MMGB']+all_values[ligand]['strain']
        name='%s\t\t' % ligand
        entry=''.join(['%0.2f\t' % round(float(all_values[ligand][x]), 2) for x in keyorder])
        entry=''.join([name, entry, '\n'])
        ohandle.write(entry)
        print entry
    ohandle.close()
    if errorwrite==True:
        ohandle=open('%s/estimates_error.txt' % dir, 'w')
        header='name     MMGB+str   err    MMGB   err    strain    err\n'
        ohandle.write(header)
        print header
        for ligand in sorted_ligands:
            comboerr=numpy.sqrt(all_errors[ligand]['strain']**2+all_errors[ligand]['MMGB']**2)
            line='%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ligand, round(all_values[ligand]['MMGB+str'], 2), round(comboerr,2), round(all_values[ligand]['MMGB'], 2), round(all_errors[ligand]['MMGB'],2), round(all_values[ligand]['strain'],2), round(all_errors[ligand]['strain'], 2))
            ohandle.write(line)
            print line
    ohandle.close()


def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-d','--dir',dest='dir',type='string')
    parser.add_option('-e', action="store_true", dest="error", help="using -e will print errors on estimates")
    (options, args) = parser.parse_args()

    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    if options.error==True:
        print_table(dir=options.dir, errorwrite=True)
    else:
        print_table(dir=options.dir)

