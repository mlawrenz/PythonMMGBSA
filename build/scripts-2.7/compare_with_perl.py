import os
import glob
import numpy
import optparse



def main(dir):
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
                if 'EEL' in line:
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
        all_errors[ligand]['strain']=error=numpy.sqrt(components['ligcpx'][ligand]['err']**2+components['ligsolv'][ligand]['err']**2)
    ligands=all_values.keys()
    sorted_ligands=sorted(ligands, key=lambda x: all_values[x]['MMGB'])
    ohandle=open('%s/sorted_results.tbl' % dir, 'w')
    keyorder=['MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB', 'E_surf', 'E_lig'] 
    entry=''.join(['%s\t' % x for x in keyorder])
    entry=''.join(['name\t', entry, '\n'])
    ohandle.write(entry)
    print entry
    for ligand in sorted_ligands:
        name='%s\t' % ligand
        entry=''.join(['%0.3f\t' % round(float(all_values[ligand][x]), 3) for x in keyorder])
        entry=''.join([name, entry, '\n'])
        ohandle.write(entry)
        print entry


def parse_cmdln():
    parser=optparse.OptionParser()
    parser.add_option('-d','--dir',dest='dir',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(dir=options.dir)

