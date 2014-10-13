import os

class OptionList(object):
   """
   Just a container to hold the command-line options. 
   """
   pass

# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods, using
                        minimized snapshots of MD simulation snapshots run in 
                        implicit or explicit water. The free energy calculations 
                        themselves are evaluated with GB implicit solvation term.'''
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', action='version',
                    version='%%(prog)s %s' % __version__)
parser.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                    help='Print all available options in the input file.',
                    default=False)
group = parser.add_argument_group('Miscellaneous Options')
group.add_argument('-O', '--overwrite', default=False, action='store_true',
                  help='Allow output files to be overwritten', dest='overwrite')
group.add_argument('-prefix', dest='prefix', default='_MMPBSA_',
                  metavar='<file prefix>',
                  help='Prefix for intermediate files.')
group = parser.add_argument_group('Input and Output Files', '''These options
                        specify the input files and optional output files.''')
group.add_argument('-i', dest='input_file', metavar='FILE',
                   help='MM/PBSA input file.')
group.add_argument('-xvvfile', dest='xvvfile', help='XVV file for 3D-RISM.',
                  default=os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa',
                                       'spc.xvv'))
group.add_argument('-o', dest='output_file', default='FINAL_RESULTS_MMPBSA.dat',
                  metavar='FILE', help='Output file with MM/PBSA statistics.')
group.add_argument('-do', dest='decompout', metavar='FILE',
                   default='FINAL_DECOMP_MMPBSA.dat',
                   help='Output file for decomposition statistics summary.')
