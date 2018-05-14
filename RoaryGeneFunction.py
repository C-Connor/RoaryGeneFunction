import os
import sys
import csv
from Bio import SeqIO
import tempfile
import pipes
import time
import readline
import glob
import argparse

#parse command line arguments
cl_parser = argparse.ArgumentParser(description="Calculates gene frequency in population from Roary matrix. Optionally assigns function.")
cl_parser.add_argument("-c","--csv", help="Path to gene_presence_absence.csv file produced by roary",dest="csv",type=str,required=True)
cl_parser.add_argument("-o","--output",help="Specify output path",default=os.getcwd()+os.sep,dest="output",type=str)
cl_parser.add_argument("-a","--assign_func",help="Assign generic cell function to genes (may prompt for manual input)",dest="func",action="store_true")
cl_parser.add_argument("-l","--function_list",help="Path to previously created gene function list, use with -a/--assign_func",dest="list",type=str)

cl_args = cl_parser.parse_args()

#check user input
if not os.path.isfile(cl_args.csv):
    print('Could not find gene_presence_absence.csv file.')
    print('Exiting')
    sys.exit()

if cl_args.func and cl_args.list != None and not os.path.isfile(cl_args.list):
    print('Could not find previously created gene function list')
    print('Omitting -l/--function_list will result in creation of blank function list.')
    print('Exiting')
    sys.exit()

if not os.path.exists(cl_args.output):
    os.makedirs(cl_args.output)

if cl_args.func and cl_args.list == None:
    tempfunctionlist = tempfile.NamedTemporaryFile(suffix = '.txt', prefix = 'TempFunctionList', dir=cl_args.output, bufsize=0, delete = True)
    tempfunctionlist.write('hypothetical protein\tHypothetical Protein\n')
    tempfunctionlist.seek(0)

    
#function dictionary
if cl_args.func:
    functiondict = dict()
    if cl_args.list == None:
        reader = csv.reader(tempfunctionlist, delimiter = '\t')
    else:
        indict = open(cl_args.list,'r')
        reader = csv.reader(indict, delimiter = '\t')
    for tabline in reader:
        functiondict[tabline[0]] = [tabline[1]][0]
    if cl_args.list == None:
        tempfunctionlist.close()
    else:
        indict.close()

#create function response dictionary
responseoptions = {'ad':'Adhesion',
                   'adhesion':'Adehsion',
                   'ar':'Anaerobic Respiration',
                   'anaerobic respiration':'Anaerobic Respiration',
                   'ab':'Antimicrobial Resistance',
                   'antimicrobial resistance':'Antimicrobial Resistance',
                   'cd':'Cell Division',
                   'cell division':'Cell Division',
                   'cw':'Cell Wall',
                   'cell wall':'Cell Wall',
                   'co':'Conjugal Transfer',
                   'conjugal transfer':'Conjugal Transfer',
                   'ei':'Environmental Interaction',
                   'environmental interaction':'Environmental Interaction',
                   'fl':'Flagella',
                   'flagella':'Flagella',
                   'gp':'Genetic Processing',
                   'genetic processing':'Genetic Processing',
                   'im':'Inner Membrane Proteins',
                   'inner membrane proteins':'Inner Membrane Proteins',
                   'ir':'Iron Acquisition',
                   'iron acquisition':'Iron Acquisition',
                   'me':'Metabolism',
                   'metabolism':'Metabolism',
                   'mg':'Mobile Genetic Elements',
                   'mobile genetic elements':'Mobile Genetic Elements',
                   'tr':'Molecular Transport',
                   'molecular transport':'Molecular Transport',
                   'om':'Outer Membrane Proteins',
                   'outer membrane proteins':'Outer Membrane Proteins',
                   'ph':'Phage',
                   'phage':'Phage',
                   'pl':'Plasmid',
                   'plasmid':'Plasmid',
                   'pu':'Protein of Unknown Function',
                   'protein of unknown function':'Protein of Unknown Function',
                   'ri':'Ribosome',
                   'ribosome':'Ribosome',
                   'vi':'Virulence',
                   'virulence':'Virulence',
                   'misc':'Miscellaneous',
                   'hyp':'Hypothetical Protein'}

#create storage lists for gene names and functions
outputcontents = list()

#read in csv file and calculate gene presence percentage
#output gene name and function to lists
with open(cl_args.csv, 'r') as csvin:
    reader = csv.reader(csvin, delimiter=',')
    reader.next() #skip header line
    for tabline in reader: #file contains 14 columns before gene presence/absence, presence is marked by genome name: convert to True with bool then convert to 1 with float
        tabline[14:len(tabline)] = [float(bool(x)) for x in tabline[14:len(tabline)]]
        presence_percentage = float(sum(tabline[14:len(tabline)])) / float(len(tabline)-14)
        if cl_args.func:
            if not tabline[2] in functiondict.keys():
                print('Please assign gene annotation to cellular process.')
                print('Case Insensitive Options: (for full list of keys, type FullKeys)')
                print('')
                print('AD, Adhesion')
                print('AR, Anaerobic Respiration')
                print('AB, Antimicrobial Resistance')
                print('CD, Cell Division')
                print('CO, Conjugal Transfer')
                print('CW, Cell Wall')
                print('EI, Enivronmental Interaction')
                print('FL, Flagella')
                print('GP, Genetic Processing')
                print('IM, Inner Membrane Proteins')
                print('IR, Iron Acquisition')
                print('ME, Metabolism')
                print('MG, Mobile Genetic Elements')
                print('TR, Molecular Transport')
                print('OM, Outer Membrane Proteins')
                print('PH, Phage')
                print('PL, Plasmid')
                print('PU, Protein of Unknown Function')
                print('RI, Ribosome')
                print('VI, Virulence')
                print('MISC, Miscellaneous')
                print('HYP, Hypothetical Protein')
                print('')
                #print('If unsure type BLAST to run blast search with nucleotide sequence')
                print(tabline[2])
                response = pipes.quote(raw_input('Enter option: '))
                while True:
                    if response.lower() in responseoptions.keys():
                        functiondict[tabline[2]] = [responseoptions.get(response.lower())][0]
                        break
                    elif response == 'FullKeys':
                        for i in sorted(responseoptions.items()):
                            print i
                    #elif response == 'BLAST':
                        #somehow run blast maybe
                    print('')
                    print(tabline[2])
                    response = pipes.quote(raw_input('Please enter option: '))
            outputcontents.append([tabline[0],tabline[2],functiondict.get(tabline[2]),str(presence_percentage)])
        else:
            outputcontents.append([tabline[0],tabline[2],str(presence_percentage)])

#output file paths
outputfile = os.path.join(cl_args.output, 'Gene_Frequency_'+time.strftime('%Y%m%d_%H%M')+'.txt')

#write out gene names and functions to file handles as tab separated values
with open(outputfile,'w') as out_file:
    out_file.writelines('\t'.join(i)+'\n' for i in outputcontents)

#write out new function dictionary
if cl_args.func:
    outputdict = os.path.join(cl_args.output,'Function_List_'+time.strftime('%Y%m%d_%H%M')+'.txt')
    with open(outputdict,'w') as out_dict:
        for prokka, function in sorted(functiondict.items()):
            out_dict.write('%s\t%s\n' % (prokka, function))
print('Complete')
print('Ouput in: '+cl_args.output)
