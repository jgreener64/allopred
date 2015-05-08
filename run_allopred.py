#
#    AlloPred - Predict allosteric pockets on proteins
#    Copyright (C) 2014 Joe Greener
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


"""AlloPred - Predict allosteric pockets on proteins
Copyright (C) 2014 Joe Greener

Usage:
python $ALLOPRED_DIR/run_allopred.py in_file act_res.txt

Arguments are the input file prefix and the active site residue file.
The working directory must contain the input PDB file (in_file.pdb) and the fpocket output directory (in_file_out).
The environmental variables $ALLOPRED_DIR and $SVM_LIGHT_DIR have to be defined - they refer to the AlloPred and SVM-light directories respectively.
See the AlloPred README for more information.
"""


# PARAMETERS
#
# these can be adjusted, but bear in mind the SVM is optimised on the default parameters
#
# proportion to increase the spring constant by at a chosen pocket - default is 1.5
# this parameter represents the reduction in flexibility of a pocket on modulator binding
spring_frac = 1.5
# atoms to use for normal mode calculation - default is 'calpha' but other options are 'backbone' and 'all'
# choosing 'backbone' will increase computation time; choosing 'all' will increase computation time significantly
nma_atoms = 'calpha'
# cutoff in Angstrom for ANM normal mode calculation; default is 15
nma_cutoff = 15


# check whether to print instructions and exit
import sys
if len(sys.argv) != 3:
    print __doc__
    sys.exit()


# check whether environmental variables are set
import os
try:
    allopred_dir = os.environ['ALLOPRED_DIR']
except KeyError:
    print 'Environmental variable ALLOPRED_DIR cannot be found - see the AlloPred README'
    sys.exit()
try:    
    svm_light_dir = os.environ['SVM_LIGHT_DIR']
except KeyError:
    print 'Environmental variable SVM_LIGHT_DIR cannot be found - see the AlloPred README'
    sys.exit()


import re
from operator import itemgetter
from prody import *


# import custom functions
from functions.pocket_methods import *
from functions.nma_methods import *


# read arguments
input_prefix, act_res_filepath = sys.argv[1:]
pdb_file = input_prefix+'.pdb'
pocket_folder = input_prefix+'_out'
print '- Running AlloPred -'
print 'PDB file:', pdb_file
print 'Active site residue file:', act_res_filepath
print 'Pocket directory:', pocket_folder


# read PDB file
print 'Reading PDB file'
pro = parsePDB(pdb_file)
calphas = pro.select('calpha')
nma_selection = pro.select(nma_atoms)


# read active site residues from file
print 'Reading active site residue file'
act_res_file = open(act_res_filepath)
act_res_lines = act_res_file.readlines()
act_res_file.close()
active_res = re.split(r',',act_res_lines[0])
print 'Active site residues:', act_res_lines[0].rstrip()


# extract active site residue indices
act_sel_list = []
for res in active_res:
    info = re.split(r':',res)
    if len(info) == 2:
        act_sel_list.append('(chain '+info[1].strip()+' and resnum '+info[0].strip()+')')
act_selector = str(' or '.join(act_sel_list))
no_act_res = len(act_sel_list)
if no_act_res > 0:
    try:
        act_res_indices = calphas.select(act_selector).getResindices()
    except AttributeError:
        print 'No residues found matching your active residues'
        sys.exit()
else:
    print 'No active residues could be detected - please ensure that the format is correct, including chain labels'
    sys.exit()
# check all active residues are included
if len(act_res_indices) != len(active_res):
    print 'One or more active site residues could not be read, or is missing in the PDB file; continuing anyway'


# fpocket information
pocket_info_list = [
    'Score :',
    'Druggability Score :',
    'Number of Alpha Spheres :',
    'Total SASA :',
    'Polar SASA :',
    'Apolar SASA :',
    'Volume :',
    'Mean local hydrophobic density :',
    'Mean alpha sphere radius :',
    'Mean alp. sph. solvent access :',
    'Apolar alpha sphere proportion :',
    'Hydrophobicity score:',
    'Volume score:',
    'Polarity score:',
    'Charge score :',
    'Proportion of polar atoms:',
    'Alpha sphere density :',
    'Cent. of mass - Alpha Sphere max dist:',
    'Flexibility :'
]


print 'Extracting pocket information'
# extract pocket information into a dictionary
pockets_info = extract_info('',input_prefix,pocket_info_list)
# extract pocket residues into a dictionary
pockets_res, pockets_res_name = extract_pockets('',input_prefix,calphas)
# extract distances of each pocket to the active site into a dictionary
pockets_dist = dist_to_active(calphas,act_res_indices,pockets_res)


# RUN NMA
# calculate normal modes without perturbation
print 'Calculating NMA without perturbation'
anm = ANM('normal')
anm.buildHessian(nma_selection,cutoff=nma_cutoff,gamma=1)
anm.calcModes(n_modes=None)
anm_used, sel_active = sliceModel(anm,nma_selection,act_selector)


# iterate over each pocket and perturb the normal modes
print 'Calculating NMA with perturbation'
results = {}
for pocket_no in pockets_res:
    # calculate normal modes with perturbation to residues in pocket
    anm_mod = ANM('modified')
    # do not perturb active residues directly
    indices_to_perturb = [i for i in pockets_res[pocket_no] if i not in act_res_indices]
    anm_mod.buildHessian(nma_selection,cutoff=nma_cutoff,gamma=GammaResidue(nma_selection,res_nos=indices_to_perturb,frac_change=spring_frac))
    anm_mod.calcModes(n_modes=None)
    anm_mod_used, sel_mod_active = sliceModel(anm_mod,nma_selection,act_selector)
    # measure the change in the normal modes
    diff_100 = quantify_difference(anm_used,anm_mod_used,no_modes=100)
    diff_200 = quantify_difference(anm_used,anm_mod_used,no_modes=200)
    diff_all = quantify_difference(anm_used,anm_mod_used,no_modes=len(anm))
    # add results to results dictionary
    results[pocket_no] = [diff_100,diff_200,diff_all]


# form title line
title_string = 'fpocket_rank\tpocket_size\tdist_to_active\tno_pockets\tC_100\tE_100\tC_200\tE_200\tC_all\tE_all'
clean_string = title_string
for info_id in pocket_info_list:
    title_string += '\t'+info_id
    clean_string += '\tfpocket_'+re.sub('[:.\-\s+]','',info_id)
title_string += '\tresidues'
clean_string += '\tresidues'


# form output line for each pocket
lines = []
for pocket_no in range(len(pockets_res)):
    # write fpocket rank
    line = str(pocket_no)
    # write pocket size
    line += '\t'+str(len(pockets_res[pocket_no]))
    # write distance to active site
    line += '\t'+str(pockets_dist[pocket_no])
    # write number of pockets
    line += '\t'+str(len(pockets_res))
    # write combined NM output for 100 modes
    line += '\t'+str(results[pocket_no][0])
    line += '\t'+str(results[pocket_no][0]/len(pockets_res[pocket_no]))
    # write combined NM output for 200 modes
    line += '\t'+str(results[pocket_no][1])
    line += '\t'+str(results[pocket_no][1]/len(pockets_res[pocket_no]))
    # write combined NM output for all modes
    line += '\t'+str(results[pocket_no][2])
    line += '\t'+str(results[pocket_no][2]/len(pockets_res[pocket_no]))
    # write pocket information
    for info_id in pocket_info_list:
        line += '\t'+str(pockets_info[pocket_no][info_id])
    # write residues in pocket
    line += '\t'+','.join(pockets_res_name[pocket_no])
    lines.append(line)


# FORM SVM INPUT FILE
# optimal SVM features and the ranges in the training set
svm_features = [
    ['Number of Alpha Spheres :_r','raw'],
    ['E_200','ranked'],
    ['Score :_r','raw'],
    ['E_all','ranked'],
    ['dist_to_active_r','raw'],
    ['pocket_size_r','raw'],
    ['fpocket_rank_r','raw']
]

svm_features_range = {
    'Number of Alpha Spheres :_r':(36.0, 2276.0),
    'Score :_r':(-6.782, 70.118),
    'dist_to_active_r':(0.1, 110.6),
    'pocket_size_r':(2.0, 275.0),
    'fpocket_rank_r':(0.0, 83.0),
}


# split output lines for SVM sorting
print 'Forming SVM input file'
breaks = []
for line in lines:
    split = re.split(r'\t+',line.rstrip())
    # convert to float for sorting; ignore residue names
    split = [float(i) for i in split[:-1]]
    # convert to int for readability
    split[0] = int(split[0])
    breaks.append(split)
title_break = re.split(r'\t+',title_string.rstrip())


# form dictionary for SVM features
pocket_info = {}
for i in range(len(breaks)):
    pocket_info[str(breaks[i][0])] = {}
for feature in svm_features:
    # find column to sort by
    index = 'not_found'
    for column in title_break:
        if feature[1] == 'ranked' and column == feature[0]:
            index = title_break.index(column)
            break
        elif feature[1] == 'raw' and column+'_r' == feature[0]:
            index = title_break.index(column)
            break
    if index != 'not_found':
        # sort by correct column
        breaks.sort(key=itemgetter(index),reverse=True)
        # cycle through pockets
        for i in range(len(breaks)):
            if feature[1] == 'ranked':
                # append fractional ranking
                pocket_info[str(breaks[i][0])][feature[0]] = round((i + 1) / float(len(breaks)),4)
            elif feature[1] == 'raw':
                # scale value using min and max to be between 0 and 1
                scaled_value = (breaks[i][index] - svm_features_range[feature[0]][0]) / float(svm_features_range[feature[0]][1] - svm_features_range[feature[0]][0])
                pocket_info[str(breaks[i][0])][feature[0]] = round(scaled_value,4)


# write SVM input file
svm_in_file = open(input_prefix+'.svm','w')
for j in range(len(pocket_info)):
    line = '0 '
    for i in range(len(svm_features)):
        if svm_features[i][0] in pocket_info[str(j)]:
            # feature ID is an integer for SVMlight input
            line += ' '+str(i+1)+':'+str(pocket_info[str(j)][svm_features[i][0]])
    svm_in_file.write(line+'\n')
svm_in_file.close()


# RUN SVM
print 'Running SVM'
svm_classify_cmd = svm_light_dir+'/svm_classify '+input_prefix+'.svm '+allopred_dir+'/svm_model.txt '+input_prefix+'.pred'
os.system(svm_classify_cmd)


# FORM OUTPUT FILE
# read and sort prediction file
print 'Forming output file'
svm_pred_file = open(input_prefix+'.pred')
svm_pred_lines = svm_pred_file.readlines()
svm_pred_file.close()
ordering = []
for i in range(len(svm_pred_lines)):
    ordering.append([i,svm_pred_lines[i].rstrip()])
ordering.sort(key=itemgetter(1),reverse=True)
# remove prediction file
os.system('rm '+input_prefix+'.pred')


# write tab-delimited output file
out_file = open(input_prefix+'.out','w')
out_file.write('# - AlloPred output -\n')
out_file.write('#\n')
out_file.write('# PDB file: '+pdb_file+'\n')
out_file.write('# Active site residue file: '+act_res_filepath+'\n')
out_file.write('# Active site residues: '+act_res_lines[0].rstrip()+'\n')
out_file.write('# Pocket directory: '+pocket_folder+'\n')
out_file.write('#\n')
out_file.write('# C_n is the NMA effect over n modes; E_n is the NMA effect per perturbed residue over n modes\n')
out_file.write('#\n')
out_file.write('# AlloPred_rank\t'+clean_string+'\n')
out_file.write('#\n')


# write pocket lines in order of AlloPred ranking
for i in range(len(ordering)):
    line = str(i)+'\t'+lines[ordering[i][0]]
    out_file.write(line+'\n')
out_file.close()
print 'Output file written'


# finish
print 'AlloPred finished successfully'

