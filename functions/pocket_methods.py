import re
import numpy as np
from prody import *


def extract_pockets(filepath,pdbid,selection):
    """Returns dictionaries where the key is pocket ID (starting at zero) and the value is a list of residue indices or numbers/chains in the pocket."""
    # dictionary where key is pocket ID (starting at zero) and value is list of residue indices in pocket
    pockets_res = {}
    # dictionary where key is pocket ID (starting at zero) and value is list of residue numbers/chains in pocket
    pockets_res_name = {}
    pocket_count = count_pockets(filepath,pdbid)
    for pocket_no in range(pocket_count):
        pockets_res[pocket_no] = []
        pockets_res_name[pocket_no] = []
        in_file = open(filepath+pdbid+'_out/pockets/pocket'+str(pocket_no)+'_atm.pdb')
        for line in in_file:
            if line[:4] == 'ATOM':
                # add residue name
                name = line[22:26].lstrip()+':'+line[21]
                if name not in pockets_res_name[pocket_no]:
                    pockets_res_name[pocket_no].append(name)
                # add residue ID
                selector = 'chain '+line[21]+' and resnum '+line[22:26].lstrip()
                try:
                    index = selection.select(selector).getResindices()[0]
                    if index not in pockets_res[pocket_no]:
                        pockets_res[pocket_no].append(index)
                # selector might not have a resindex
                except AttributeError:
                    pass
                # cannot select negative resnums
                except atomic.select.SelectionError:
                    pass
        in_file.close()
    return pockets_res, pockets_res_name


def extract_info(filepath,pdbid,info_id_list):
    """Returns a dictionary where the key is pocket ID (starting at zero) and the value is a dictionary of information points."""
    pockets_info = {}
    pocket_file = open(filepath+pdbid+'_out/'+pdbid+'_info.txt')
    pocket_lines = pocket_file.readlines()
    pocket_file.close()
    # create inner dictionaries
    counter = 0
    for line in pocket_lines:
        if line[:6] == 'Pocket':
            pockets_info[counter] = {}
            counter += 1
    # populate inner dictionaries
    for info_id in info_id_list:
        counter = 0
        for line in pocket_lines:
            if line.lstrip()[:len(info_id)] == info_id:
                split = re.split(r'\s+',line.rstrip())
                pockets_info[counter][info_id] = float(split[-1])
                counter += 1
    return pockets_info


def count_pockets(filepath,pdbid):
    """Counts the number of pockets found by fpocket."""
    pocket_count = 0
    pocket_file = open(filepath+pdbid+'_out/'+pdbid+'_info.txt')
    for line in pocket_file:
        if line[:6] == 'Pocket':
            pocket_count += 1
    pocket_file.close()
    return pocket_count


def dist_to_active(calphas,active_indices,pockets_res):
    """Returns a dictionary where the key is pocket index and the value is distance to the active site in Angstrom."""
    active_selector = 'resindex '+' or resindex '.join([str(i) for i in active_indices])
    active_coords = calphas.select(active_selector).getCoords()
    active_site = sum(active_coords)/float(len(active_indices))
    # dictionary where key is pocket index and value is distance to active site in A
    pockets_dist = {}
    for pocket in pockets_res:
        pocket_selector = 'resindex '+' or resindex '.join([str(i) for i in pockets_res[pocket]])
        pocket_coords = calphas.select(pocket_selector).getCoords()
        pocket_center = sum(pocket_coords)/float(len(pockets_res[pocket]))
        diff = pocket_center - active_site
        dist = np.sqrt(diff.dot(diff))
        pockets_dist[pocket] = round(dist,1)
    return pockets_dist

