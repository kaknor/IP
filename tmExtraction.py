 #!/usr/bin/env python
import os
import re
import time
import fileinput

timestr = time.strftime("%Y-%m-%d")

general_file_data = {}
fasta_genes = {}

topology = []


def topology_parser(file_name):
    with open(os.path.join('inputTopologyFiles', file_name + '.csv'), 'r') as f:
        lines = f.read().splitlines(True)
        for line in lines:
            line = re.sub(r'\t+', '\t', line)
            l = line.split('\t')
            topo = re.sub(r'Topology=', '', l[5])
            length = int(re.sub(r'len=', '', l[1]))
            #0 = id; 1 = len; 2 = ; 3 =; 4 =; 5 = topology
            general_file_data[l[0]] = (length, l[2], l[3], l[4], topo)
            topology.append((l[0].strip(' '), topo))
    create_csv(file_name)

def fasta_parser(file_name):
    #if fasta matching is wrong pb is most likely there
    fasta_pattern = r'Invalid IDs:|PBANKA_[MITAP]*[0-9]*|length=[0-9]*|\n|'
    with open(os.path.join('../stageGenes-wp', file_name), 'r') as f:
        lines = f.read().split('>')
        for line in lines:
            l = line.split('|')
            fasta_genes[l[0].strip(' ')] = (re.sub(fasta_pattern, '', l[-1]))

def get_bounds(s):
    part = s.partition('-')
    return (int(part[0]), int(part[2]))
    
def get_tmd(name, topo):
    # 0 = precedent lbound; 1 = suiv upbound; 2 = current lbound; 3 = current upbound
    tmd = [()]
    tab = re.findall("[io]([0-9]+-[0-9]+)", topo)
    if not tab:
        raise ValueError('ERROR WITH TMD, the topology might be wrong for this gene ' + name)
    for d in tab:
        left, right = get_bounds(d)
        if tmd[-1] is ():
            tmd.append((1, left - 1, left, right))
        else:
            tmd.append((tmd[-1][3] + 1, left - 1, left, right))
    return tmd

def get_prot_seq(fasta_genes, name, bounds):
    res = ''
    if general_file_data[name][0] == bounds[1]:
        res = fasta_genes[name][bounds[0]:]
    else:
        res = fasta_genes[name][bounds[0]:bounds[1] + 1]
    return res

def create_line(f, fasta_genes, i, tmd):
    f.write(i[0] + '\t')
    f.write(i[1].strip('\n') + '\t')
    f.write(str(tmd[2]) + ' - ' + str(tmd[3]) + '\t')
    f.write(get_prot_seq(fasta_genes, i[0], (tmd[2], tmd[3])) + '\t')
    f.write(str(tmd[0]) + ' - ' + str(tmd[1]) + '\t')
    f.write(get_prot_seq(fasta_genes, i[0], (tmd[0], tmd[1])) + '\t')

def create_last_ntmd(f, fasta_genes, i, lbound):
    l = general_file_data[i[0]][0]
    f.write(i[0] + '\t')
    f.write(i[1].strip('\n') + '\t')
    f.write('N/A\t')
    f.write('N/A\t')
    seq = get_prot_seq(fasta_genes, i[0], (lbound, l))
    f.write(str(lbound + 1) + ' - ' + str(l)  + '\t')
    f.write(seq + '\t')
    
def create_csv(file_name):
    with open('outputTopologyFiles/' + file_name + '_' + timestr + '.csv', 'w') as f:
        f.write('[ID]' + '\t' + '[TOPOLOGY]' + '\t' + '[TMD #]' + '\t' + '[TMD protein sequence]' + '\t'
                + '[NON TMD #]' + '\t' + '[NON TMD protein sequence]\n')
        for i in topology:
            tmd = get_tmd(i[0], i[1])[1:]
            for d in tmd:
                create_line(f, fasta_genes, i, d)
                f.write('\n')
            if not tmd:
                print(i[0])
            upbound = tmd[-1][3]
            if general_file_data[i[0]][0] > upbound:
                create_last_ntmd(f, fasta_genes, i, upbound)
                f.write('\n')

fasta_parser('Pf3D7GenesByTransmembraneDomains')
topology_parser('TopologyFpPf3D7')

# fasta_parser('liverStageGenes-wp')
# topology_parser('TopologyPbLs')

# fasta_parser('testFpPf3D7')
# topology_parser('test')
