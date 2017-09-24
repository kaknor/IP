#!/usr/bin/env python
import os
import re
import time

timestr = time.strftime("%Y-%m-%d")

general_file_data = {}
topology = []

fasta_genes = []

def topology_parser(file_name):
    with open(os.path.join('inputTopologyFiles', file_name + '.csv'), 'r') as f:
        lines = f.read().splitlines(True)
        for line in lines:
            l = line.split('\t')
            topo = re.sub(r'Topology=', '', l[5])
            length = int(re.sub(r'len=', '', l[1]))
            #0 = id; 1 = len; 2 = ; 3 =; 4 =; 5 = topology
            general_file_data[l[0]] = (length, l[2], l[3], l[4], topo)
            topology.append((l[0].strip(' '), topo))
    create_csv(file_name)

def fasta_parser(file_name):
    with open(os.path.join('../stageGenes-wp', file_name), 'r') as f:
        lines = f.read().split('>')
        for line in lines:
            l = line.split('|')
            fasta_genes.append((l[0].strip(' '), re.sub('\n', '', re.sub(r'length=[0-9]*\n', '', l[-1]))))

def get_bounds(s):
    part = s.partition('-')
    return (int(part[0]), int(part[2]))
    

# def get_topo(name, topo):
#     tmd = []
#     ntmd = [()]
#     for i in range(0, len(topo)):
#         if topo[i] == 'o':
#             left, right = get_bounds(topo[i + 1:], 'i')
#             if (left != '' and right != ''):
#                 ntmd.append((int(left), int(right)))
#             else:
#                 tmd.append(((tmd[-1][1][1], general_file_data[name][0]), ()))
#         if topo[i] == 'i':
#             left, right = get_bounds(topo[i + 1:], 'o')
#             if (left != '' and right != ''):
#                 tmd.append((ntmd[-1], (int(left), int(right))))
#             else:
#                 tmd.append(((ntmd[-1]), (ntmd[-1][1], general_file_data[name][0])))
#     return tmd, ntmd

def get_topo(name, topo):
    tmd = [()]
    for d in re.findall("[io]([0-9]*-[0-9]*)", topo):
        left, right = get_bounds(d)
        if tmd[-1] is ():
            tmd.append((1, left - 1, left, right))
        else:
            tmd.append((tmd[-1][3] + 1, right - 1, left, right))
    return tmd

def get_prot_seq(d_fasta, name, bounds):
    return d_fasta[name][bounds[0]:bounds[1]]

def form_line(f, d_fasta, i, tmdom, ntmdom):
    f.write(i[0] + '\t')
    f.write(i[1].strip('\n'))
    if tmdom:
        f.write('\t' + str(tmdom[0]) + ' - ' + str(tmdom[1]) + '\t' + get_prot_seq(d_fasta, i[0], tmdom))
    if ntmdom:
        f.write('\t' + str(ntmdom[0]) + ' - ' + str(ntmdom[1]) + '\t' + get_prot_seq(d_fasta, i[0], ntmdom))
    f.write('\n')
            
def create_csv(file_name):
    d_fasta = dict(fasta_genes)
    with open('outputTopologyFiles/' + file_name + '_' + timestr + '.csv', 'w') as f:
        f.write('[ID]' + '\t' + '[TOPOLOGY]' + '\t' + '[TMD #]' + '\t' + '[TMD protein sequence]' + '\t'
                + '[NON TMD #]' + '\t' + '[NON TMD protein sequence]\n')
        for i in topology:
            tmd = get_topo(i[0], i[1])
            print(tmd)
            # for tmdom, ntmdom in zip(tmd, ntmd):
            #     form_line(f, d_fasta, i, tmdom, ntmdom)
            # if len(tmd) == len(ntmd):
            #     continue
            # if len(tmd) > len(ntmd):
            #     form_line(f, d_fasta, i, tmd[-1], [])
            # else:
            #     form_line(f, d_fasta, i, [], ntmd[-1])
            
def print_topology():
    for i in topology:
        print(i)

def print_fasta():
    for i in fasta_genes:
        print(i)

#fasta_parser('liverStageGenes-wp')
#topology_parser('TopologyPbLs')
fasta_parser('testls')
topology_parser('test')
