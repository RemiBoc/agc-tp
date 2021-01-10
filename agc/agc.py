#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
from __future__ import division

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    @brief : Crée un générateur de séquences de longueur l >= minseqlen.
    @param amplicon_file : string, lien vers le fichier fasta
    @param minseqlen : int, longueur minimale des séquences
    @returns : generator, générateur de séquences.
    """

    isfile(amplicon_file)
    with gzip.open(amplicon_file) as fast_a:
        seq = ''
        for line in fast_a:
            if (line.startswith(">") and len(seq) >= minseqlen):
                yield seq
                seq = ''
            elif not line.startswith(">"):
                seq += line.replace(" ", "").replace("\n", "")
            else:
                seq = ''
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    @brief : Fait appel au générateur fourni par read_fasta pour créer un générateur des séquences
    uniques ayant une occurrence O>=mincount ainsi que leur occurrence. Retourne les séquences
    par ordre décroissant d’occurrence.
    @param amplicon_file : string, lien vers le fichier fasta
    @param minseqlen : int, longueur minimale des séquences
    @param mincount : int, comptage minimal des séquences.
    @returns : generator, générateur de séquences uniques et d'occurences.
    """

    fasta_file = read_fasta(amplicon_file, minseqlen)
    sequences = {}
    for seq in fasta_file:
        if seq in sequences:
            sequences[seq] += 1
        else:
            sequences[seq] = 1

    sequences = sorted(sequences.items(), key=lambda item: item[1], reverse=True)

    for seq, count in sequences:
        if count >= mincount:
            yield [seq, count]


def get_chunks(sequence, chunk_size):
    """
    @brief : Prend une séquence et une longueur de segment chunk_size et crée une liste d'au moins
    4 sous-séquences de taille chunk_size non chevauchantes.
    @param sequence : string, séquence
    @param chunk_size : int, taille des segments
    @returns : list, liste de sous séquences.
    """

    chunks_list = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)
                   if len(sequence[i:i+chunk_size]) == chunk_size]
    return chunks_list if len(chunks_list) >= 4 else None

def get_unique(ids):
    """
    Renvoie les ids uniques.
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """
    Renvoie les élements communs.
    """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    @brief : Prend une séquence et une longueur de k-mer et crée un générateur
    de tous les mots de longueur k présents dans cette séquence
    @param sequence : string, sequence
    @param kmer_size : int, longueur des kmer
    @returns, generator, générateur de séquences de taille kmer_size
    """

    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i:(i + kmer_size)]

def get_identity(alignment_list):
    """
    @brief : Calcule le pourcentage d’identité entre deux séquences
    @param alignment_list : list, alignement sous forme de liste
    @returns : float, pourcentage d'itentité entre les deux séquences
    """

    equ_char_count = 0

    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            equ_char_count += 1

    return equ_char_count / len(alignment_list[0]) * 100

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    @brief : Récupère les séquences non chimériques du fichier donné.
    @param amplicon_file : string, lien du fichier d'entrée.
    @param minseqlen : int, longueur minimale des séquences.
    @param mincount : int, comptage minimal des séquences.
    @param chunk_size : int, taille des segments.
    @kmer_size : int, taille des kmers.
    @returns : generator, générateur de séquences non chimériques.
    """

    kmer_dict = {}
    non_chimera = []
    perc_id_matrix = []
    id_seq = 0

    for seq, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        chunks = get_chunks(seq, chunk_size)[:4]

        mates = [search_mates(kmer_dict, sub_seq, kmer_size) for sub_seq in chunks]
        parents = []

        for mate, _ in enumerate(mates):
            parents = common(parents, mates[mate])

        if len(parents) >= 2:

            perc_id_matrix = [[] for _ in range(len(chunks))]
            for parent in parents:
                parent_chunks = get_chunks(non_chimera[parent], chunk_size)
                for index, chunk in chunks:
                    alignment = nw.global_align(chunk, parent_chunks[index])
                    identity = get_identity(alignment)
                    perc_id_matrix[index].append(identity)

        if not detect_chimera(perc_id_matrix):
            kmer_dict = get_unique_kmer(kmer_dict, seq, id_seq, kmer_size)
            non_chimera.append(seq)
            id_seq += 1
            yield [seq, count]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    @brief : Regroupement glouton de sequences.
    @param amplicon_file : string, lien du fichier d'entrée.
    @param minseqlen : int, longueur minimale des séquences.
    @param mincount : int, comptage minimal des séquences.
    @param chunk_size : int, taille des segments.
    @kmer_size : int, taille des kmers.
    @returns : list, liste d'OTU.
    """

    otu_list = []
    chimeras = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    for i in range(len(enumerate(chimeras))):
        otu = True
        for j in range(i + 1, len(chimeras)):
            if get_identity(nw.global_align(chimeras[i][0], chimeras[j][0])) > 97 and \
                chimeras[i][1] >= chimeras[j][1]:
                otu_list.append(chimeras[j])
                otu = False
                break
        if otu:
            otu_list.append(chimeras)
    return otu_list[0]


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(otu_list, output_file):
    """
    @brief : affiche les OTU au format >OTU_{numéro partant de 1}
    occurrence:{nombre d’occurrence à la déréplication}{séquence au format fasta}.
    @param OTU_list : list, liste d'OTU.
    @param output_file : string, lien vers le fichier de sortie.
    """

    with open(output_file, "w") as fasta_file:
        for otu in range(len(otu_list)):
            fasta_file.write(">OTU_{} occurrence:{}\n{}\n"
                             .format(otu, otu_list[1], fill(otu_list[0][0])))


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    @brief : Ajoute les kmers de la séquence d'entrée au dictionnaire donné.
    @param kmer_dict : dict, dictionnaire à remplir.
    @param sequence : string, séquence à découper.
    @param id_seq : int, ID de la séquence.
    @param kmer_size : int, taille des kmers à découper.
    @returns : dict, kmer_dict mis à jour.
    """

    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    @brief : Renvoie les 8 séquences contenues dans kmer_dict les plus proches de
    la séquence donnée.
    @param kmer_dict : dict, dictionnaire de kmers.
    @param sequence : string, sequence à analyser.
    @param kmer_size : int, taille des kmers.
    @returns : string, séquences les plus proches de la séquence donnée.
    """

    counter = Counter([id for kmer in cut_kmer(sequence, kmer_size)
                       if kmer in kmer_dict for _ in kmer_dict[kmer]])
    mates = [mate[0] for mate in counter.most_common(8)]
    return mates


def detect_chimera(perc_identity_matrix):
    """
    @brief : Détecte si une séquence est une chimère.
    @param perc_identity_matrix :
    """

    seq_1 = []
    seq_2 = []
    std_perc = 0
    for perc_id in perc_identity_matrix:

        std_perc += statistics.stdev(id)

        if perc_id[0] not in seq_1:
            seq_1.append(perc_id[0])
        if perc_id[1] not in seq_2:
            seq_2.append(perc_id[1])

    if (len(seq_1) >= 2 or len(seq_2) >= 2) and std_perc/len(perc_identity_matrix) > 5.0:
        return True
    return False


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    otu_list = abundance_greedy_clustering(args.amplicon_file,
                                           args.minseqlen,
                                           args.mincount,
                                           args.chunk_size,
                                           args.kmer_size)
    write_OTU(otu_list, args.output_file)


if __name__ == '__main__':
    main()
