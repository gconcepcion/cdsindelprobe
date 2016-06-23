"""
simple utils to mine single bp indels out of alignments
"""
import os
import sys
import logging
from collections import namedtuple
from itertools import groupby

from pandas import DataFrame as df
from matplotlib import pyplot as plt
import pysam

from pbcore.io import FastaTable


def setup_log(alog, level=logging.INFO,
              str_formatter=
              '[%(levelname)s] %(asctime)-15s [%(name)s %(funcName)s %(lineno)d] %(message)s'):
    """setup log handler

    :param alog: a log instance
    :param level: (int) Level of logging debug
    :param file_name: (str, None) if None, stdout is used, str write to file
    :param log_filter: (LogFilter, None)
    :param str_formatter: (str) log formatting str
    """
    alog.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    alog.addHandler(handler)


def make_bed(indel_list):
    """

    :param indel_list:
    :return:
    """
    outfile = os.path.join(os.getcwd(), 'indels.bed')

    stringified = []

    for refname, refstart, refend, itype in indel_list:
        stringified.append((refname, str(refstart), str(refend), itype))

    with open(outfile, 'w') as file_out:
        for line in stringified:
            file_out.write("\t".join(line) + "\n")
    file_out.close()


def parse_reference(refpath):
    """

    :param refpath:
    :return:
    """

    return FastaTable(refpath)


def parse_bamfile(log, bamfile):
    """

    :param log:
    :param bamfile:
    :return:
    """
    insertions = 0

    deletions = 0
    sheets = []

    alignments = pysam.AlignmentFile(bamfile)
    log.info('{x} coding sequences were aligned against the reference genome'.format(
        x=alignments.mapped))

    for iterate, alignment in enumerate(alignments):
        log.debug("Investigating CDS: {x} alignment against {y}".format(
            x=alignment.query_name, y=alignment.reference_name))
        cigar_sheet, icount, dcount = get_cigar_sheet(alignment)

        log.debug("Found {x} insertions; {y} deletions;".format(
            x=icount, y=dcount))

        insertions += icount
        deletions += dcount

        if cigar_sheet is not None:
            log.debug('{x} contains {i} insertion(s) and {d} deletion(s) relative to {z}'
                      .format(x=alignment.query_name,
                              i=icount,
                              d=dcount,
                              z=alignment.reference_name))
        sheets.extend(cigar_sheet)

    unique_indels = list(set(sheets))

    make_bed(unique_indels)
    log.info("{a} CDS alignments scanned".format(a=alignments.mapped))
    log.info("{i} 1 bp Insertions detected".format(i=insertions))
    log.info("{d} 1 bp Deletions detected".format(d=deletions))
    alignments.close()

    return unique_indels


def get_cigar_sheet(alignment):
    """ sift through alignment cigar string for single bp indels

    :param pysam.AlignedSegment:
    :return: cigar_sheets, # of single bp insertions, # of single bp deletions
    """
    aln = alignment
    sheets = []
    pairs = aln.aligned_pairs
    icount = 0
    dcount = 0
    aln_pos = 0

    for operation, length in aln.cigar:
        if operation == 0:
            aln_pos += length
        elif operation == 1:
            # we're only looking at single bp frameshifts
            if length == 1:
                icount += 1
                aln_pos += length
                genomic_position = int(
                    [item[1] for item in pairs if item[0] == aln_pos][0])
                sheet = (aln.reference_name, genomic_position,
                         genomic_position, 'insertion')
                sheets.append(sheet)
            else:
                aln_pos += length

        elif operation == 2:  # skip 'l' number of coordinates in reference
            if length == 1:
                dcount += 1
                genomic_position = int(
                    [item[1] for item in pairs if item[0] == aln_pos][0])
                sheet = (aln.reference_name, genomic_position,
                         genomic_position, 'deletion')
                sheets.append(sheet)
        elif operation == 3:
            pass
        elif operation == 4:
            pass
    return sheets, icount, dcount


def check_homopolymers(log, reffile, bamfile, indels):
    """

    :param log:
    :param reffile: path to fasta
    :param bamfile: path to bamfile with raw reads mapped to reference
    :param indels: list of indel positions to check
    :return:
    """
    ref = parse_reference(reffile)
    alignments = pysam.AlignmentFile(bamfile)
    hp_flank = 100
    summary = []
    log.info("Checking Indels for homopolymer regions")

    for refname, rstart, rend, itype in indels:
        refseq = ref[refname]
        if itype == 'insertion':
            indel_base = refseq.sequence[int(rstart) - 2]
        elif itype == 'deletion':
            # for deletion I think I need to subtract one
            indel_base = refseq.sequence[int(rstart) - 1]

        hp_summary = []
        reads = alignments.fetch(refname, int(rstart), int(rend) + 1)
        for read in reads:
            pairs = read.aligned_pairs
            for read_pos, ref_pos in pairs:
                if ref_pos == int(rstart):
                    if read_pos is not None:
                        position = read_pos
                    else:
                        position = [read_pos for read_pos, ref_pos
                                    in pairs if ref_pos == int(rstart) - 1][0]
            if position != None:
                seq = read.query_sequence[
                    position - hp_flank: position + hp_flank]
                if seq != '':
                    query_hp_len = get_hp_len(seq, hp_flank, indel_base)
                    hp_summary.append((indel_base, query_hp_len))

        plot_hps(hp_summary, refname, rstart)
        log.debug("Plot written {x}_{z}".format(x=refname, z=rstart))
        indel_summary = refname, rstart, indel_base, itype
        summary.append(indel_summary)

    sorted_summary = sorted(summary, key=lambda x: (x[0], x[1]))

    for each in sorted_summary:
        print each
    return


def get_hp_len(seq, hp_flank, base):
    """

    :param seq: subset of read with indel in center
    :param hp_flank: # bp's flanking indel on either side
    :param base: indel base
    :return:
    """
    hp_tuple = namedtuple("homopolymer", "bp len pos")
    pos = 0
    tup_list = []

    for basepair, group in groupby(str(seq)):
        hp_length = len(list(group))
        query_tuple = hp_tuple(basepair, hp_length, pos)
        tup_list.append(query_tuple)
        pos += hp_length


    positions = [hp.pos for hp in tup_list]
    seq_pos = 0
    last_i = 0
    for tup in tup_list:
        if seq_pos < hp_flank:
            seq_pos += tup.len
        elif seq_pos > hp_flank:
            seq_pos -= last_i
            break
        elif seq_pos == hp_flank:
            break
        last_i = tup.len

    index = 0
    for item, tup in enumerate(tup_list):
        if tup.pos == seq_pos:
            index = item

    tup_item = tup_list[index - 1]
    tup_length = tup_item.len
    tup_bp = tup_item.bp
    tup_pos = tup_item.pos

    if tup_bp != base:
        idx = positions.index(tup_pos)
        new_idx = idx + 1
        try:
            tup_item = tup_list[new_idx]
            tup_length = tup_item.len
        except IndexError as error:
            print error
            tup_length = 0

    return tup_length


def plot_hps(hp_summary, name, rstart):
    """ Make homopolymer distribution plots

    :param hp_summary:
    :param name:
    :param rstart:
    :return:
    """
    plotname = "{x}_{p}".format(x=name, p=rstart)
    dataframe = df.from_records(hp_summary)

    try:

        plt.figure()
        dataframe.hist()
        plt.title(plotname)
        plt.savefig('{n}.png'.format(n=plotname))

    except Exception as error:
        print 'Exception: {e}'.format(e=error)
