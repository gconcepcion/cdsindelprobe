"""
A simple tool to probe for single bp indels in bamfiles.
"""
import os
import sys
import logging
import argparse

from cdsindelprobe import utils as u

log = logging.getLogger()

__author__ = 'gconcepcion'
__version__ = 0.1


def get_parser():
    """

    :return:
    """
    desc = "Tool to probe for single bp frameshifts in coding region alignments"
    parser = argparse.ArgumentParser(version=__version__, description=desc)
    parser.add_argument('cds_bam', type=str,
                        help="path to bamfile containing CDS aligned to a reference")
    parser.add_argument('refpath', type=str,
                        help="path to reference fasta to parse")
    parser.add_argument('--raw_bam', dest='raw_bam',
                        help='Optional check indel positions for homopolymers. '
                             'Must supply path to bam file containing raw reads mapped to '
                             'reference')
    parser.add_argument('--debug', action='store_true')
    return parser


def run(log, cds_bam, reffile, raw_bam):
    """

    :param log:
    :param cds_bam:
    :param reffile:
    :param raw_bam:
    :return:
    """
    unique_indels = u.parse_bamfile(log, cds_bam)

    if raw_bam is not None:
        u.check_homopolymers(log, reffile, raw_bam, unique_indels)

    return


def main():
    """

    :return:
    """
    parser = get_parser()
    args = parser.parse_args()
    bamfile = args.cds_bam
    reffile = args.refpath
    rawreads = args.raw_bam

    if args.debug:
        u.setup_log(log, level=logging.INFO)
    else:
        u.setup_log(log, level=logging.INFO)

    if os.path.exists(bamfile):
        run(log, bamfile, reffile, rawreads)
    else:
        log.info("Can't find sam file to parse! {x}".format(x=bamfile))

    return


if __name__ == '__main__':
    sys.exit(main())
