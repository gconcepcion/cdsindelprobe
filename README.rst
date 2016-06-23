cdsindelprobe
=============

Overview
--------

This is a basic python utility to probe Coding Sequence Alignments for single basepair
indel events against a reference genome. The main output is an `indels.bed` track that you can
load into IGV along side the alignment to quickly identify CDS regions potentially harboring
single basepair indel frameshifts.

Installation
------------

After cloning the repository, installation is simple

.. code-block:: bash

   $ virtualenv ~/venvs/cdsindelprobe
   $ source ~/venvs/cdsindelprobe/bin/activate
   (cdsindelprobe) $ cd /path/to/cdsindelprobe
   (cdsindelprobe) $ pip install -r REQUIREMENTS.txt
   (cdsindelprobe) $ python setup.py install

Running
-------

There are only two required inputs.

.. code-block::

   CDS_BAM = a bamfile containing coding sequences mapped to a reference
   REF_FASTA = path to reference.fasta used in `CDS_BAM` alignment

.. code-block:: bash

   (cdsindel_probe) $ probe-indels.py $CDS_BAM $REF_FASTA

If you have an alignment of raw reads to the consensus, you can optionally pass the parameter `--raw_bam` with a path to
a bamfile containing raw reads aligned to the reference. This will result in generation of homopolymer
distribution plots for each single basepair indel identified.

.. code-block:: bash

   (cdsindel_probe) $ probe-indels.py $CDS_BAM $REF_FASTA --raw_bam /path/to/raw.bamerence.fasta --raw_bam /path/to/raw.bam





this is a work in progress...