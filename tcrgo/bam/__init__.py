"""Relative imports"""
from .bamdict import BAMDict 
from .umi import UMI 
from .barcode import Barcode 
from .read import Read 
import pysam

IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile