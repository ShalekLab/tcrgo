"""Relative imports"""
from .bamdict import BAMDict 
from .umi import UMI 
from .barcode import Barcode 
from .read import Read
from .reference import ReferenceDict, Reference
import pysam

IndexedReads = pysam.libcalignmentfile.IndexedReads
AlignedSegment = pysam.libcalignedsegment.AlignedSegment
BAM = pysam.libcalignmentfile.AlignmentFile