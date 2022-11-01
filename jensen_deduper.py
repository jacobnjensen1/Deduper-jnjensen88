#!/usr/bin/env python

# fix - strand position: leftmost + [DNMrightS]  ((\d+)[DMN])|((\d+)S$) Look at match[1] and match[3]
# fix + strand position: leftmost - [leftS]  (\d+)S.+ Don't need to index

import re
import argparse

def get_args():
  parser = argparse.ArgumentParser(description="Deduplicates a sorted sam file where all reads are uniquely mapped")
  parser.add_argument("-f", "--file", help="input sam file", required=True)
  parser.add_argument("-o", "--outfile", help="output fasta file", required=True)
  parser.add_argument("-u", "--umi", help="text file containing list of UMIs", required=True)
  return parser.parse_args()

args = get_args()

QNAME_INDEX = 0
BITSTRING_INDEX = 1
RNAME_INDEX = 2
POS_INDEX = 3
CIGAR_INDEX = 5

with open(args.umi) as UMIfile:
  umis = {entry.strip() for entry in UMIfile.readlines()}

def getUMI(record):
  return record[QNAME_INDEX].split(":")[-1]

def hasGoodUMI(record):
  return getUMI(record) in umis

def isPlusStrand(record):
  return int(record[BITSTRING_INDEX]) & 16 != 16

def adjustedPosition(record):
  if isPlusStrand(record):
    adjustSearch = re.findall(r"(\d+)S.+", record[CIGAR_INDEX])
    if len(adjustSearch) == 0:
      toAdjustBy = 0
    else:
      toAdjustBy = int(adjustSearch[0])
    # print(toAdjustBy)
    return int(record[POS_INDEX]) - toAdjustBy
  else: # is on minus strand
    adjustSearch = re.findall(r"((\d+)[DMN])|((\d+)S$)", record[CIGAR_INDEX])
    if len(adjustSearch) == 0:
      toAdjustBy = 0
    else:
      toAdjustBy = 0
      for entry in adjustSearch:
        if entry[1]:
          toAdjustBy += int(entry[1])
        if entry[3]:
          toAdjustBy += int(entry[3])
    # print(toAdjustBy)
    return int(record[POS_INDEX]) + toAdjustBy
    

# print(isPlusStrand(["a", "0"])) #should be true
# print(isPlusStrand(["a", "15"]))    #true
# print(isPlusStrand(["a", "16"]))    #false
# print(isPlusStrand(["a", "17"]))    #false
# print(isPlusStrand(["a", "83"]))    #false
# print(isPlusStrand(["a", "163"]))   #true

# adjustedPosition(["a", "0", "b", "45", "1", "1S40M"])
# adjustedPosition(["a", "0", "b", "45", "1", "10S40M20S"])
# adjustedPosition(["a", "0", "b", "45", "1", "40M"])
# adjustedPosition(["a", "0", "b", "45", "1", "40M1S"])

# adjustedPosition(["a", "16", "b", "45", "1", "1S40M1S"])
# adjustedPosition(["a", "16", "b", "45", "1", "40M10D1S"])
# adjustedPosition(["a", "16", "b", "45", "1", "40M10D"])
# adjustedPosition(["a", "16", "b", "45", "1", "50M1000N40M10D"])
# adjustedPosition(["a", "16", "b", "45", "1", "10S50M1000N40M10D10S"])

# print(hasGoodUMI(["LKJLKJ:AACGCCAT"]))
# print(hasGoodUMI(["LKJLKJ:AACGCNAT"]))

removedDuplicateCount = 0
removedBadUMICount = 0
foundInChrom = set()  #holds tuples of (UMI, pos, strand)
currentChrom = ""
with open(args.file) as inFile, open(args.outfile, "w") as outFile:
  while True:
    line = inFile.readline()
    # print(line)
    if line == "":
      break
    if line.startswith("@"):
      outFile.write(line)
      continue    
    line = line.strip().split()     #LINE IS A LIST AFTER THIS POINT
    if not hasGoodUMI(line):
      removedBadUMICount += 1
      continue
    umi = getUMI(line)
    pos = adjustedPosition(line)
    strand = "+" if isPlusStrand(line) else "-"
    recordTuple = (umi, pos, strand)
    if line[RNAME_INDEX] != currentChrom:
      foundInChrom.clear()
      currentChrom = line[RNAME_INDEX]
    if recordTuple not in foundInChrom:
      foundInChrom.add(recordTuple)
      outFile.write("\t".join(line))
      outFile.write("\n")
    else:
      removedDuplicateCount += 1

print(f"removed due to UMI: {removedBadUMICount}")
print(f"removed due to duplication: {removedDuplicateCount}")