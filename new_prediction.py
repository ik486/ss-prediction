# Protein Secondary structure prediction
#
# Copyright (c) 2013-2014 Ignatius Kunjumon <ignatius.kunjumon@gmail.com>
#
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA  02111-1307  USA
#
# $Id: $
#


import sys
import re

aminoAcidNames = {
        "A": ["Ala","Alanine"],
        "R": ["Arg","Arginine"],
        "D": ["Asp","Aspartic Acid"],
        "N": ["Asn","Asparagine"],
        "C": ["Cys","Cysteine"],
        "E": ["Glu","Glutamic Acid"],
        "Q": ["Gln","Glutamine"],
        "G": ["Gly","Glycine"],
        "H": ["His","Histidine"],
        "I": ["Ile","Isoleucine"],
        "L": ["Leu","Leucine"],
        "K": ["Lys","Lyscine"],
        "M": ["Met","Methionine"],
        "F": ["Phe","Phenylalanine"],
        "P": ["Pro","Proline"],
        "S": ["Ser","Serine"],
        "T": ["Thr","Threonine"],
        "W": ["Trp","Tryptophan"],
        "Y": ["Tyr","Tyrosine"],
        "V": ["Val","Valine"]
    }



class Q3Stat:
	def __init__(self):
		self.read_dssp_q3()
		return

	def read_dssp_q3(self):
		fd = open("DSSP_Q3_Data.txt")
		data = fd.read()
		fd.close()
		data = data.split("\n")
		data.pop()
		self.seqDict = {}
		self.ssDict  = {}
		for d in data:
			pdbid, chainid, seq, ss = d.split(":")
			key = pdbid+":"+chainid
			self.seqDict[key] = seq
			self.ssDict[key]  = ss
		return 

	def getQ3Stat(self, cidList, iN=4):
		out = {}
		for cid in cidList:
			seqStr = self.seqDict[cid]
			ssStr  = self.ssDict[cid]
			n = len(ssStr)
			for i in range(n-iN):
				aa = seqStr[i:i+iN]
				bb =  ssStr[i:i+iN]
				if not out.has_key(aa):  out[aa] ={}
				if not out[aa].has_key(bb):  out[aa][bb] =0
				out[aa][bb] += 1
		output = {}
		for key1 in sorted(out.keys()):
			result = []
			total = 0
			for key2 in out[key1].keys():
				total += out[key1][key2]
			for i in range(iN):
				h_val = 0
				e_val = 0
				c_val = 0
				for key2 in out[key1].keys():
					if key2[i] == 'H': h_val += out[key1][key2]
					elif key2[i] == 'E': e_val += out[key1][key2]
					else:
						c_val += out[key1][key2]
				mydict = {}
				if total == 0:
					mydict['H'] = 0.0
					mydict['E'] = 0.0
					mydict['C'] = 0.0
				else:
					mydict['H'] = float(h_val)/total
					mydict['E'] = float(e_val)/total
					mydict['C'] = float(c_val)/total
				result.append(mydict)
			output[key1] = result
		return output


def calculate_probability(ipos, seq_length, seq, ss_type, stat):
	N = seq_length -1
	if ipos == 0:
		seq1 = seq[:4]
		stat1 = stat.get(seq1, None)
		if stat1 is None: return 0.0
		p1 = stat1[0].get(ss_type, 0.0)
		return p1
	elif ipos == 1:
		seq1 = seq[:4]
		stat1 = stat.get(seq1, None)
		if stat1 is None: return 0.0
		p1 = stat1[1].get(ss_type, 0.0)
		#----------------------------------------
		seq2 = seq[1:5]
		stat2 = stat.get(seq2, None)
		if stat2 is None: return 0.0
		p2 = stat2[0].get(ss_type, 0.0)
		#----------------------------------------
		return p1*p2
	elif ipos == 2:
		seq1 = seq[:4]
		stat1 = stat.get(seq1, None)
		if stat1 is None: return 0.0
		p1 = stat1[2].get(ss_type, 0.0)
		#----------------------------------------
		seq2 = seq[1:5]
		stat2 = stat.get(seq2, None)
		if stat2 is None: return 0.0
		p2 = stat2[1].get(ss_type, 0.0)
		#----------------------------------------
		seq3 = seq[2:6]
		stat3 = stat.get(seq3, None)
		if stat3 is None: return 0.0
		p3 = stat3[0].get(ss_type, 0.0)
		#----------------------------------------
		return p1*p2*p3
	elif ipos == N:
		seq4 = seq[-4:]
		stat4 = stat.get(seq4, None)
		if stat4 is None: return 0.0
		p4 = stat4[3].get(ss_type, 0.0)
		return p4
	elif ipos == N-1:
		seq4 = seq[-4:]
		stat4 = stat.get(seq4, None)
		if stat4 is None: return 0.0
		p4 = stat4[2].get(ss_type, 0.0)
		#----------------------------------------
		seq3 = seq[-5:-1]
		stat3 = stat.get(seq3, None)
		if stat3 is None: return 0.0
		p3 = stat3[3].get(ss_type, 0.0)
		return p4*p3
	elif ipos == N-2:
		seq4 = seq[-4:]
		stat4 = stat.get(seq4, None)
		if stat4 is None: return 0.0
		p4 = stat4[1].get(ss_type, 0.0)
		#----------------------------------------
		seq3 = seq[-5:-1]
		stat3 = stat.get(seq3, None)
		if stat3 is None: return 0.0
		p3 = stat3[2].get(ss_type, 0.0)
		#----------------------------------------
		seq2 = seq[-6:-2]
		stat2 = stat.get(seq2, None)
		if stat2 is None: return 0.0
		p2 = stat2[3].get(ss_type, 0.0)
		return p4*p3*p2
	else:
		seq1 = seq[ipos:ipos+4]
		stat1 = stat.get(seq1, None)
		if stat1 is None: return 0.0
		p1 = stat1[0].get(ss_type, 0.0)
		#----------------------------------------
		seq2 = seq[ipos-1:ipos+3]
		stat2 = stat.get(seq2, None)
		if stat2 is None: return 0.0
		p2 = stat2[1].get(ss_type, 0.0)
		#----------------------------------------
		seq3 = seq[ipos-2:ipos+2]
		stat3 = stat.get(seq3, None)
		if stat3 is None: return 0.0
		p3 = stat3[2].get(ss_type, 0.0)
		#----------------------------------------
		seq4 = seq[ipos-3:ipos+1]
		stat4 = stat.get(seq4, None)
		if stat4 is None: return 0.0
		p4 = stat4[3].get(ss_type, 0.0)
		#----------------------------------------
		return p1*p2*p3*p4

def calculate_one(q3stat, pdbStat, seq):
	ss_list = ['H', 'E', 'C']
	stat = pdbStat
	seq_length = len(seq)
	out = []
	pre = []
	for i in range(seq_length):
		ipos = i
		mylist = []
		mydict = {}
		for ss_type in ss_list:
			p = calculate_probability(ipos, seq_length, seq, ss_type, stat)
			mylist.append([p, ss_type])
			mydict[ss_type] = p
		mylist.sort()
		pre.append(mylist[-1][1])
		out.append(mydict)
	result = "".join(pre)
	print result


def sampling(level = 2):
	fd = open("new_q3.data")
	data = fd.read()
	fd.close()
	data = data.split("\n")
	data.pop()
	out = []
	for d in data:
		vals = d.split(",")[:6]
		vals = ",".join(vals)
		vals = eval(vals)
		pdbid = vals[0]
		per = vals[4]
		pdbid = pdbid[:4]+":"+pdbid[4]
		out.append([per, pdbid])
	out.sort()
	N = len(out)
	samples = []
	for i in range(level):
		j1 = N*i/level
		j2 = N*(i+1)/level
		if j2 == N: j2 += 1
		samples.append(out[j1:j2])
	output = []
	for sample in samples:
		pdbidList = []
		for per, pdbid in sample:
			pdbidList.append(pdbid)
		output.append(pdbidList)
	return output


def predict(N, seq):
	sampleList = sampling(N)
	for sample in sampleList:
		q3stat = Q3Stat()
		pdbidList = sample
		pdbStat = q3stat.getQ3Stat(pdbidList)
		calculate_one(q3stat, pdbStat, seq)

if __name__ == '__main__':
        if len(sys.argv) < 3:
		print
                print "python", sys.argv[0], "<NO_OF_SAMPLING> SEQ"
		print "\tNO_OF_SAMPLING = 1, 2, 3, 4, ... 512 ..."
		print "\tExample"
		print "\t\t"+sys.argv[0],"4 RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGAIGPWENL"
                sys.exit(0)
	no = int(sys.argv[1])
	seq = sys.argv[2]

	amino_letters = aminoAcidNames.keys()
	amino_letters.sort()
	out = []
	for s in seq: 
		if s not in amino_letters: 
			out.append(s)
	if len(out):
		print "following letters are not valid amino acids:", ",".join(out)
		sys.exit(0)
	predict(no, seq)

