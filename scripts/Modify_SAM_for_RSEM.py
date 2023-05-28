import sys
import os.path
import re

def filter_alignments(buf,total_reads,out_fh):
	# Criterion 1: Minimal mismatch(es)/edit distance
	min_mismatch = 3
	for x in buf:
		a = x.split()
		nm = int(re.findall("[0-9]+",a[11])[0])
		if min_mismatch > nm:
			min_mismatch = nm

	# Criterion 2: Maximal alignments length
	max_match = 0
	for x in buf:
		a = x.split()
		nm = int(re.findall("[0-9]+",a[11])[0])
		if min_mismatch == nm:
			md = int(a[5][:-1])
			if md > max_match:
				max_match = md
	buf_tmp = []
	for x in buf:
		a = x.split()
		nm = int(re.findall("[0-9]+",a[11])[0])
		md = int(a[5][:-1])
		if nm==min_mismatch and md==max_match:
			buf_tmp.append("\t".join(a))

	if len(buf_tmp)>0: total_reads[0] += 1
	for x in buf_tmp:
		out_fh.write(x + "\n")
	buf_tmp[:] = []
	buf[:] = []
	

# Variable arguments
argvs = sys.argv
sam_file = argvs[1]
out_file = argvs[2]
buf = []
rid = None
total_reads = [0]
out_fh = open(out_file,"w")
with open(sam_file) as sam_fh:
	for line in sam_fh:
		line = line.rstrip()
		if line[0]=="@":
			out_fh.write(line + "\n")
		else:
			a = line.split()
			if int(a[1]) & 0x4:
				continue
			b = re.findall("(\d*M)",a[5])
			s = []
			if len(b)==1:
				s = a[5].split(b[0])
				if len(s)==2:
					nm_tag = ""
					md_tag = ""
					for tag in a[11:]:
						if tag[0:2] == "NM": nm_tag = tag
						elif tag[0:2] == "MD": md_tag = tag 
					nm = re.findall("[0-9]+",nm_tag)
					if int(nm[0]) < 3: # Mismatch <= 2
						a[5] = b[0]
						match_start = 0
						if len(s[0])>0:
							match_start += int(s[0][:-1])
						match_end = match_start + int(b[0][:-1]) 
						a[9] = a[9][match_start:match_end]
						a[10] = a[10][match_start:match_end]
						a[11] = nm_tag
						a[12] = md_tag
						if(a[0] != rid):
							filter_alignments(buf,total_reads,out_fh)
						buf.append("\t".join(a[0:13]))
						rid = a[0]
filter_alignments(buf,total_reads,out_fh)
out_fh.close()
print "Total mapped reads: " + str(total_reads[0])






