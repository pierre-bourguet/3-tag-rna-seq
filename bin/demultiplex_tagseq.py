from glob import glob
import argparse
import subprocess as sp
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--input", dest = "d_input",
		type=str, required=True,
		                    help="Input dir with the R1 and R2 files")
parser.add_argument("--output", dest = "d_output",
		type=str, required=True,
		                    help="Output dir")
parser.add_argument("--barcodes", dest = "barcodes",
		type=str, required=True,
		                    help="Barcode table: WellPosition<TAB>Name<TAB>Sequence, with a header line")
args = parser.parse_args()

base_dir = args.d_input
output_dir = args.d_output

sp.call(["mkdir","-p", output_dir])

# Identify R1 & R2 files
R1 = glob(base_dir + "/*R1*.fastq*") 
if len(R1) != 1: 
	raise Exception("Identified R1 fastq files number {} != 1".format(len(R1)))   
R1 = R1[0]

R2 = glob(base_dir + "/*R2*.fastq*") 
if len(R2) != 1: 
	raise Exception("Identified R2 fastq files number {} != 1".format(len(R2)))
R2 = R2[0]

# Loading barcodes
barcodes = [x.split("\t") for x in open(args.barcodes,"r").read().split("\n")[1:] if x.strip()]
out_fn = {}

for b in barcodes:
	fn = "{}/{}_{}.fastq".format(output_dir, b[0], b[2])
	out_fn[b[2]] = open(fn, "w")

no_barcode_match_R1 = open("{}/no_barcode_match_R1.fastq".format(output_dir), "w")
no_barcode_match_R2 = open("{}/no_barcode_match_R2.fastq".format(output_dir), "w")

print("Start reading the reads file")
if R1[-2:] == "gz":
	read_R1 = gzip.open(R1, "rt")
else:
	read_R1 = open(R1, "r")

if R2[-2:] == "gz":
	read_R2 = gzip.open(R2, "rt")
else:
	read_R2 = open(R2, "r")

cnt_barcodes = {}
cnt_umis = {}

index = 0
lines_R1 = [read_R1.readline() for i in range(4)]
lines_R2 = [read_R2.readline() for i in range(4)]
while ((len(lines_R1[3]) > 0) and (len(lines_R2[3]) > 0)):
	index = index + 1
	if (index % 10000000) == 0:
		print("Line {}\n".format(index))
	barcode = lines_R2[1][:7]
	cnt_barcodes[barcode] = cnt_barcodes.get(barcode, 0) + 1
	if barcode in out_fn:
		umi = lines_R2[1][7:-1][:8]
		cnt_umis[umi] = cnt_umis.get(umi, 0) + 1
		read_name = lines_R1[0][:-1].split(" ")
		read_name[0] = "{}:{}".format(read_name[0], umi)
		read_name = " ".join(read_name)
		out_fn[barcode].write("{}+{}\n".format(read_name, barcode))
		for i in range(1,4):
			out_fn[barcode].write(lines_R1[i])
	else: # In case barcode is not part of our list
		for i in range(4):
			no_barcode_match_R1.write(lines_R1[i])
			no_barcode_match_R2.write(lines_R2[i])


	# Read next read
	for i in range(4):
		lines_R1[i] = read_R1.readline()
		lines_R2[i] = read_R2.readline()

read_R1.close()
read_R2.close()

# Close fastq files
for b in barcodes:
	out_fn[b[2]].close()

no_barcode_match_R1.close()
no_barcode_match_R2.close()

# Output statistics
fout = open("{}/umis_stats.txt".format(output_dir), "w")
for u in cnt_umis.keys():
	fout.write("{}\t{}\n".format(u, cnt_umis[u]))
fout.close()

fout = open("{}/barcodes_stats.txt".format(output_dir), "w")
for b in cnt_barcodes.keys():
	fout.write("{}\t{}\n".format(b, cnt_barcodes[b]))
fout.close()
