from glob import glob
import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", dest = "d_input",
		type=str, required=True,
		                    help="Input dir with fastq file per well")
parser.add_argument("--order_file", dest = "order_fn",
		type=str, required=True,
		                    help="file that containg the names of each sample")
args = parser.parse_args()

order_file = [line.split("\t") for line in open(args.order_fn, "r").read().split("\n")[1:]]
order = {}
for line in order_file:
	well = line[0] + line[1]
	fastq_fn = glob("{}/{}_*.fastq".format(args.d_input, well))
	if len(fastq_fn) != 1:
		raise Exception("Expected to find 1 file for well {}, but found {}".format(well, len(fastq_fn)))
	fastq_fn = fastq_fn[0]
	new_dir = "{}/{}_{}".format(args.d_input, fastq_fn.split("/")[-1].split(".")[0], line[2])
	print("{}\t{}\t{}\n".format(well, fastq_fn.split("/")[-1],new_dir))
	sp.call(["mkdir", new_dir])
	sp.call(["mv", fastq_fn, new_dir+"/"])


