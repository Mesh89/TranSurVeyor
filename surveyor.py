import argparse, os
import pysam, pyfaidx
from random_pos_generator import RandomPositionGenerator
import numpy as np

MAX_READS = 1000
GEN_DIST_SIZE = 100000
MAX_ACCEPTABLE_IS = 20000

cmd_parser = argparse.ArgumentParser(description='TranSurVeyor, a transposition caller.')
cmd_parser.add_argument('bamFile', help='Input bam file.')
cmd_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
cmd_parser.add_argument('reference', help='Reference genome in FASTA format.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--samtools', default='samtools', help='Samtools path.')
cmd_parser.add_argument('--bwa', default='bwa', help='BWA path.')
cmd_parser.add_argument('--maxSCDist', type=int, default=5, help='Max distance (in bp) for two clips to be considered '
                                                                  'representing the same breakpoint.')
cmd_parser.add_argument('--maxTRAsize', type=int, default=10000, help='Maximum size for which TranSurVeyor will try to '
                                                                     'predict both breakpoints.')
cmd_parser.add_argument('--maxIPdist', type=int, default=1000, help='Maximum distance between forward and reverse '
                                                                    'insertion points.')


cmd_args = cmd_parser.parse_args()


# Check if workdir is empty

if os.listdir(cmd_args.workdir):
    print "Workdir is not empty!"
    exit(1)


# Create config file in workdir

config_file = open(cmd_args.workdir + "/config.txt", "w")
config_file.write("threads %d\n" % cmd_args.threads)
config_file.write("max_sc_dist %d\n" % cmd_args.maxSCDist)
config_file.write("max_tra_size %d\n" % cmd_args.maxTRAsize)
config_file.write("max_ip_dist %d\n" % cmd_args.maxIPdist)

TRANSURVEYOR_PATH = os.path.dirname(os.path.realpath(__file__))

# Find read length

read_len = 0
bam_file = pysam.AlignmentFile(cmd_args.bamFile)
for i, read in enumerate(bam_file.fetch(until_eof=True)):
    if i > MAX_READS: break
    read_len = max(read_len, read.query_length)
config_file.write("read_len %d\n" % read_len)


# Generate general distribution of insert sizes

contig_map = open("%s/contig_map" % cmd_args.workdir, "w")
num_contigs = 0

reference_fa = pyfaidx.Fasta(cmd_args.reference)
i = 1
for k in reference_fa.keys():
    contig_map.write("%s %d\n" % (k, i));
    i += 1
    num_contigs += 1
contig_map.close();

rand_pos_gen = RandomPositionGenerator(reference_fa)
random_positions = []
for i in range(1,1000001):
    if i % 100000 == 0: print i, "random positions generated."
    random_positions.append(rand_pos_gen.next())

with open("%s/random_pos.txt" % cmd_args.workdir, "w") as random_pos_file:
    for random_pos in random_positions:
        random_pos_file.write("%s %d\n" % random_pos)

general_dist = []
avg_depth = 0
samplings = 0
rnd_i = 0
while len(general_dist) < GEN_DIST_SIZE:
    chr, pos = random_positions[rnd_i]
    rnd_i += 1

    if pos > len(reference_fa[chr])-10000:
        continue

    samplings += 1
    i = 0
    for read in bam_file.fetch(contig=chr, start=pos, end=pos+10000):
        if read.reference_start < pos:
            avg_depth += 1
        if read.is_proper_pair and not read.is_secondary and not read.is_supplementary and \
        0 < read.template_length < MAX_ACCEPTABLE_IS and 'S' not in read.cigarstring and 'S' not in read.get_tag('MC'):
            if i > 100: break
            i += 1
            general_dist.append(read.template_length)

mean_is = np.mean(general_dist)
stddev_is = np.std(general_dist)
avg_depth = float(avg_depth)/samplings;

print "Average depth:", avg_depth

general_dist = [x for x in general_dist if abs(x-mean_is) < 5*stddev_is]

mean_is = int(np.mean(general_dist))
lower_stddev_is = int(np.sqrt(np.mean([(mean_is-x)**2 for x in general_dist if x < mean_is])))
higher_stddev_is = int(np.sqrt(np.mean([(x-mean_is)**2 for x in general_dist if x > mean_is])))

min_is, max_is = mean_is-3*lower_stddev_is, mean_is+3.5*higher_stddev_is
print mean_is, min_is, max_is

config_file.write("min_is %d\n" % min_is)
config_file.write("avg_is %d\n" % mean_is)
config_file.write("max_is %d\n" % max_is)


workspace = cmd_args.workdir + "/workspace"

config_file.write("read_len %d\n" % read_len)
config_file.write("avg_depth %f\n" % avg_depth)
config_file.close();

if not os.path.exists(workspace):
    os.makedirs(workspace)

read_categorizer_cmd = TRANSURVEYOR_PATH + "/reads_categorizer %s %s" % (cmd_args.bamFile, cmd_args.workdir);
print "Executing:", read_categorizer_cmd
os.system(read_categorizer_cmd)


clip_consensus_builder_cmd = TRANSURVEYOR_PATH + "/clip_consensus_builder %s %s" % (cmd_args.workdir, cmd_args.reference)
print "Executing:", clip_consensus_builder_cmd
os.system(clip_consensus_builder_cmd)


dc_remapper_cmd = TRANSURVEYOR_PATH + "/dc_remapper %s %s" % (cmd_args.workdir, cmd_args.reference)
print "Executing:", dc_remapper_cmd
os.system(dc_remapper_cmd)


add_filtering_info_cmd = TRANSURVEYOR_PATH + "/add_filtering_info %s %s" % (cmd_args.bamFile, cmd_args.workdir)
print "Executing:", add_filtering_info_cmd
os.system(add_filtering_info_cmd)
