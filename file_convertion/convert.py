import Bio
from Bio import SeqIO

# In here first install Bio 
# give the input file as your source file [line 12]
# * sim.contig.ans and sim.contig.nt files generate as output_file
# edit line [22] and [23] according to the source file 

# relevant index of [id name revese/forward length ] for sim.contig.ans
# relevant index of [id sequence] for sim.contig.nt 

input_file = "/sample_data/simulatedDS/454.10species.fasta"

f = open("sim.contig.ans", "w")
seq = open("sim.contig.nt", "w")

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    arr = name.split('|')
    f.write(arr[0] +'\t'+ arr[len(arr)-1] +'\t'+ arr[5] +'\t'+ arr[4]+'\n')
    seq.write(">" + arr[0] + "\n" +sequence + "\n")
f.close()