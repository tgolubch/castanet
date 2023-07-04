import os
import re
from app.utils.shell_cmds import shell
'''
Reading grouped HBV reads in experiments/HBV13_2/consensus_dbs/pHBV.fasta
Maps against pre-compiled database of HBV genomes in test/full_hbv.fasta#
'''


class Consensus:
    def __init__(self, SeqName, RefPathogen, GroupedAlignmentsFname) -> None:
        self.Bwa = './bwa-mem2-2.2.1_x64-linux/bwa-mem2'
        self.SeqName = SeqName
        self.RefPathogen = RefPathogen
        self.PathogenSeqDbs = "./pathogen_seq_dbs/"
        self.GroupedAlignmentsFname = GroupedAlignmentsFname
        self.GroupedAlignmentsFpath = f'experiments/{SeqName}/grouped_alignments/'
        self.RemapDir = f'{self.GroupedAlignmentsFpath}remaps/'
        self.ConsensusDir = f"{self.GroupedAlignmentsFpath}consensus_seqs/"
        self.ConsensusCandidates = f"{self.ConsensusDir}consensus_candidates_{self.RefPathogen}.fasta"
        if not os.path.isdir(self.PathogenSeqDbs):
            shell(f"mkdir {self.PathogenSeqDbs}")
        if not os.path.isdir(self.RemapDir):
            shell(f"mkdir {self.RemapDir}")
        if not os.path.isdir(self.ConsensusDir):
            shell(f"mkdir {self.ConsensusDir}")
            shell(f"touch {self.ConsensusCandidates}")
        else:
            shell(f"> {self.ConsensusCandidates}")

    def get_ref_db(self):
        # RM < TODO MAKE DYNAMIC
        return "pathogen_seq_dbs/hbv/all/all.fasta"

    def get_ref_seq(self, refseq):
        # RM < TODO MAKE DYNAMIC
        return f'hbv/{refseq.replace(".fasta", "")}/{refseq}'

    def reparse(self, refdb):
        '''Re-map to bam, stream SAM to new parse functions'''
        shell(f"{self.Bwa} mem {refdb} {self.GroupedAlignmentsFpath}{self.GroupedAlignmentsFname}.fasta | samtools view -F4 -Sb - | samtools sort - | samtools view -F2048 -F4 - > {self.GroupedAlignmentsFpath}{self.GroupedAlignmentsFname}.sam")
        shell(
            f"python3 -m app.src.parse_bam -Mode reparse -SeqName {self.SeqName}")

    def make_consensus(self, fname):
        print(f"Generating candidate consensus for: {fname}")
        # Check if we have a db for this pathogen
        refseq = self.get_ref_seq(fname)
        # map vs genome
        shell(f"{self.Bwa} mem {self.PathogenSeqDbs}{refseq} {self.RemapDir}{fname} | samtools view -F4 -Sb - | samtools sort - 1> {self.RemapDir}{fname}.bam")
        # build consensus
        shell(
            f"viral_consensus -i {self.RemapDir}{fname}.bam -r {self.PathogenSeqDbs}/{refseq} -o - >> {self.ConsensusCandidates}")

    def qa(self):
        # QA Consensus Genomes
        COV_THRESH = 0.8
        good_seqs = {}
        with open(self.ConsensusCandidates, "r") as f:
            for l in f:
                if l[0] == ">":
                    '''Header'''
                    find_acc = re.findall(
                        r"[A-Z]{1,2}[0-9]{5,6}|[A-Z]{4}[0-9]{6,8}|[A-Z]{2}_[0-9]{6}", l)
                    if find_acc:
                        assert len(
                            find_acc) > 1, f"Can't find two Accession IDs in the fasta header for line {l}"
                        assert find_acc[0] == find_acc[
                            1], f"Accession IDs were mismatched between input bam and ref for line {l}"
                        acc = find_acc[0]
                else:
                    '''Seq'''
                    cov = round((len(l) - l.count('N')) / len(l) * 100)
                    if cov * 0.01 > COV_THRESH:
                        good_seqs[acc] = l
                        print(
                            f"Accepting consensus for hit: {acc}, as genome coverage is {cov}%.")
                    else:
                        print(
                            f"Rejecting consensus for hit: {acc}, as genome coverage is {cov}%.")

        for key in good_seqs.keys():
            with open(f"{self.ConsensusDir}/consensus_{key}.fasta", "w") as file:
                file.write(f">{key}_consensus\n{good_seqs[key]}\n")

    def main(self):
        self.reparse(self.get_ref_db())
        if len(os.listdir(self.RemapDir)) == 0:
            print(
                f"Failed to remap grouped alignments for grouping: {self.RefPathogen} to reference genomes")
            return
        [self.make_consensus(i) for i in os.listdir(
            self.RemapDir) if not "bam" in i]
        self.qa()


if __name__ == "__main__":
    cls = Consensus(SeqName="HBV13_2",
                    RefPathogen="HBV",
                    GroupedAlignmentsFname="pHBV")
    cls.main()
