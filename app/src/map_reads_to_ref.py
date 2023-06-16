from app.utils.shell_cmds import shell 
from app.utils.system_messages import end_sec_print

def run_map(p, bwa_path='./bwa-mem2-2.2.1_x64-linux/bwa-mem2', api_entry=True):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
            "RefStem": p.RefStem
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"{bwa_path} index {p['ExpDir']}{p['RefStem']}")
    shell(f"{bwa_path} mem {p['ExpDir']}{p['RefStem']} {p['ExpDir']}{p['SeqName']}_[12]_clean.fastq | samtools view -F4 -Sb - | samtools sort - 1> {p['ExpDir']}{p['SeqName']}.bam")
    end_sec_print(f"INFO: BWA mapping complete")