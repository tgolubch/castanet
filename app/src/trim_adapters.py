from app.utils.shell_cmds import shell

def run_trim(p, trim_path='java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar', api_entry=True):
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
            "NThreads": 4,
            "AdaptP": p.AdaptP
        }
    else: 
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"{trim_path} PE -threads {p['NThreads']} {p['ExpDir']}{p['SeqName']}_1_filt.fastq {p['ExpDir']}{p['SeqName']}_2_filt.fastq {p['ExpDir']}{p['SeqName']}_1_clean.fastq {p['ExpDir']}{p['SeqName']}_1_trimmings.fq {p['ExpDir']}{p['SeqName']}_2_clean.fastq {p['ExpDir']}{p['SeqName']}_2_trimmings.fq ILLUMINACLIP:{p['AdaptP']}:2:10:7:1:true MINLEN:80")
