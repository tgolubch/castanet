from app.utils.shell_cmds import shell

def run_post_filter(p, api_entry=True):
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
            "ExpName": p.ExpName
        }
    else: 
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"samtools view -h {p['ExpDir']}{p['SeqName']}.bam | python3 -m app.src.parse_bam -Mode filter -ExpDir {p['ExpDir']} -SeqName {p['SeqName']} -FilterFile {p['ExpName']}_reads_to_drop.csv")
    print(f"\nINFO: Post filter complete\n")