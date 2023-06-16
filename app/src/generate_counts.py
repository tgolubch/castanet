from app.utils.shell_cmds import shell

def run_counts(p, api_entry=True):
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"""for BamFilePath in $(ls {p['ExpDir']}/*.bam); do BamPath=$BamFilePath; BamName=$(basename "BamPath%%.bam"); BamName=$(sed s'/_dedup//' <<< ${{BamName}}); samtools view -F2048 -F4 ${{BamPath}} | python3 src/parse_bam_positions.py {p['SeqName']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g; done > {p['ExpDir']}PosCounts.csv""")
    print("\nINFO: Counts generated, no errors reported by Samtools\n")