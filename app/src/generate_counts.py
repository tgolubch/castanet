from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print


def run_counts(p, api_entry=True):
    if not api_entry:
        p = {
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName,
        }
    else:
        p["ExpDir"] = f"{p['ExpDir']}/"
    shell(f"""for BamFilePath in $(ls experiments/{p['ExpName']}/*.bam); do BamPath=$BamFilePath; BamName=$(basename "BamPath%%.bam"); BamName=$(sed s'/_dedup//' <<< ${{BamName}}); samtools view -F2048 -F4 ${{BamPath}} | python3 -m app.src.parse_bam -Mode parse -SeqName {p['SeqName']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g; done > experiments/{p['ExpName']}/{p['SeqName']}_PosCounts.csv""")
    end_sec_print("INFO: Counts generated")
