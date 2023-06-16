from app.utils.shell_cmds import shell

def run_kraken(p, api_entry=True):
    '''Call Kraken2 to remove unwanted reads'''
    if not api_entry:
        p = {
            "KrakenDbDir": "kraken2_human_db/",
            "NThreads": 4,
            "ExpDir": p.ExpDir,
            "SeqName": p.SeqName
        }
    shell(f'kraken2 --db {p["KrakenDbDir"]} --threads {p["NThreads"]} {p["ExpDir"]}{p["SeqName"]}_1.fastq.gz > {p["ExpDir"]}{p["SeqName"]}_1.kraken')
