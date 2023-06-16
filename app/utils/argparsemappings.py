def map_args_filter_keep_reads(ap):
    '''Map arg parser arguments to API arguments'''
    opts = []
    for field in [ap.r, ap.x, ap.rT, ap.xT, ap.lineagefile]:
        try:
            opts.append(field)
        except AttributeError:
            opts.append("")
    return {
        "input_file":  ap.i,
        "kraken": ap.k,
        "RetainIds": opts[0], 
        "ExcludeIds": opts[1], 
        "RetainNames": opts[2],
        "ExcludeNames": opts[3],
        "LineageFile": opts[4]
    }