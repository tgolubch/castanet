from fastapi import FastAPI, BackgroundTasks
from fastapi.encoders import jsonable_encoder
from app.utils.system_messages import banner
from app.utils.utility_fns import make_exp_dir
from app.src.preprocess import run_kraken
from app.src.filter_keep_reads import FilterKeepReads
from app.src.trim_adapters import run_trim
from app.src.map_reads_to_ref import run_map
from app.src.generate_counts import run_counts
from app.src.analysis import Analysis
from app.src.post_filter import run_post_filter
from app.utils.api_classes import (E2e_data, Preprocess_data, Filter_keep_reads_data,
                                   Trim_data, Mapping_data, Count_map_data, Analysis_data,
                                   Post_filter_data)

description = """
CASTANET is software for analysis of targeted metagenomics sequencing data, originally by tgolubch (https://github.com/tgolubch) and refactored to Python3 by mayne941 (https://github.com/Mayne941).
"""

tags_metadata = [
    {
        "name": "End to end pipeline",
        "description": "Run an end-to-end Castanet job",
    },
    {
        "name": "Individual pipeline functions",
        "description": "Run individual functions from the Castanet pipeline",
    },
    {
        "name": "Convenience functions",
        "description": "Supplementary functions for data processing",
    },
    {
        "name": "Dev endpoints",
        "description": "Developer tools"
    }
]

banner()

app = FastAPI(
    title="Castanet",
    version="0.2",
    description=description,
    contact={
        "name": "Nuffield Department of Medicine, University of Oxford",
        "url": "https://www.ndm.ox.ac.uk/",
    },
    license_info={  # RM < TODO CHECK LICENSE
        "name": "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "url": "https://creativecommons.org/licenses/by-nc/4.0/",
    },
    openapi_tags=tags_metadata
)

'''Dev Endpoints'''


@app.get("/", tags=["Dev endpoints"])
async def read_root():
    return {"response": "API is healthy. Append the current URL to include '/docs/' at the end to visit the GUI."}


'''Consumer endpoints'''


@app.post("/end_to_end/", tags=["End to end pipeline"])
async def end_to_end(payload: E2e_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_kraken(payload)
    do_filter_keep_reads(payload)
    run_trim(payload)
    run_map(payload)
    run_counts(payload)
    run_analysis(payload)
    if payload["PostFilt"]:
        run_post_filter(payload)
    return "Task complete. See terminal output for details."


@app.post("/preprocess/", tags=["Individual pipeline functions"])
async def preprocess(payload: Preprocess_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_kraken(payload)
    return "Task complete. See terminal output for details."


@app.post("/filter_keep_reads/", tags=["Individual pipeline functions"])
async def filter_keep_reads(payload: Filter_keep_reads_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    do_filter_keep_reads(payload)
    return "Task complete. See terminal output for details."


def do_filter_keep_reads(payload):
    cls = FilterKeepReads(payload)
    cls.main()


@app.post("/trim_data/", tags=["Individual pipeline functions"])
async def trim_data(payload: Trim_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_trim(payload)
    return "Task complete. See terminal output for details."


@app.post("/mapping/", tags=["Individual pipeline functions"])
async def mapping(payload: Mapping_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_map(payload)
    return "Task complete. See terminal output for details."


@app.post("/generate_counts/", tags=["Individual pipeline functions"])
async def count_mapped(payload: Count_map_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_counts(payload)
    return "Task complete. See terminal output for details."


@app.post("/analysis/", tags=["Individual pipeline functions"])
async def analysis(payload: Analysis_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_analysis(payload)
    return "Task complete. See terminal output for details."


def run_analysis(payload):
    cls = Analysis(payload)
    cls.main()


@app.post("/post_filter/", tags=["Individual pipeline functions"])
async def post_filter(payload: Post_filter_data):
    payload = jsonable_encoder(payload)
    make_exp_dir(payload["ExpName"])
    run_post_filter(payload)
    return "Task complete. See terminal output for details."

'''Conveninece endpoints'''
...
