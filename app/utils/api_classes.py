from fastapi import Query
from typing import Union, Literal
from pydantic import BaseModel, Field, FilePath, DirectoryPath

'''Primitives'''

class Data_ExpDir(BaseModel):
    ExpDir: DirectoryPath = Query('./data/',
                                     description="Path to retrieve and store data.")
    
class Data_SeqName(BaseModel):
    SeqName: str = Query('mysequence',
                          description="Base filename for your input sequences. Naming convention is mysequence_1, ..._2.")

class Data_AdaptP(BaseModel):
    AdaptP: str = Query('Trimmomatic-0.39/adapters/all.fa',
                        description='Location of your Trimmomatic adapter sequences - may be in your Trimmomatic path, but a backup is in the data dir.')

class Data_RefStem(BaseModel):
    RefStem: str = Query("rmlst_virus_extra_ercc_2018.fasta",
                    description="Path to mapping file, in fasta format. Assumes lives in data dir.")

class Data_PostFilt(BaseModel):
    PostFilt: bool= Query(False,
                         description="Post hoc filter BAM file to remove reads marked as contaminations")

class Data_Samples(BaseModel):
    Samples: str = Query("data/samples.csv",
                         description="CSV file containing sample data for annotations during analysis phase. Absolute path required.")

class Data_Probes(BaseModel):
    Probes: str = Query("data/probelengths_rmlst_virus_extra_ercc.csv",
                        description="CSV file containing probe length mappings. Absolute path required.")
    
class Data_ExpName(BaseModel):
    ExpName: str = Query('myexperiment',
                         description="Name your experiment/batch")
    
class Data_KrakenDir(BaseModel):
    KrakenDbDir: DirectoryPath = Query('kraken2_human_db/',
                                     description="Path to Kraken2 database for filtering human/other unwanted species reads.")
    
class Data_NThreads(BaseModel):
    NThreads: int = Query(4,
                                     description="Number of threads to doing parallel jobs on, with supported dependencies. Cannot exceed n logical CPU cores.")


'''Endpoint objects'''

class E2e_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_AdaptP, Data_RefStem, 
                Data_PostFilt, Data_Probes, Data_KrakenDir, Data_NThreads):
    ...

class Preprocess_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_KrakenDir, Data_NThreads):
    ...

class Filter_keep_reads_data(Data_ExpDir, Data_SeqName, Data_ExpName):
    ...

class Trim_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_AdaptP):
    ...

class Mapping_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_RefStem):
    ...

class Count_map_data(Data_ExpDir, Data_SeqName, Data_ExpName):
    ...

class Analysis_data(Data_ExpName, Data_Samples, Data_Probes):
    ...

class Post_filter_data(Data_ExpDir, Data_SeqName, Data_ExpName):
    ...