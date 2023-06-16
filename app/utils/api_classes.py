from fastapi import Query
from typing import Union, Literal, Optional
from pydantic import BaseModel, Field, FilePath, DirectoryPath, validator

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
    
class Data_ExpName(BaseModel):
    ExpName: str = Query('myexperiment',
                         description="Name your experiment/batch")
    
class Data_KrakenDir(BaseModel):
    KrakenDbDir: DirectoryPath = Query('kraken2_human_db/',
                                     description="Path to Kraken2 database for filtering human/other unwanted species reads.")
    
class Data_NThreads(BaseModel):
    NThreads: int = Query(4,
                                     description="Number of threads to doing parallel jobs on, with supported dependencies. Cannot exceed n logical CPU cores.")

class Data_FilterFilters(BaseModel):
    LineageFile: Union[None, str] = Query('data/ncbi_lineages_2023-06-15.csv.gz',
                                      description="(OPTIONAL) Path to CSV file containing lineages of all NCBI taxa.")

    ExcludeIds: Union[None, str] = Query("1969841",
                                      description="(OPTIONAL) Exclude these NCBI TaxID/s from filter keep reads step. Comma separate without spaces. Set to 9606 to exclude Human.")
    
    RetainIds: Union[None, str] = Query("Delete this line if not using",
                                      description="(OPTIONAL) Exclude these NCBI TaxID/s from filter keep reads step. Comma separate without spaces.")

    RetainNames: Union[None, str] = Query("Delete this line if not using",
                                      description="(OPTIONAL) Retain these species names from filter keep reads step. Comma separate without spaces. Will be ignored if no Linneage file specified.")

    ExcludeNames: Union[None, str] = Query("Homo,Alteromonas,Achromobacter",
                                      description="(OPTIONAL) Exclude these species names from filter keep reads step. Comma separate without spaces. Will be ignored if no Linneage file specified.")

class Data_AnalysisExtras(BaseModel):
    Probes: str = Query("data/probelengths_rmlst_virus_extra_ercc.csv",
                        description="CSV file containing probe length mappings. Absolute path required.")
    Samples: str = Query("data/samples.csv",
                         description="CSV file containing sample data for annotations during analysis phase. Absolute path required.")
    KeepDups: bool = Query(True,
                           description='(OPTIONAL) If true, do not reassign duplicates to the sample with the majority in each duplicate cluster (Default: True).')
    Clin: Optional[str] = Query("Delete this line if not using",
                                description='(OPTIONAL) Path to CSV file containing clinical data (must have at least following fields: pt, clin_int; the field "sampleid" if present will be ignored). Other fields will be ignored.')
    DepthInf: str = Query("Delete this line if not using",
                          description='(OPTIONAL, For regenerating full CSV with new clinical info): Path to previously generated CSV file of read depth per position for each probe, for all samples in this batch. Must contain the following fields: sampleid, target_id, depth_mean, depth_std, depth_25pc, depth_median, depth_75pc, prop_target_covered, prop_target_covered_mindepth2, prop_target_covered_mindepth5, prop_target_covered_mindepth10, udepth_mean, udepth_std, udepth_25pc, udepth_median, udepth_75pc, uprop_target_covered, uprop_target_covered_mindepth2, uprop_target_covered_mindepth5, uprop_target_covered_mindepth10')

'''Endpoint objects'''

class E2e_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_AdaptP, Data_RefStem, 
                Data_PostFilt, Data_AnalysisExtras, Data_KrakenDir, Data_NThreads, Data_FilterFilters):
    ...

class Preprocess_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_KrakenDir, Data_NThreads):
    ...

class Filter_keep_reads_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_FilterFilters):
    ...

class Trim_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_AdaptP, Data_NThreads):
    ...

class Mapping_data(Data_ExpDir, Data_SeqName, Data_ExpName, Data_RefStem):
    ...

class Count_map_data(Data_ExpDir, Data_SeqName, Data_ExpName):
    ...

class Analysis_data(Data_ExpDir, Data_ExpName, Data_AnalysisExtras):
    ...

class Post_filter_data(Data_ExpDir, Data_SeqName, Data_ExpName):
    ...