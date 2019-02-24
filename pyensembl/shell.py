# Copyright (c) 2015-2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Manipulate pyensembl's local cache.

    %(prog)s {install, delete, delete-sequence-cache} [--release XXX --species human...]

To install particular Ensembl human release(s):
    %(prog)s install --release 75 77

To install particular Ensembl mouse release(s):
    %(prog)s install --release 75 77 --species mouse

To delete all downloaded and cached data for a particular Ensembl release:
    %(prog)s delete-all-files --release 75 --species human

To delete only cached data related to transcript and protein sequences:
    %(prog)s delete-index-files --release 75

To list all installed genomes:
    %(prog)s list

To install a genome from source files:
    %(prog)s install \
 --reference-name "GRCh38" \
 --gtf URL_OR_PATH \
 --transcript-fasta URL_OR_PATH \
 --protein-fasta URL_OR_PATH
"""

from __future__ import print_function, division, absolute_import

import argparse
import logging
import logging.config
import pkg_resources
import os

from .ensembl_release import EnsemblRelease, MAX_ENSEMBL_RELEASE
from .genome import Genome
from .species import Species

logging.config.fileConfig(pkg_resources.resource_filename(__name__, "logging.conf"))
logger = logging.getLogger(__name__)


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument(
    "--overwrite",
    default=False,
    action="store_true",
    help="Force download and indexing even if files already exist locally")


root_group = parser.add_mutually_exclusive_group()

release_group = root_group.add_argument_group()
release_group.add_argument(
    "--release",
    type=int,
    nargs="+",
    default=[],
    help="Ensembl release version(s) (default=%d)" % MAX_ENSEMBL_RELEASE)

release_group.add_argument(
    "--species",
    default=[],
    nargs="+",
    help="Which species to download Ensembl data for (default=human)")

release_group.add_argument(
    "--custom-mirror",
    default=None,
    help="URL and directory to use instead of the default Ensembl FTP server")

path_group = root_group.add_argument_group()

path_group.add_argument(
    "--reference-name",
    type=str,
    default=None,
    help="Name of the reference, e.g. GRCh38")

path_group.add_argument(
    "--annotation-name",
    default=None,
    help="Name of annotation source (e.g. refseq)")

path_group.add_argument(
    "--annotation-version",
    default=None,
    help="Version of annotation database")

path_group.add_argument(
    "--gtf",
    type=str,
    default=None,
    help="URL or local path to a GTF file containing annotations.")

path_group.add_argument(
    "--transcript-fasta",
    type=str,
    action='append',
    default=[],
    help="URL or local path to a FASTA files containing the transcript "
    "data. This option can be specified multiple times for multiple "
    "FASTA files.")

path_group.add_argument(
    "--protein-fasta",
    type=str,
    default=[],
    action="append",
    help="URL or local path to a FASTA file containing protein data.")

path_group.add_argument(
    "--annotation_name",
    type=str,
    default=[],
    action="append",
    help="assembly name?.")

path_group.add_argument(
    "--shared-prefix",
    default="",
    help="Add this prefix to URLs or paths specified by --gtf, --transcript-fasta, --protein-fasta")

actions=("install",
        "delete-all-files",
        "delete-index-files",
        "list")
parser.add_argument(
    "action",
    type=lambda arg: arg.lower().strip(),
    choices=actions,
    help=(
        "\"install\" will download and index any data that is  not "
        "currently downloaded or indexed. \"delete-all-files\" will delete all data "
        "associated with a genome annotation. \"delete-index-files\" deletes "
        "all files other than the original GTF and FASTA files for a genome. "
        "\"list\" will show you all installed Ensembl genomes."))

# def collect_all_installed_ensembl_releases():
#     # genomes = []
#     import pandas as pd
#     # dataframe that contains all the info (and can be exported as a tsv).
#     dspecies=pd.DataFrame(columns=['latin name','release','synonymn','assembly'])  
#     genomei=0  
#     for (species, release,assembly) in Species.all_species_release_pairs():
#         dspecies.loc[genomei,'assembly']=assembly
#         dspecies.loc[genomei,'release']=release
#         dspecies.loc[genomei,'synonymn']=species
#         dspecies.loc[genomei,'latin name']=species
#         genomei+=1
#     return dspecies
    #     species_obj=Species.register(
    #     latin_name=species,
    #     synonyms=[species],
    #     reference_assemblies={
    #         assembly: (release, release),
    #     })            
    # species_obj = Species.register(
    #     latin_name="saccharomyces_cerevisiae",
    #     synonyms=["yeast"],
    #     reference_assemblies={
    #         "R64-1-1": (76, MAX_ENSEMBL_RELEASE),
    #     })        
    # genome = EnsemblRelease(species=species_obj,release=release)
    # if genome.required_local_files_exist():
    #     genomes.append(genome)
    # return sorted(genomes, key=lambda g: (g.species.latin_name, g.release))


# def all_combinations_of_ensembl_genomes(args):
#     """
#     Use all combinations of species and release versions specified by the
#     commandline arguments to return a list of EnsemblRelease or Genome objects.
#     The results will typically be of type EnsemblRelease unless the
#     --custom-mirror argument was given.
#     """
#     species_list = args.species if args.species else ["human"]
#     release_list = args.release if args.release else [MAX_ENSEMBL_RELEASE]
#     genomes = []
#     for species in species_list:
#         # Otherwise, use Ensembl release information
#         for version in release_list:
#             ensembl_release = EnsemblRelease(species=species,release=version)
#             if not args.custom_mirror:
#                 genomes.append(ensembl_release)
#             else:
#                 # if we're using a custom mirror then we expect the provided
#                 # URL to be a directory with all the same filenames as
#                 # would be provided by Ensembl
#                 gtf_url = os.path.join(
#                     args.custom_mirror,
#                     os.path.basename(ensembl_release.gtf_url))
#                 transcript_fasta_urls = [
#                     os.path.join(
#                         args.custom_mirror,
#                         os.path.basename(transcript_fasta_url))
#                     for transcript_fasta_url in ensembl_release.transcript_fasta_urls
#                 ]
#                 protein_fasta_urls = [
#                     os.path.join(
#                         args.custom_mirror,
#                         os.path.basename(protein_fasta_url))
#                     for protein_fasta_url in
#                     ensembl_release.protein_fasta_urls
#                 ]
#                 reference_name = ensembl_release.reference_name
#                 genome = Genome(
#                     reference_name=reference_name,
#                     annotation_name="ensembl",
#                     annotation_version=version,
#                     gtf_path_or_url=gtf_url,
#                     transcript_fasta_paths_or_urls=transcript_fasta_urls,
#                     protein_fasta_paths_or_urls=protein_fasta_urls)
#                 genomes.append(genome)
#     return genomes

# def collect_selected_genomes(args):
#     # If specific genome source URLs are provided, use those
#     if args.gtf or args.transcript_fasta or args.protein_fasta:
#         if args.release:
#             raise ValueError(
#                 "An Ensembl release cannot be specified if "
#                 "specific paths are also given")
#         if not args.reference_name:
#             raise ValueError("Must specify a reference name")
#         if not args.annotation_name:
#             raise ValueError("Must specify the name of the annotation source")

#         return [Genome(
#             reference_name=args.reference_name,
#             annotation_name=args.annotation_name,
#             annotation_version=args.annotation_version,
#             gtf_path_or_url=os.path.join(args.shared_prefix, args.gtf),
#             transcript_fasta_paths_or_urls=[
#                 os.path.join(args.shared_prefix, transcript_fasta)
#                 for transcript_fasta in args.transcript_fasta],
#             protein_fasta_paths_or_urls=[
#                 os.path.join(args.shared_prefix, protein_fasta)
#                 for protein_fasta in args.protein_fasta])]
#     else:        
#         return all_combinations_of_ensembl_genomes(args)

def get_ensembl_genome_object(args):
    """
    Use all combinations of species and release versions specified by the
    commandline arguments to return a list of EnsemblRelease or Genome objects.
    The results will typically be of type EnsemblRelease unless the
    --custom-mirror argument was given.
    """
    species = args.species[0]  
    version = args.release[0] if args.release else [MAX_ENSEMBL_RELEASE]
    assembly=args.reference_name
    species_obj=Species.register(
    latin_name=species,
    synonyms=[species],
    reference_assemblies={
        assembly: (version, version),
    })
    ensembl_release = EnsemblRelease(species=species_obj,release=version)
    if not args.custom_mirror:
        return ensembl_release

def args2genome(args):
    # AKA collect_selected_genome (no 's')
    if args.gtf or args.transcript_fasta or args.protein_fasta:
        if args.release:
            raise ValueError(
                "An Ensembl release cannot be specified if "
                "specific paths are also given")
        if not args.reference_name:
            raise ValueError("Must specify a reference name")
        if not args.annotation_name:
            raise ValueError("Must specify the name of the annotation source")

        return Genome(
            reference_name=args.reference_name,
            annotation_name=args.annotation_name,
            annotation_version=args.annotation_version,
            gtf_path_or_url=os.path.join(args.shared_prefix, args.gtf),
            transcript_fasta_paths_or_urls=[
                os.path.join(args.shared_prefix, transcript_fasta)
                for transcript_fasta in args.transcript_fasta],
            protein_fasta_paths_or_urls=[
                os.path.join(args.shared_prefix, protein_fasta)
                for protein_fasta in args.protein_fasta])
    else:
        return get_ensembl_genome_object(args)

def get_genome_files(args):
    from glob import glob
    from os.path import dirname,basename
    from .species import normalize_species_name

    species = normalize_species_name(args.species[0])  
    species_=species[0].upper()+species[1:]

    release = args.release[0] if args.release else [MAX_ENSEMBL_RELEASE]
    release="ensembl"+str(release)
    assembly=args.reference_name

    # here's how I get the .cache directory eg. '/home/user/.cache/pyensembl'
    import datacache
    pyensembl_cache_dir=f"{dirname(datacache.get_data_dir())}/pyensembl" #FIXME if genomes are installed at other places than .cache
    # all the assemblies
    genome_=f"{pyensembl_cache_dir}/{assembly}/{release}/{species_}"
    genome_files=glob(f"{genome_}*")
    return genome_files,genome_

def run():
    args = parser.parse_args()
    if args.action == "list":
        # TODO: how do we also identify which non-Ensembl genomes are
        # installed?
        # dspecies = collect_all_installed_ensembl_releases()
        print(Species.dspecies)
        # for genome in genomes:
        #     # print every directory in which downloaded files are located
        #     # in most case this will be only one directory
        #     filepaths = genome.required_local_files()
        #     directories = {os.path.split(path)[0] for path in filepaths}
        #     print("-- %s: %s" % (genome, ", ".join(directories)))
    else:
        if args.action in actions:
            if args.species:
                species = args.species[0]  
            else: 
                raise ValueError(f"Species argument provided is {args.species}")    
            logger.info("Running '%s' for %s", args.action, args.species)
            genome_files,genome_=get_genome_files(args)
            if len(genome_files)>1:
                genome = args2genome(args)
                if args.action == "delete-all-files":
                    genome.download_cache.delete_cache_directory()
                elif args.action == "delete-index-files":
                    genome.delete_index_files()
            else:
                from os import makedirs
                from os.path import dirname
                # makedirs(dirname(genome_),exist_ok=True)
                # open(genome_, 'w').close()
                # print(genome_)
                genome = args2genome(args)
                # except:
                if args.action == "install":   
                    genome.download(overwrite=args.overwrite)
                    genome.index(overwrite=args.overwrite)
        else:
            raise ValueError("Invalid action: %s" % args.action)
