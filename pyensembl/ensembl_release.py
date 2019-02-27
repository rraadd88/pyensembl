# Copyright (c) 2015. Mount Sinai School of Medicine
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
Contains the EnsemblRelease class, which extends the Genome class
to be specific to (a particular release of) Ensembl.
"""
from weakref import WeakValueDictionary

from .genome import Genome
from .ensembl_release_versions import check_release_number, MAX_ENSEMBL_RELEASE
from .species import check_species_object
import logging
from .ensembl_url_templates import (
    ENSEMBL_FTP_SERVER,
    make_gtf_url,
    make_fasta_url
)

def get_dtype2url(name,release,test=False):
    import ftplib
    host='ftp.ensembl.org'
    ftp=ftplib.FTP(host,user='anonymous', passwd='anonymous')
    dtype=['gtf','cdna','ncrna','pep']
    subds=[f'gtf/{name}/',
            f'fasta/{name}/cdna/',
            f'fasta/{name}/ncrna/',
            f'fasta/{name}/pep/']
    dtype2subd=dict(zip(dtype,subds))
    dtype2url={}
    for dtype in dtype2subd:
        subd=dtype2subd[dtype]
        ext=f'/pub/release-{int(release)}/{subd}'
        if test:
            print(ext)
        ftp.cwd(ext)
        if subd.startswith('gtf'):
            file_fmt='gtf'
        elif subd.startswith('fasta'):
            file_fmt='fa'        
        fns=[p for p in ftp.nlst() if p.endswith(f'.{file_fmt}.gz') and (not 'abinitio' in p) and (not 'patch' in p) and (not 'scaff' in p)]
        fns=[fn for fn in fns if not '.chr.' in fn]
        if len(fns)>1:
            logging.warning(f'two files instead of one found for {dtype}')
            if test:
                print(fns)                
        url=f"ftp://{host}{ext}{fns[0]}"
        if test:
            print(url)
        dtype2url[dtype]=url
#         print(f'wget {dlp} -O pub --no-verbose --passive-ftp -T 3600')
    if test:
        print(dtype2url)
    return dtype2url

class EnsemblRelease(Genome):
    """
    Bundles together the genomic annotation and sequence data associated with
    a particular release of the Ensembl database.
    """
    @classmethod
    def normalize_init_values(cls, release, species, server):
        """
        Normalizes the arguments which uniquely specify an EnsemblRelease
        genome.
        """
        release = check_release_number(release)
        species = check_species_object(species)
        return (release, species, server)

    # Using a WeakValueDictionary instead of an ordinary dict to prevent a
    # memory leak in cases where we test many different releases in sequence.
    # When all the references to a particular EnsemblRelease die then that
    # genome should also be removed from this cache.
    _genome_cache = WeakValueDictionary()

    @classmethod
    def cached(
            cls,
            species=None,
            release=MAX_ENSEMBL_RELEASE,
            server=ENSEMBL_FTP_SERVER):
        """
        Construct EnsemblRelease if it's never been made before, otherwise
        return an old instance.
        """
        init_args_tuple = cls.normalize_init_values(release, species, server)
        if init_args_tuple in cls._genome_cache:
            genome = cls._genome_cache[init_args_tuple]
        else:
            genome = cls._genome_cache[init_args_tuple] = cls(*init_args_tuple)
        return genome

    def __init__(
            self,
            species=None,
            release=MAX_ENSEMBL_RELEASE,
            server=ENSEMBL_FTP_SERVER):
        # self.release, self.species, self.server = self.normalize_init_values(
        #     release=release, species=species, server=server)
        self.release, self.species, self.server = release, species,server
        try:
            # pyensembl way
            self.gtf_url = make_gtf_url(
                ensembl_release=self.release,
                species=self.species.latin_name,
                server=self.server)

            self.transcript_fasta_urls = [
                make_fasta_url(
                    ensembl_release=self.release,
                    species=self.species.latin_name,
                    sequence_type="cdna",
                    server=server),
                make_fasta_url(
                    ensembl_release=self.release,
                    species=self.species.latin_name,
                    sequence_type="ncrna",
                    server=server)
            ]

            self.protein_fasta_urls = [
                make_fasta_url(
                    ensembl_release=self.release,
                    species=self.species.latin_name,
                    sequence_type="pep",
                    server=self.server)]
        except:
            # my way
            logging.info('using ftplib to download genome')
            try:
                dtype2url=get_dtype2url(name=species.latin_name,release=release,
                                        test=False) 
                
            except:
                # if transient connection problem
                dtype2url=get_dtype2url(name=species.latin_name,release=release,
                                        test=False)            

            self.gtf_url = dtype2url['gtf']
            self.transcript_fasta_urls = [dtype2url['cdna'],dtype2url['ncrna']]
            self.protein_fasta_urls = [dtype2url['pep']]
        #common    
        self.reference_name = self.species.which_reference(self.release)
        
        Genome.__init__(
            self,
            reference_name=self.reference_name,
            annotation_name="ensembl",
            annotation_version=self.release,
            gtf_path_or_url=self.gtf_url,
            transcript_fasta_paths_or_urls=self.transcript_fasta_urls,
            protein_fasta_paths_or_urls=self.protein_fasta_urls)

    def install_string(self):
        return "pyensembl install --release %d --species %s" % (
            self.release,
            self.species.latin_name)

    def __str__(self):
        return "EnsemblRelease(release=%d, species='%s')" % (
            self.release,
            self.species.latin_name)

    def __eq__(self, other):
        return (
            other.__class__ is EnsemblRelease and
            self.release == other.release and
            self.species == other.species)

    def __hash__(self):
        return hash((self.release, self.species))

    def to_dict(self):
        return {
            "release": self.release,
            "species": self.species,
            "server": self.server
        }

    @classmethod
    def from_dict(cls, state_dict):
        """
        Deserialize EnsemblRelease without creating duplicate instances.
        """
        return cls.cached(**state_dict)


def cached_release(release, species):
    """
    Create an EnsemblRelease instance only if it's hasn't already been made,
    otherwise returns the old instance.
    Keeping this function for backwards compatibility but this functionality
    has been moving into the cached method of EnsemblRelease.
    """
    return EnsemblRelease.cached(release=release, species=species)
