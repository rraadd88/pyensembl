# Copyright (c) 2015-2016. Mount Sinai School of Medicine
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

from __future__ import print_function, absolute_import, division

from serializable import Serializable

from .ensembl_release_versions import MAX_ENSEMBL_RELEASE


class Species(Serializable):
    """
    Container for combined information about a species name, its synonyn names
    and which reference to use for this species in each Ensembl release.
    """

    # as species instances get created, they get registered in these
    # dictionaries
    _latin_names_to_species = {}
    _common_names_to_species = {}
    _reference_names_to_species = {}

    @classmethod
    def register(cls, latin_name, synonyms, reference_assemblies):
        """
        Create a Species object from the given arguments and enter into
        all the dicts used to look the species up by its fields.
        """
        species = Species(
            latin_name=latin_name,
            synonyms=synonyms,
            reference_assemblies=reference_assemblies)
        cls._latin_names_to_species[species.latin_name] = species
        for synonym in synonyms:
            if synonym in cls._common_names_to_species:
                raise ValueError("Can't use synonym '%s' for both %s and %s" % (
                    synonym,
                    species,
                    cls._common_names_to_species[synonym]))
            cls._common_names_to_species[synonym] = species
        for reference_name in reference_assemblies:
            if reference_name in cls._reference_names_to_species:
                raise ValueError("Can't use reference '%s' for both %s and %s" % (
                    reference_name,
                    species,
                    cls._reference_names_to_species[reference_name]))
            cls._reference_names_to_species[reference_name] = species
        return species

    @classmethod
    def all_registered_latin_names(cls):
        """
        Returns latin name of every registered species.
        """
        return list(cls._latin_names_to_species.keys())

    @classmethod
    def all_registered_assemblies(cls):
        """
        Returns latin name of every registered species.
        """
        return list(cls._reference_names_to_species.keys())

    @classmethod
    def all_species_release_pairs(cls):
        """
        Generator which yields (species, release) pairs
        for all possible combinations.
        """
        for assembly in cls.all_registered_assemblies():
            for species_name in cls.all_registered_latin_names():
                species = cls._latin_names_to_species[species_name]
                for _, release_range in species.reference_assemblies.items():
                    for release in range(release_range[0], release_range[1] + 1):
                        yield species_name, release,assembly

    def __init__(self, latin_name, synonyms=[], reference_assemblies={}):
        """
        Parameters
        ----------
        latin_name : str

        synonyms : list of strings

        reference_assemblies : dict
            Mapping of names of reference genomes onto inclusive ranges of
            Ensembl releases Example: {"GRCh37": (54, 75)}
        """
        self.latin_name = latin_name.lower().replace(" ", "_")
        self.synonyms = synonyms
        self.reference_assemblies = reference_assemblies
        self._release_to_genome = {}
        for (genome_name, (start, end)) in self.reference_assemblies.items():
            for i in range(start, end + 1):
                assert i not in self._release_to_genome, \
                    "Ensembl release %d already has an associated genome"
                self._release_to_genome[i] = genome_name

    def which_reference(self, ensembl_release):
        if ensembl_release not in self._release_to_genome:
            raise ValueError("No genome for %s in Ensembl release %d" % (
                self.latin_name, ensembl_release))
        return self._release_to_genome[ensembl_release]

    def __str__(self):
        return (
            "Species(latin_name='%s', synonyms=%s, reference_assemblies=%s)" % (
                self.latin_name, self.synonyms, self.reference_assemblies))

    def __eq__(self, other):
        return (
            other.__class__ is Species and
            self.latin_name == other.latin_name and
            self.synonyms == other.synonyms and
            self.reference_assemblies == other.reference_assemblies)

    def to_dict(self):
        return {"latin_name": self.latin_name}

    @classmethod
    def from_dict(cls, state_dict):
        return cls._latin_names_to_species[state_dict["latin_name"]]

    def __hash__(self):
        return hash((self.latin_name,
                     tuple(self.synonyms),
                     frozenset(self.reference_assemblies.items())))


def normalize_species_name(name):
    """
    If species name was "Homo sapiens" then replace spaces with underscores
    and return "homo_sapiens". Also replace common names like "human" with
    "homo_sapiens".
    """
    lower_name = name.lower().strip()

    # if given a common name such as "human", look up its latin equivalent
    if lower_name in Species._common_names_to_species:
        return Species._common_names_to_species[lower_name].latin_name

    return lower_name.replace(" ", "_")


def find_species_by_name(species_name):
    # for (species, release) in Species.all_species_release_pairs():
    #     print(species, release)
    latin_name = normalize_species_name(species_name)
    if latin_name not in Species._latin_names_to_species:
        raise ValueError("Species not found: %s" % species_name)
    return Species._latin_names_to_species[latin_name]


def check_species_object(species_name_or_object):
    """
    Helper for validating user supplied species names or objects.
    """
    if isinstance(species_name_or_object, Species):
        return species_name_or_object
    elif isinstance(species_name_or_object, str):
        return find_species_by_name(species_name_or_object)
    else:
        raise ValueError("Unexpected type for species: %s : %s" % (
            species_name_or_object, type(species_name_or_object)))

# start of a PR (#207) from @rraadd88 
def collect_all_genomes():
    """
    data aware generation of Species object.
    searches in .cache dir and generates a Species object

    Also generates a tsv file with all the genome info.
    Such file can be used to install sets of genomes at once (logic is bit like conda environment profile).
    It would be relatively easy to code export and import for such a file.
    """

    def str2num(s,cat=False,force=True):
        """
        Converts string to integer
        eg. ensembl92 to 92

        :param s: string
        :param cat: Whether to concatenate detected integers. eg. 20,23 to 2023
        :param force: If True, ignores decimal point error. 
        """
        import re    
        if '.' in s and not force:
            raise ValueError(f"A string can only be converted to integeres, found a '.' in {s}")
        n=re.findall(r'\d+',s)
        if len(n)==0:
            raise ValueError(f"No digits found in string {s}")        
        elif len(n)==1:
            return int(n[0])
        else:
            if cat:
                return int(''.join(n))
            else:
                return n

    from glob import glob
    from os.path import dirname,basename,exists
    import numpy as np
    import pandas as pd
    from pyensembl.species import normalize_species_name,Species
            
    # here's how I get the .cache directory eg. '/home/user/.cache/pyensembl'
    import datacache
    pyensembl_cache_dir=f"{dirname(datacache.get_data_dir())}/pyensembl" #FIXME if genomes are installed at other places than .cache

    # all the assemblies
    assemblies=[basename(p) for p in glob(f"{pyensembl_cache_dir}/*")]
    # dataframe that contains all the info (and can be exported as a tsv).
    dspecies=pd.DataFrame(columns=['latin name','release','synonymn','assembly'])
    # assempy to release min max dict needed as an input to create Species object
    assembly2releasesminmax={}
    # following loop populates the dataframe 
    genomei=0
    for assembly in assemblies:
        releases=[basename(p) for p in glob(f"{pyensembl_cache_dir}/{assembly}/*")]
        for release in releases:
            releasei=str2num(release) #FIXME is realease is a float
            genome_dir=f"{pyensembl_cache_dir}/{assembly}/{release}"
            genome_files=glob(f"{genome_dir}/*")
            is_genome_installed=True if len(genome_files)>4 else False #FIXME need more than 4 (.gz) files to be strict
            if is_genome_installed:
                dspecies.loc[genomei,'assembly']=assembly
                dspecies.loc[genomei,'release']=releasei
                dspecies.loc[genomei,'synonymn']=basename(genome_files[0]).split('.')[0]
                dspecies.loc[genomei,'latin name']=normalize_species_name(dspecies.loc[genomei,'synonymn'])
                genomei+=1
    # following loop generates the Species object
    for spc in dspecies['latin name'].unique():
        assembly2releases={}
        for assembly in dspecies.loc[(dspecies['latin name']==spc),'assembly'].unique():
            d=dspecies.loc[((dspecies['latin name']==spc) & (dspecies['assembly']==assembly)),:]
            assembly2releases[assembly]=d['release'].min(),d['release'].max() #FIXME if MAX_ENSEMBL_RELEASE very important and has to be used
        Species.register(
        latin_name=spc,
        synonyms=dspecies.loc[(dspecies['latin name']==spc),'synonymn'].unique().tolist(),
        reference_assemblies=assembly2releases)
        Species.dspecies=dspecies
    return Species
Species=collect_all_genomes()

# from .species import normalize_species_name,collect_all_genomes
# Species=collect_all_genomes()
# Species.register(
#         latin_name=normalize_species_name(args.species),
#         synonyms=[args.species],
#         reference_assemblies={
#         args.annotation_name: (args.reference_name, args.annotation_version),
#         })                
# for (species, release) in Species.all_species_release_pairs():
#     print(species, release)
# end of a PR (#207) from @rraadd88 