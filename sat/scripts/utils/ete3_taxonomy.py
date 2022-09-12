from ete3 import NCBITaxa
from collections import Counter


ncbi = NCBITaxa()


class Taxon:
    def __init__(self, taxonID, taxonomy_levels="", populate_class=True, is_name=False):
        """
        If populate_class is true, will add the following slots to the object:
        - lineage
        - canonical_lineage

        taxonomy_levels should be a list of taxonomy levels.
        """
        if not is_name:
            self.taxonID = taxonID
        else:
            taxonID = self.name_to_taxID(taxonID)
            self.taxonID = taxonID

        if populate_class:
            if taxonomy_levels == "":
                msg = "If you want to prepopulate the Taxon object, you need to provide"
                msg += " the taxonomy_levels upon instantiation!"
                raise ValueError(msg)
            self.lineage = self.get_lineage(self.taxonID)
            self.canonical_lineage = self.get_cannonical_lineage(
                self.taxonID, taxonomy_levels
            )

    def name_to_taxID(self, name: str):
        """
        Given name, returns taxonID. If the name taxonID can't be found,
        will try to strip any following numbers from the name. Otherwise,
        will return X.
        """
        try:
            return ncbi.get_name_translator([name])[name][0]
        except:
            # Try to repair by stripping any numbers from the end of the name
            name = name.rstrip("0123456789").rstrip(" ")
            try:
                return ncbi.get_name_translator([name])[name][0]
            except:
                return "X"

    def get_level(self, taxonID):
        """
        Given a single taxonID, returns the taxonomic level.
        """
        level = list(ncbi.get_rank([taxonID]).values())

        # Unknown taxonID would yield [], which can't be indexed by [0] to get the
        # string
        if level == []:
            level = "UNKNOWN"
        else:
            level = level[0]
        return level

    def get_name(self, taxonID):
        """
        Given a single taxonID, returns the name of the taxon.
        """
        name = list(ncbi.get_taxid_translator([taxonID]).values())

        # Unknown taxonID would yield [], which can't be indexed by [0] to get the
        # string
        if name == []:
            name = "UNKNOWN"
        else:
            name = name[0]

        name = name.replace(" ", "_")
        return name

    def get_lineage(self, taxonID):
        """
        Given a taxonID, returns a list of all taxonIDs in its lineage.
        """
        try:
            return ncbi.get_lineage(taxonID)
        except ValueError:
            print("Cannot find taxonID " + str(taxonID))
            return None

    def get_cannonical_lineage(
        self,
        taxonID,
        desired_levels: str = [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        terminal_as_species=True,
    ):
        """
        taxonID: interger taxonID

        desired_levels: list of desired levels - e.g. ['superkingdom', 'species']

        terminal_as_species:
            If species is the last level in desired levels but the taxonID is a level
            below species, will save the terminal/last level in place of the species.
            This is helpful in some cases where a viral taxonID is an unranked level
            more specific than species and multiple viruses share the same species
            taxonID.

        This function outputs the lineage as the NAMES of each taxon.
        """

        # Check input
        if type(desired_levels) != list:
            msg = "Desired levels should be a list! It should be ordered coherently."
            raise ValueError(msg)

        # Get a list of lineage and a list of their ranks
        lineage = self.get_lineage(taxonID)
        if lineage is None:
            return None
        levels = [self.get_level(taxonID) for taxonID in lineage]

        # Override species as terminal if specified. If species is the last level
        # requested but isn't the last level of the lineage, replace the species taxonID
        # with the terminal taxonID
        if terminal_as_species is True and "species" in levels:
            if desired_levels[-1] == "species":
                species_i = levels.index("species")
                if species_i + 1 < len(levels):
                    lineage[species_i] = lineage[-1]

        # Make lookup dict for easy conversion of level to taxonID
        taxonID_lookup = dict(zip(levels, lineage))

        cannonical_lineage = []
        # Iterate over each of the levels in prefix dictionary
        for level in desired_levels:

            if level not in desired_levels:
                continue

            # If the level isn't here, it is unknown
            if level not in taxonID_lookup:
                cannonical_lineage.append("")
                continue

            # Get the taxon name
            taxonID = taxonID_lookup[level]
            name = self.get_name(taxonID)

            # Report it out with the lineage
            cannonical_lineage.append(name)

        return cannonical_lineage


def taxonID_list_to_lineage_counts(taxonIDs, taxonomy_levels, observed_taxa=dict()):
    """
    Takes in a list or set of taxonIDs and returns a nested dictionary of structure:
    taxonomy_level --> taxon --> count.

    e.g.

    result['family']['polyomaviridae'] would yield 3 if there were 3 taxonIDs in the
    input that are from the family polyomaviridae.

    observed_taxa lets you input a dictionary of structure taxonID:Taxon() object to
    avoid searching the taxa database multiple times.
    """

    # Load taxon objects
    taxon_objects = []
    for taxonID in taxonIDs:

        if taxonID in observed_taxa:
            taxon_objects.append(observed_taxa[taxonID])
        else:
            taxon_object = Taxon(taxonID, taxonomy_levels)
            taxon_objects.append(taxon_object)
            observed_taxa[taxonID] = taxon_object

    # Count each taxon object at each level
    taxa_counts = dict()
    for i, level in enumerate(taxonomy_levels):

        for taxon_object in taxon_objects:

            lineage = taxon_object.canonical_lineage
            if lineage is None:
                continue

            taxa = lineage[i]

            if level not in taxa_counts:
                taxa_counts[level] = Counter()

            taxa_counts[level][taxa] += 1

    return taxa_counts, observed_taxa


if __name__ == "__main__":
    msg = "This script has utilities and functions. Don't call it directly!"
    raise ValueError(msg)
