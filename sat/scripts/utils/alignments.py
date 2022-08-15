from .structure import pdb_to_structure_object, structure_to_pLDDT, structure_to_seq
from .ete3_taxonomy import get_cannonical_lineage

from collections import Counter
import numpy as np
import matplotlib as plt


# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Alignment_group:
    """
    Holds all alignments and other information for a specific query.

    Holds a list of Alignment_objects, the query structure,
    the query pLDDT, etc. Can also calculate and plot the coverage of all alignments on
    the query.
    """

    def __init__(self, query):
        self.query = query
        self.alignments = []

    def add_alignment(self, alignment):
        """
        Makes sure the query and target aren't the same.
        """
        if alignment.target == alignment.query:
            pass
        else:
            self.alignments.append(alignment)

    def add_structure(self, structure_path):
        """
        Loads structure as a biopython structure object.
        Loads the pLDDT at each position as a dict of structure
        position:pLDDT. Position is 1-indexed, as are
        the alignments.
        """
        self.structure = pdb_to_structure_object(structure_path)
        self.pLDDTs = structure_to_pLDDT(self.structure)
        self.seq = structure_to_seq(self.structure)
        self.query_len = len(self.seq)

    def add_scaled_alignment_coverage(self):
        """
        Counts the number of alignments at every position of the query. Returns
        a np array containing the coverage at each position of the query. The
        array is one-dimensional... the first value of the array is the
        coverage at the first residue of the protein (e.g. 1-indexed), and so
        on.

        Rather than an alignment contributing +1 for each position covered, it
        will contribute + (1 * 1/percent_AA_identity)... this will prioritize
        alignments at low sequence identity.
        """

        # initialize counter at every position
        c = Counter()
        for i in range(self.query_len):
            i += 1  # adjusting to 1 indexing
            c[i] = 0

        # Iterate through each alignment and add
        for alignment in self.alignments:

            # Working with 1-indexed positions, so want end to be "inclusive"
            for i in range(alignment.qstart, alignment.qend + 1):
                c[i] += 1 * 1 / float(alignment.fident)

        # Return array of coverages at each position. Array length is equal to query_len
        cov_array = np.array(list(c.values()))
        if len(cov_array) != self.query_len:
            msg = "The length of the coverage array isn't equal to query len for some "
            msg += "reason! Something is wrong."
            raise ValueError(msg)

        self.scaled_alignment_coverage = cov_array

    def add_alignment_coverage(self):
        """
        Counts the number of alignments at every position of the query. Returns
        a np array containing the coverage at each position of the query. The
        array is one-dimensional... the first value of the array is the
        coverage at the first residue of the protein (e.g. 1-indexed), and so
        on.
        """

        # initialize counter at every position
        c = Counter()
        for i in range(self.query_len):
            i += 1  # adjusting to 1 indexing
            c[i] = 0

        # Iterate through each alignment and add
        for alignment in self.alignments:

            # Working with 1-indexed positions, so want end to be "inclusive"
            for i in range(alignment.qstart, alignment.qend + 1):
                c[i] += 1

        # Return array of coverages at each position. Array length is equal to query_len
        cov_array = np.array(list(c.values()))
        if len(cov_array) != self.query_len:
            msg = "The length of the coverage array isn't equal to query_len for some "
            msg += "reason! Something is wrong."
            raise ValueError(msg)

        self.alignment_coverage = cov_array

    def plot_coverage(self):

        if (
            "alignment_coverage" not in self.__dict__
            and "scaled_alignment_coverage" not in self.__dict__
        ):
            msg = "Cannot find alignment_coverage or scaled_alignment_coverage "
            msg += "in the object!"
            raise ValueError(msg)

        fig, ax = plt.subplots()
        if "alignment_coverage" in self.__dict__:
            ax.plot(self.alignment_coverage, label="Alignment Coverage")
        if "scaled_alignment_coverage" in self.__dict__:
            ax.plot(
                self.scaled_alignment_coverage
                * np.max(self.alignment_coverage)
                / np.max(self.scaled_alignment_coverage),
                label="Scaled Alignment Coverage",
            )
        ax.legend(bbox_to_anchor=(1.0, 1.0))
        ax.set_title(self.query)
        ax.set_ylabel("Alignment Coverage")
        ax.set_xlabel("Query AA Position")
        return fig


class Alignment_object:
    """
    Holds information for a single alignment. Note that the alignment
    coordinates are 1-indexed.

    alignment_field and alignment_entry should be lists of the same size, where
    the alignment_entry at index i corresponds to the alignment_field at index i.
    """

    def __init__(self, alignment_entries: list, alignment_fields: list):
        if len(alignment_fields) != len(alignment_entries):
            msg = "alignment_fields and alignment_entries must be the same length!"
            raise ValueError(msg)

        for (alignment_field, alignment_entry) in zip(
            alignment_fields, alignment_entries
        ):
            self.__dict__[alignment_field] = alignment_entry

        self.alignment_fields = alignment_fields

    def trim_by_plddt(self, pLDDTs: dict, threshold: int):
        """
        Trims the alignment to remove query sequences at the ends if they have
        a pLDDT below threshold.

        pLDDTs: Dict of structure position:pLDDT of the query
        Threshold: int value that marks the minimum pLDDT to avoid trimming
        at the ends.
        """
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)

        qstart_pLDDT = pLDDTs[self.qstart]
        qend_pLDDT = pLDDTs[self.qstart]

        while qstart_pLDDT < threshold:
            # If none of the alignments meet the pLDDT threshold, qstart == qend
            # In this case, mark qstart and qend as None
            if self.qstart == self.qend:
                self.qstart = None
                self.qend = None
                return

            self.qstart += 1
            qstart_pLDDT = pLDDTs[self.qstart]

        while qend_pLDDT < threshold:
            self.qend -= 1
            qend_pLDDT = pLDDTs[self.qend]

    def write_output(self, fields: list, sep: str = "\t"):
        """
        Returns the desired fields as a tab-delimited string followed by a \\n.
        """
        out = []
        for field in fields:
            val = self.__dict__[field]
            if type(val) == list:
                val = [str(v) for v in val]
                val = sep.join(val)
            elif type(val) in {float, int}:
                val = str(val)
            out.append(val)
        return sep.join(out) + "\n"

    def add_taxid(self, query_or_target, location):
        if query_or_target not in {"query", "target"}:
            msg = "Must specify if you want the taxonID of the query or of the target"
            raise ValueError(msg)

        field_name = f"{query_or_target}_taxid"

        if location == 0:
            self.__dict__[field_name] = ""
            return
        elif location == 1:
            taxid = self.__dict__[query_or_target].rstrip(".pdb").split("__")[-1]
            if "domain" in taxid:
                while "_" in taxid:
                    taxid = taxid[:-1]
            self.__dict__[field_name] = taxid
        elif location == 2:
            if query_or_target == "query":
                msg = "For query, location option 2 is not allowed."
                raise ValueError(msg)
            taxid = self.taxid
            if "domain" in taxid:
                while "_" in taxid:
                    taxid = taxid[:-1]
            self.__dict__[field_name] = taxid

    def add_query_lineage(
        self,
        taxonomy_levels: list = [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    ):
        """
        Returns the taxonomy of the query. The query taxid must be built into the query
        name. - e.g. query_taxid_location must be set to 1.
        """
        taxid = self.query_taxid
        if taxid == "":
            self.query_lineage = ""
        else:
            self.query_lineage = get_cannonical_lineage(taxid, taxonomy_levels)

    def add_target_lineage(
        self,
        taxonomy_levels: list = [
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    ):
        """
        Returns the taxonomy of the target. The target taxid can either be present
        in the target name, or in the taxid field.
        """
        taxid = self.target_taxid
        if taxid == "":
            self.target_lineage = ""
        else:
            self.target_lineage = get_cannonical_lineage(taxid, taxonomy_levels)


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
def parse_alignment(alignment_file_path, alignment_fields):
    """
    Given an input path, parses the alignment file into a dictionary of
    structure query:Alingment_group, where each Alignment_group contains
    a list of Alignment objects.
    """
    alignments_dict = dict()
    with open(alignment_file_path) as infile:

        for line in infile:
            line = line.rstrip("\n").split("\t")

            if len(line) != len(alignment_fields):
                msg = "The line and alignment_fields don't have the same "
                msg += "number of entries!"
                raise ValueError(msg)

            alignment = Alignment_object(line, alignment_fields)

            # Save to output dictionary
            if alignment.query not in alignments_dict:
                alignments_dict[alignment.query] = Alignment_group(alignment.query)
                alignments_dict[alignment.query].add_alignment(alignment)
            else:
                alignments_dict[alignment.query].add_alignment(alignment)

    return alignments_dict
