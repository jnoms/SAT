from .structure import pdb_to_structure_object, structure_to_pLDDT, structure_to_seq
from .ete3_taxonomy import Taxon

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt


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
            for i in range(int(alignment.qstart), int(alignment.qend) + 1):
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
            for i in range(int(alignment.qstart), int(alignment.qend) + 1):
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

    def filter_alignments(self, filter_field, max_val, min_val):
        """
        This function iterates through each alignment in this alignment_object and
        only keeps alignments that have a value in the slot filter_field that is
        <= a max_val and >= a min_val.

        This function acts IN PLACE.
        """

        if max_val < min_val:
            msg = (
                f"max_val can't be higher than min_val! max_val: {max_val}, "
                f"min_val:{min_val}"
            )
            raise ValueError(msg)

        filtered_alignments = []
        for alignment in self.alignments:
            try:
                val = alignment.__dict__[filter_field]
            except KeyError:
                msg = f"Cannot find the field {filter_field} in the alingments."
                msg += " Something is wrong!"
                raise ValueError(msg)

            # Handle scientific notation and convert to float
            # Update the field as well
            if isinstance(val, str):
                try:
                    val = float(val.replace("E", "e"))
                except ValueError:
                    msg = (
                        f"Detected a string in the field {filter_field} which can't be "
                        f"converted to a string. The val is {val}."
                    )
                    raise ValueError(msg)

            # Override the object with the formatted field
            alignment.__dict__[filter_field] = val

            if val <= max_val and val >= min_val:
                filtered_alignments.append(alignment)

        self.alignments = filtered_alignments

    def keep_top_N_alignments(self, filter_field, N):
        """
        This function sorts the list of alingments (in self.alignments) such that the
        alignment with a higher value in the filter_field slot is first, and so on. It
        then keeps the first N alignments, tossing out the others. This acts in place!

        If N is set to 0, will return all alingments and won't filter.
        """
        if self.alignments == []:
            return

        if N == 0:
            return

        if filter_field not in self.alignments[0].__dict__:
            msg = f"Cannot find the field {filter_field} in the alingments."
            msg += " Something is wrong!"
            raise ValueError(msg)

        self.alignments.sort(key=lambda x: x.__dict__[filter_field], reverse=True)
        self.alignments = self.alignments[:N]


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

    def __eq__(self, other):
        """
        Lets you compare two alignment objects to one another. Helpful for tests.

        To use: one_alignment_obj.__eq__(another_alignment_obj) will return True
        if the slots contain the same values.
        """
        return self.__dict__ == other.__dict__

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

        If the field is a list, it will delimit each item in the list with the sep.
        """
        out = []
        for field in fields:
            val = self.__dict__[field]

            # If the slot is a list, collapse it by the delimiter.
            if isinstance(val, list):
                val = [str(v) for v in val]
                val = sep.join(val)
            elif isinstance(val, float) or isinstance(val, int):
                val = str(val)

            # If the slot is a Taxon object, output the collapsed canonical_lineage.
            elif isinstance(val, Taxon):
                val = val.canonical_lineage
                val = sep.join(val)
            out.append(val)
        return sep.join(out) + "\n"

    def add_query_taxonID(self, delimiter="__", pos=-1):
        """
        This will search the query column string for the taxonID and add it to the
        query_taxonID slot of the alignment_object
        """
        query_taxonID = self.query.split(delimiter)[pos].rstrip(".pdb")
        if "domain" in query_taxonID:
            while "_" in query_taxonID:
                query_taxonID = query_taxonID[:-1]
        self.query_taxonID = query_taxonID

    def add_target_taxonID(self, delimiter="__", pos=-1, field_priority=True):
        """
        Adds the target taxonID. Here, if field_priority is set to true, it will first
        check if there is a taxid slot which was from the foldseek column... this column
        would have the target taxonID. If field_priority is false, will ignore that
        and take it from the string of the target_name
        """
        if field_priority and hasattr(self, "taxid"):
            self.target_taxonID = self.taxid
        else:
            target_taxonID = self.target.split(delimiter)[pos].rstrip(".pdb")
            if "domain" in target_taxonID:
                while "_" in target_taxonID:
                    target_taxonID = target_taxonID[:-1]
            self.target_taxonID = target_taxonID


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
