# from collections import Counter

# from .ete3_taxonomy import taxon_list_to_lineage_counts


# ------------------------------------------------------------------------------------ #
# Classes
# ------------------------------------------------------------------------------------ #
class Cluster_information:
    """
    Cluster_information class is a way to easily parse cluster information files from
    either foldseek or the similar files from generate_superclusters.

    The fields cluster_rep and cluster_member are REQUIRED.

    Fields can be specified as a comma-delimited string or a list in parse_cluster_file
    function. The fields can also be left blank, in which case the function will try to
    parse the headers. The first field must be cluster_rep if you want to use headers
    present in the file.s

    The organization of information on these objects is as follows:
    - the .clusters slot holds a dictionary of cluster_rep: Cluster() object.
    Each Cluster() object has the slot .cluster_members, which holds a list of
    Cluster_member() objects. Each Cluster_member() object holds a slot from all
    fields in the input cluster information file.

    Overview of storage structure:
    Cluster_info() object --> .clusters [dictionary --> cluster_rep:Cluster()]
        --> Cluster() object --> .cluster_members [list of Cluster_member() objects]
            --> Cluster_member() --> stores all cluster fields in slots.

    If you only need basic services such as converting a cluster representative to a
    set of its members, or finding the cluster representative for a member, you can
    use the following slots:
    - .cluster_rep_to_members: Dictionary of each rep to a set of their members
    - .cluster_member_to_rep: Dictionary of each member to its rep
    """

    def parse_cluster_file(self, cluster_file, cluster_file_fields):
        """
        cluster_file is the path to the cluster_file.
        cluster_file_fields is a comma-delimited string or a list containing the
        column names of the columns in the cluster_file. You can also pass an empty
        string if the fields are present as the first row of the cluster_file - in this
        case, cluster_rep must be the first column. The fields cluster_rep and
        cluster_member are MANDATORY. All other fields will be saved in the
        .clusters --> Cluster() --> Cluster_member() object slots.
        """

        # Format alignment fields from comma-delimited string to list if neceesary
        if cluster_file_fields != "" and not isinstance(cluster_file_fields, list):
            cluster_file_fields = cluster_file_fields.split(",")

        cluster_dict = dict()

        with open(cluster_file) as infile:

            # Check if the first line is a header - if so, and alignment_fields is
            # an empty string, use those as the alignment_fields
            if cluster_file_fields == "":
                for line in infile:
                    if line.startswith("cluster_rep"):
                        cluster_file_fields = line.rstrip("\n").split("\t")
                        break
                    else:
                        msg = "cluster_file_fields has not been passed to"
                        msg += " parse_cluster_file which is only allowed when the"
                        msg += " first line has headers! (e.g. first line should start"
                        msg += " with 'cluster_rep')"
                        raise ValueError(msg)

            # Record alignment fields for later use
            self.input_cluster_file_fields = cluster_file_fields

            # Make sure the required cluster fields are there
            for required_field in ["cluster_rep", "cluster_member"]:
                if required_field not in cluster_file_fields:
                    msg = f"Cannot find a required file field, {required_field},"
                    msg += " in the cluster file!"
                    raise ValueError(msg)

            # Get a set of all members and all cluster_reps - can be useful later
            all_reps = set()
            all_members = set()

            # Also make a cluster_rep_to_members dict - this is a simple conversion of the
            # cluster rep to the string representation of the cluster members. Also
            # make the reverse
            cluster_rep_to_members = dict()
            cluster_member_to_rep = dict()

            # Store cluster objects in a dictionary with cluster_rep as their key
            cluster_dict = dict()

            for line in infile:

                # It's okay if there is a header, but need to remove it
                if line.startswith("cluster_rep"):
                    continue

                line = line.rstrip("\n").split("\t")

                if len(line) != len(cluster_file_fields):
                    msg = "The line and cluster_file_fields don't have the same "
                    msg += "number of entries!"
                    msg += f"Current line is: {line}.\n"
                    msg += f"Alignment fields are: {cluster_file_fields}"
                    raise ValueError(msg)

                # Parse fields into a dictionary
                line_dict = dict()
                for (cluster_entry, cluster_field) in zip(line, cluster_file_fields):
                    line_dict[cluster_field] = cluster_entry

                # Initilize a cluster object if one is not present
                cluster_rep = line_dict["cluster_rep"]
                if cluster_rep not in cluster_dict:
                    cluster = Cluster(cluster_rep)
                    cluster_dict[cluster_rep] = cluster

                # Add information to the cluster object
                # Here, will be storing the line in a Cluster_member object. Then,
                # will be storing all of the cluster members within the cluster object's
                # cluster_members slot
                cluster_member = Cluster_member(line_dict)
                cluster_dict[cluster_rep].cluster_members.append(cluster_member)

                # Also store all cluster members and cluster reps into a set - may be
                # useful later
                all_reps.add(cluster_rep)
                all_members.add(cluster_member.cluster_member)

                # Update a simple lookup dict to convert rep-->members and
                # member-->rep
                if cluster_rep not in cluster_rep_to_members:
                    cluster_rep_to_members[cluster_rep] = set()
                cluster_rep_to_members[cluster_rep].add(cluster_member.cluster_member)
                if cluster_member.cluster_member not in cluster_member_to_rep:
                    cluster_member_to_rep[cluster_member.cluster_member] = cluster_rep
                else:
                    msg = "Have found that a cluster_member seems to have two reps!"
                    msg += f" Member: {cluster_member.cluster_member},"
                    msg += " Current Rep in the dictionary:"
                    msg += f" {cluster_member_to_rep[cluster_member.cluster_member]}"
                    msg += f" Current rep in loop: {cluster_rep}"
                    raise ValueError(msg)

        # Store the dictionary of cluster_rep --> cluster objects (which in turn have
        # cluster member objects) as the .clusters slot
        self.clusters = cluster_dict
        self.all_reps = all_reps
        self.all_members = all_members
        self.cluster_rep_to_members = cluster_rep_to_members
        self.cluster_member_to_rep = cluster_member_to_rep


class Cluster:
    """
    The purpose of this class is to store a lsit of cluster_member() objects. This
    class is initialized with the string name of the cluster_rep. The cluster_member()
    objects can be stored in the .cluster_members slot, which is a list.
    """

    def __init__(self, cluster_rep):
        self.cluster_rep = cluster_rep
        self.cluster_members = list()

    def __str__(self):
        return self.cluster_rep

    def __repr__(self):
        return self.cluster_rep


class Cluster_member:
    """
    The purpose of this object is to store the lines from the cluster file. So, each
    field in each line gets saved as a slot here.
    """

    def __init__(self, cluster_member_dict):

        for required_field in ["cluster_rep", "cluster_member"]:
            if required_field not in cluster_member_dict:
                msg = f"Cannot find a required file field, {required_field},"
                msg += " in the cluster file! Current fields are "
                msg += f"{cluster_member_dict.keys()}"
                raise ValueError(msg)

        for key, val in cluster_member_dict.items():
            self.__dict__[key] = val

    def __str__(self):
        return self.cluster_member

    def __repr__(self):
        return self.cluster_member


# ------------------------------------------------------------------------------------ #
# Functions
# ------------------------------------------------------------------------------------ #
