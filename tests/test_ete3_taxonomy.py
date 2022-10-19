from sat.scripts.utils.ete3_taxonomy import taxon_list_to_lineage_counts, Taxon

from collections import Counter

import pytest


@pytest.mark.ete3
def test_taxon_list_to_lineage_counts_multi():
    taxonIDs = {"2025360", "2200830", "1010"}
    taxonomy_levels = "superkingdom,phylum,class,order,family,genus,species".split(",")
    taxons = []
    for taxonID in taxonIDs:
        taxon = Taxon(taxonID, taxonomy_levels)
        taxons.append(taxon)

    expected = {
        "superkingdom": Counter({"Viruses": 2, "Bacteria": 1}),
        "phylum": Counter({"Nucleocytoviricota": 2, "Bacteroidetes": 1}),
        "class": Counter({"Pokkesviricetes": 2, "Sphingobacteriia": 1}),
        "order": Counter({"Chitovirales": 2, "Sphingobacteriales": 1}),
        "family": Counter({"Poxviridae": 2, "Sphingobacteriaceae": 1}),
        "genus": Counter(
            {"Centapoxvirus": 1, "Sphingobacterium": 1, "Orthopoxvirus": 1}
        ),
        "species": Counter(
            {"NY_014_poxvirus": 1, "Sphingobacterium_mizutaii": 1, "Akhmeta_virus": 1}
        ),
    }

    observed, _ = taxon_list_to_lineage_counts(taxons, taxonomy_levels)
    for level, c in observed.items():
        for taxon, count in c.items():
            assert expected[level][taxon] == count


@pytest.mark.ete3
def test_taxonID_list_to_lineage_counts_empty():
    taxonIDs = {}
    taxonomy_levels = "superkingdom,phylum,class,order,family,genus,species".split(",")
    taxons = []
    for taxonID in taxonIDs:
        taxon = Taxon(taxonID, taxonomy_levels)
        taxons.append(taxon)

    expected = {}

    observed, _ = taxon_list_to_lineage_counts(taxons, taxonomy_levels)
    for level, c in observed.items():
        for taxon, count in c.items():
            assert expected[level][taxon] == count


@pytest.mark.ete3
def test_taxonID_list_to_lineage_counts_unknown_taxa():
    taxonIDs = {"2025360", "2200830", "1010", "324123412421412421"}
    taxonomy_levels = "superkingdom,phylum,class,order,family,genus,species".split(",")
    taxons = []
    for taxonID in taxonIDs:
        taxon = Taxon(taxonID, taxonomy_levels)
        taxons.append(taxon)

    expected = {
        "superkingdom": Counter({"Viruses": 2, "Bacteria": 1, "": 1}),
        "phylum": Counter({"Nucleocytoviricota": 2, "Bacteroidetes": 1, "": 1}),
        "class": Counter({"Pokkesviricetes": 2, "Sphingobacteriia": 1, "": 1}),
        "order": Counter({"Chitovirales": 2, "Sphingobacteriales": 1, "": 1}),
        "family": Counter({"Poxviridae": 2, "Sphingobacteriaceae": 1, "": 1}),
        "genus": Counter(
            {"Centapoxvirus": 1, "Sphingobacterium": 1, "Orthopoxvirus": 1, "": 1}
        ),
        "species": Counter(
            {
                "NY_014_poxvirus": 1,
                "Sphingobacterium_mizutaii": 1,
                "Akhmeta_virus": 1,
                "": 1,
            }
        ),
    }

    observed, _ = taxon_list_to_lineage_counts(taxons, taxonomy_levels)
    for level, c in observed.items():
        for taxon, count in c.items():
            assert expected[level][taxon] == count
