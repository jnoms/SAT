# from sat.scripts.aln_merge import aln_merge_main
# from sat.scripts.utils.Foldseek_Dataset import (
#     Foldseek_Dataset,
#     compare_foldseek_datasets,
# )

# test_files_dir = "tests/test_data/foldseek_related"


# def test_aln_merge_assure_expected_number_of_alignments():
#     aln1 = f"{test_files_dir}/add_tax_virus_vs_af2db.m8"
#     aln1_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore,taxid"
#     )
#     data1 = Foldseek_Dataset()
#     data1.parse_alignment(aln1, aln1_fields)
#     data1_initial_count = data1.count_alignments()

#     aln2 = f"{test_files_dir}/add_tax_virus_vs_virus.m8"
#     aln2_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore"
#     )
#     data2 = Foldseek_Dataset()
#     data2.parse_alignment(aln2, aln2_fields)

#     data1.merge(data2)

#     assert data1_initial_count + data2.count_alignments() == data1.count_alignments()


# def test_aln_merge(tmp_path):

#     # Args
#     class args:
#         pass

#     args.aln1 = f"{test_files_dir}/add_tax_virus_vs_af2db.m8"
#     args.aln1_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore,taxid"
#     )

#     args.aln2 = f"{test_files_dir}/add_tax_virus_vs_virus.m8"
#     args.aln2_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore"
#     )

#     args.aln1_source = ""
#     args.aln2_source = ""

#     args.output = f"{tmp_path}/merged.m8"

#     # Run script
#     aln_merge_main(args)

#     # Read in expected and observed datasets
#     observed = Foldseek_Dataset()
#     observed.parse_alignment(args.output)

#     expected = Foldseek_Dataset()
#     expected.parse_alignment(f"{test_files_dir}/merged.m8")

#     # Do checks
#     compare_foldseek_datasets(expected, observed)


# def test_aln_merge_source(tmp_path):

#     # Args
#     class args:
#         pass

#     args.aln1 = f"{test_files_dir}/add_tax_virus_vs_af2db.m8"
#     args.aln1_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore,taxid"
#     )

#     args.aln2 = f"{test_files_dir}/add_tax_virus_vs_virus.m8"
#     args.aln2_fields = (
#         "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,"
#         "tend,evalue,bits,alntmscore"
#     )

#     args.aln1_source = "af2db"
#     args.aln2_source = "virus"

#     args.output = f"{tmp_path}/merged_source.m8"

#     # Run script
#     aln_merge_main(args)

#     # Read in expected and observed datasets
#     observed = Foldseek_Dataset()
#     observed.parse_alignment(args.output)

#     expected = Foldseek_Dataset()
#     expected.parse_alignment(f"{test_files_dir}/merged_source.m8")

#     # Do checks
#     compare_foldseek_datasets(expected, observed)
