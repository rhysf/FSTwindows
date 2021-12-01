# FSTwindows
Calculate FST values across sliding windows

# Example

Name_tab_Lineage.tab contains:

isolate1	ST1
isolate2	ST1
isolate3	ST2
isolate4	ST2
isolate5	ST2
isolate6	ST3

etc,

To compare ST1 to ST2:

perl VCF-multi-sample-to-FST-windows.pl -a Name_tab_Lineage.tab -v Multi_sample.VCF -c ST1 -d ST2 -r reference.fasta > out.window-10kb.tab
