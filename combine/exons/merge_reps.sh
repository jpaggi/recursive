BEGIN="python sequence/included_exons/merge_reps.py"
MOTIFS="no_motif.bed"


$BEGIN rs_exons/5_repA/$MOTIFS rs_exons/5_repB/$MOTIFS rs_exons/5_repC/$MOTIFS > rs_exons/5_merged_$MOTIFS
$BEGIN rs_exons/10_repA/$MOTIFS rs_exons/10_repB/$MOTIFS rs_exons/10_repC/$MOTIFS > rs_exons/10_merged_$MOTIFS
$BEGIN rs_exons/20_repA/$MOTIFS rs_exons/20_repB/$MOTIFS rs_exons/20_repC/$MOTIFS > rs_exons/20_merged_$MOTIFS

$BEGIN rs_exons/total_repA/$MOTIFS rs_exons/total_repB/$MOTIFS > rs_exons/total_merged_$MOTIFS


python sequence/included_exons/merge.py rs_exons/5_merged_$MOTIFS rs_exons/10_merged_$MOTIFS rs_exons/20_merged_$MOTIFS rs_exons/total_merged_$MOTIFS > rs_exons/all_merged_$MOTIFS