#include <Rdefines.h>


/* cigar_utils.c */

SEXP valid_cigar(SEXP cigar, SEXP ans_type);

SEXP explode_cigar_ops(
	SEXP cigar,
	SEXP ops
);

SEXP explode_cigar_op_lengths(
	SEXP cigar,
	SEXP ops
);

SEXP split_cigar(SEXP cigar);

SEXP cigar_op_table(SEXP cigar);

SEXP cigar_ranges(
	SEXP cigar,
	SEXP flag,
	SEXP space,
	SEXP pos,
	SEXP f,
	SEXP ops,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges,
	SEXP with_ops
);

SEXP cigar_width(
	SEXP cigar,
	SEXP flag,
	SEXP space
);

SEXP cigar_narrow(SEXP cigar, SEXP left_width, SEXP right_width);

SEXP cigar_qnarrow(SEXP cigar, SEXP left_qwidth, SEXP right_qwidth);

SEXP ref_locs_to_query_locs(
        SEXP ref_locs,
        SEXP cigar,
        SEXP pos,
        SEXP narrow_left
);

SEXP query_locs_to_ref_locs(
        SEXP query_locs,
        SEXP cigar,
        SEXP pos,
        SEXP narrow_left
);


/* transcript_utils.c */

SEXP transcript_widths(
	SEXP exonStarts,
	SEXP exonEnds
);

SEXP tlocs2rlocs(
	SEXP tlocs,
	SEXP exonStarts,
	SEXP exonEnds,
	SEXP strand,
	SEXP decreasing_rank_on_minus_strand
);

SEXP extract_transcripts(
	SEXP classname,
	SEXP x,
	SEXP exonStarts,
	SEXP exonEnds,
	SEXP strand,
	SEXP decreasing_rank_on_minus_strand,
	SEXP lkup
);

