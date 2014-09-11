#ifndef COMMANDS_H_
#define COMMANDS_H_

#include "cmd.h"
#include "cmd_mem.h"

int ctx_build(int argc, char **argv);
int ctx_sort(int argc, char **argv);
int ctx_infer_edges(int argc, char **argv);
int ctx_thread(int argc, char **argv);
int ctx_correct(int argc, char **argv);
int ctx_coverage(int argc, char **argv);
int ctx_rmsubstr(int argc, char **argv);
int ctx_breakpoints(int argc, char **argv);
int ctx_bubbles(int argc, char **argv);
int ctx_view(int argc, char **argv);
int ctx_clean(int argc, char **argv);
int ctx_pjoin(int argc, char **argv);
int ctx_contigs(int argc, char **argv);
int ctx_supernodes(int argc, char **argv);
int ctx_health_check(int argc, char **argv);
int ctx_calls2vcf(int argc, char **argv); // unfinished
int ctx_subgraph(int argc, char **argv);
int ctx_join(int argc, char **argv);
int ctx_reads(int argc, char **argv);

// int ctx_geno(int argc, char **argv); // not written yet

int ctx_unique(int argc, char **argv); // retiring
int ctx_place(int argc, char **argv); // retiring

// Experiments
int ctx_exp_abc(int argc, char **argv);

extern const char build_usage[];
extern const char sort_usage[];
extern const char view_usage[];
extern const char health_usage[];
extern const char clean_usage[];
extern const char join_usage[];
extern const char supernodes_usage[];
extern const char subgraph_usage[];
extern const char reads_usage[];
extern const char contigs_usage[];
extern const char inferedges_usage[];
extern const char thread_usage[];
extern const char pjoin_usage[];
extern const char bubbles_usage[];
extern const char correct_usage[];
extern const char breakpoints_usage[];
extern const char coverage_usage[];
extern const char rmsubstr_usage[];
extern const char calls2vcf_usage[];
// extern const char geno_usage[];

extern const char unique_usage[]; // retiring
extern const char place_usage[]; // retiring

// Experiments
extern const char exp_abc_usage[];

#endif /* COMMANDS_H_ */
