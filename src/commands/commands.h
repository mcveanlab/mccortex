#ifndef COMMANDS_H_
#define COMMANDS_H_

#include "cmd.h"
#include "cmd_mem.h"

int ctx_build(int argc, char **argv);
int ctx_infer_edges(int argc, char **argv);
int ctx_thread(int argc, char **argv);
int ctx_correct(int argc, char **argv);
int ctx_coverage(int argc, char **argv);
int ctx_rmsubstr(int argc, char **argv);
int ctx_breakpoints(int argc, char **argv);
int ctx_bubbles(int argc, char **argv);
int ctx_pview(int argc, char **argv);
int ctx_view(int argc, char **argv);
int ctx_clean(int argc, char **argv);
int ctx_pjoin(int argc, char **argv);
int ctx_supernodes(int argc, char **argv);
int ctx_health_check(int argc, char **argv);
int ctx_calls2vcf(int argc, char **argv);
// int ctx_geno(int argc, char **argv); // not written yet

int ctx_join(CmdArgs *args);
int ctx_subgraph(CmdArgs *args);
int ctx_reads(CmdArgs *args);
int ctx_contigs(CmdArgs *args);

int ctx_unique(CmdArgs *args); // retiring
int ctx_place(CmdArgs *args); // retiring
int ctx_extend(CmdArgs *args); // retiring

extern const char build_usage[];
extern const char view_usage[];
extern const char health_usage[];
extern const char clean_usage[];
extern const char join_usage[];
extern const char supernodes_usage[];
extern const char subgraph_usage[];
extern const char reads_usage[];
extern const char extend_usage[]; // retiring
extern const char contigs_usage[];
extern const char inferedges_usage[];
extern const char thread_usage[];
extern const char pview_usage[];
extern const char pjoin_usage[];
extern const char bubbles_usage[];
extern const char unique_usage[]; // retiring
extern const char place_usage[]; // retiring
extern const char correct_usage[];
extern const char breakpoints_usage[];
extern const char coverage_usage[];
extern const char rmsubstr_usage[];
extern const char calls2vcf_usage[];
// extern const char geno_usage[];

#endif /* COMMANDS_H_ */
