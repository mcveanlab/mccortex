#ifndef COMMANDS_H_
#define COMMANDS_H_

#include "cmd.h"

int ctx_build(int argc, char **argv);
int ctx_view(CmdArgs *args);
int ctx_health_check(CmdArgs *args);
int ctx_clean(CmdArgs *args);
int ctx_join(CmdArgs *args);
int ctx_supernodes(CmdArgs *args);
int ctx_subgraph(CmdArgs *args);
int ctx_reads(CmdArgs *args);
int ctx_extend(CmdArgs *args);
int ctx_contigs(CmdArgs *args);
int ctx_infer_edges(int argc, char **argv);
int ctx_thread(int argc, char **argv);
int ctx_pview(CmdArgs *args);
int ctx_pjoin(CmdArgs *args);
int ctx_bubbles(CmdArgs *args);
int ctx_diverge(CmdArgs *args);
int ctx_unique(CmdArgs *args);
int ctx_place(CmdArgs *args);
int ctx_correct(int argc, char **argv);
int ctx_breakpoints(CmdArgs *args);
int ctx_coverage(CmdArgs *args);
int ctx_rmsubstr(CmdArgs *args);

extern const char build_usage[];
extern const char view_usage[];
extern const char health_usage[];
extern const char clean_usage[];
extern const char join_usage[];
extern const char supernodes_usage[];
extern const char subgraph_usage[];
extern const char reads_usage[];
extern const char extend_usage[];
extern const char contigs_usage[];
extern const char inferedges_usage[];
extern const char thread_usage[];
extern const char pview_usage[];
extern const char pjoin_usage[];
extern const char bubbles_usage[];
extern const char diverge_usage[];
extern const char unique_usage[];
extern const char covg_usage[];
extern const char place_usage[];
extern const char correct_usage[];
extern const char breakpoints_usage[];
extern const char coverage_usage[];
extern const char rmsubstr_usage[];

#endif /* COMMANDS_H_ */
