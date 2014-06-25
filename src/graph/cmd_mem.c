#include "global.h"
#include "cmd_mem.h"
#include "util.h"
#include "hash_mem.h" // for calculating mem usage

#include "misc/mem_size.h" // in libs/misc/

void cmd_mem_args_set_memory(struct MemArgs *mem, const char *arg)
{
  if(mem->mem_to_use_set)
    cmd_print_usage("-m, --memory <M> specifed more than once");
  mem->mem_to_use = cmd_parse_arg_mem("-m, --memory <M>", arg);
  if(mem->mem_to_use == 0)
    cmd_print_usage("--memory <M> cannot be zero");
  mem->mem_to_use_set = true;
}

void cmd_mem_args_set_nkmers(struct MemArgs *mem, const char *arg)
{
  if(mem->num_kmers_set)
    cmd_print_usage("-n, --nkmers <N> specifed more than once");
  mem->num_kmers = cmd_parse_arg_mem("-n, --nkmers <M>", arg);
  if(mem->num_kmers == 0)
    cmd_print_usage("--nkmer <N> cannot be zero");
  mem->num_kmers_set = true;
}

void cmd_print_mem(size_t mem_bytes, const char *name)
{
  char mem_str[100];
  bytes_to_str(mem_bytes, 1, mem_str);
  status("[memory] %s: %s", name, mem_str);
}

// If your command accepts -n <kmers> and -m <mem> this may be useful
//  `extra_bits` is additional memory per node, above hash table+BinaryKmers
//  `use_mem_limit` if true, fill args->mem_to_use
size_t cmd_get_kmers_in_hash2(size_t mem_to_use, bool mem_to_use_set,
                              size_t num_kmers, bool num_kmers_set,
                              size_t extra_bits,
                              size_t min_num_kmer_req, size_t max_num_kmers_req,
                              bool use_mem_limit, size_t *graph_mem_ptr)
{
  size_t kmers_in_hash, graph_mem = 0, min_kmers_mem;
  char graph_mem_str[100], mem_to_use_str[100];
  char kmers_in_hash_str[100], min_num_kmers_str[100], min_kmers_mem_str[100];

  if(!use_mem_limit && min_num_kmer_req == 0 && !num_kmers_set) {
    cmd_print_usage("Cannot read from stream without -n <nkmers> set");
  }

  if(num_kmers_set)
    graph_mem = hash_table_mem(num_kmers, extra_bits, &kmers_in_hash);
  else if(use_mem_limit)
    graph_mem = hash_table_mem_limit(mem_to_use, extra_bits, &kmers_in_hash);
  else if(min_num_kmer_req > 0)
    graph_mem = hash_table_mem((size_t)(min_num_kmer_req/IDEAL_OCCUPANCY),
                               extra_bits, &kmers_in_hash);

  if(max_num_kmers_req > 0 && !num_kmers_set)
  {
    // Check if the max kmer capacity is less that requested
    size_t graph_mem2, kmers_in_hash2;
    graph_mem2 = hash_table_mem(max_num_kmers_req/IDEAL_OCCUPANCY,
                                extra_bits, &kmers_in_hash2);
    if(graph_mem2 < graph_mem) {
      graph_mem = graph_mem2;
      kmers_in_hash = kmers_in_hash2;
    }
  }

  if(kmers_in_hash < 1024)
    graph_mem = hash_table_mem(1024, extra_bits, &kmers_in_hash);
  // ^ 1024 is a very small default hash table capacity

  min_kmers_mem = hash_table_mem(min_num_kmer_req, extra_bits, NULL);

  bytes_to_str(graph_mem, 1, graph_mem_str);
  bytes_to_str(mem_to_use, 1, mem_to_use_str);
  bytes_to_str(min_kmers_mem, 1, min_kmers_mem_str);

  ulong_to_str(kmers_in_hash, kmers_in_hash_str);
  ulong_to_str(min_num_kmer_req, min_num_kmers_str);

  // Give a error/warning about occupancy
  if(kmers_in_hash < min_num_kmer_req)
  {
    die("Not enough kmers in hash: require at least %s kmers (min memory: %s)",
        min_num_kmers_str, min_kmers_mem_str);
  }
  else if(kmers_in_hash < min_num_kmer_req/WARN_OCCUPANCY)
  {
    warn("Expected hash table occupancy %.2f%% "
         "(you may want to increase -n or -m)",
         (100.0 * min_num_kmer_req) / kmers_in_hash);
  }

  if(mem_to_use_set && num_kmers_set) {
    if(num_kmers > kmers_in_hash) {
      die("-n <kmers> requires more memory than given with -m <mem> [%s > %s]",
          graph_mem_str, mem_to_use_str);
    }
    else if(use_mem_limit && graph_mem < mem_to_use) {
      status("Note: Using less memory than requested (%s < %s); allows for %s kmers",
             graph_mem_str, mem_to_use_str, kmers_in_hash_str);
    }
  }

  if(graph_mem > mem_to_use) {
    die("Not enough memory for requested graph: require at least %s [>%s]",
        graph_mem_str, mem_to_use_str);
  }

  cmd_print_mem(graph_mem, "graph");

  if(graph_mem_ptr != NULL) *graph_mem_ptr = graph_mem;

  return kmers_in_hash;
}

// Check memory against memory limit and machine memory
void cmd_check_mem_limit(size_t mem_to_use, size_t mem_requested)
{
  char memstr[100], ramstr[100];
  bytes_to_str(mem_requested, 1, memstr);

  if(mem_requested > mem_to_use)
    die("Need to set higher memory limit [ at least -m %s ]", memstr);

  // Get memory
  size_t ram = getMemorySize();
  bytes_to_str(ram, 1, ramstr);

  if(mem_requested > ram) {
    die("Requesting more memory than is available [ Reqeusted: -m %s RAM: %s ]",
        memstr, ramstr);
  }

  status("[memory] total: %s of %s RAM\n", memstr, ramstr);
}
