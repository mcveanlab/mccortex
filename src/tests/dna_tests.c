#include "global.h"
#include "all_tests.h"
#include <ctype.h>

#include "dna.h"

void test_dna_functions()
{
  test_status("[dna] Testing all dna.h functions...");

  char c;
  size_t i;

  for(c = CHAR_MIN; c < CHAR_MAX; c++) {
    TASSERT((c != 0 && strchr("ACGTacgt",c) != NULL) == char_is_acgt(c));
    TASSERT((c != 0 && strchr("ACGTNacgtn",c) != NULL) == char_is_acgtn(c));
    if(char_is_acgt(c))
      TASSERT(dna_char_complement(dna_char_complement(c)) == c);
  }

  // Loop over nucs
  for(i = 0; i < 4; i++) {
    TASSERT(dna_char_complement("TGCA"[i]) == "ACGT"[i]);
    TASSERT(dna_char_complement("tgca"[i]) == "acgt"[i]);
    TASSERT(dna_nuc_to_char(i) == "ACGT"[i]);
    TASSERT(dna_char_to_nuc(dna_nuc_to_char(i)) == i);
    TASSERT(dna_nuc_complement(dna_nuc_complement(i)) == i);
  }

  // Loop over possible char bases
  for(i = 0; i < 8; i++) {
    c = "ACGTacgt"[i];
    TASSERT2(dna_nuc_to_char(dna_char_to_nuc(c)) == toupper(c), "char: %c", c);
  }

  //
  // Test dna_reverse_complement_str
  //
  const char fwd[] = "TtCCagACgAGGCcGCaGCAGCTCTTATGC";
  const char rev[] = "GCATAAGAGCTGCtGCgGCCTcGTctGGaA";
  char str[100];
  size_t len = strlen(fwd);
  strcpy(str, fwd);

  // revcmp 0bp, remains the same
  dna_reverse_complement_str(str,0);
  TASSERT2(strcmp(str,fwd) == 0, "str: %s vs fwd: %s", str, fwd);

  // revcmp 1bp, then revert
  dna_reverse_complement_str(str,1);
  TASSERT(strcmp(str+1,fwd+1) == 0);
  TASSERT(str[0] == dna_char_complement(fwd[0]));
  dna_reverse_complement_str(str,1);
  TASSERT(strcmp(str,fwd) == 0);

  // revcmp 10bp, then revert
  dna_reverse_complement_str(str,10);
  TASSERT(strcmp(str+10,fwd+10) == 0);
  for(i = 0; i < 10; i++)
    TASSERT(str[i] == dna_char_complement(fwd[10-1-i]));
  dna_reverse_complement_str(str,10);
  TASSERT(strcmp(str,fwd) == 0);

  // revcmp the whole string
  dna_reverse_complement_str(str,len);
  TASSERT(strcmp(str,rev) == 0);
}
