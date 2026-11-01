/* recheck.c -- independent re-check of e13 survivor lines (review of the V1 result, 2026-10-31).
   Written without looking at deep.c's loop. For each "SURVIVOR ... table=..." line: run from the
   blank tape; a stage start for n is state st0, head at head0 (abs) or head0 + 2^n - 2 (rel), and
   the ones at cells >= -2 exactly {-1, ..., 2^n - 2}. At each stage start record the EXACT set of
   ones at cells <= -3 (the counter). A correct machine has the same counter at every stage with
   the same c(n) (c = 0 for n <= 2, 1 for 3..8, 2 for n >= 9: 2^9 .. 2^16 all have a ternary
   digit 2) and different counters for different c. Cells printed as "?---" are free cells the
   search never read: reading one is reported, because the machine's behaviour from there on is
   not determined. Cells printed "?move" are non-free cells never read; their only allowed value
   is the printed one.  Usage: recheck NMAX < survivors.txt   (MSVC: long is 32-bit, so long long) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
#define L 4000000LL
#define R 4000000LL
#define T(c) tp[(c) + L]
static unsigned char *tp;
static char ctr[20][8192];

int main(int argc, char **argv) {
  int nmax = argc > 1 ? atoi(argv[1]) : 15;
  tp = calloc((size_t)(L + R + 1), 1);
  char line[8192];
  int id = 0, nbad = 0, ngood = 0;
  while (fgets(line, sizeof line, stdin)) {
    char *tb = strstr(line, "table=");
    if (!tb) continue;
    id++;
    int st0 = 0, head0 = -2, abs_ok = 1, rel_ok = 0;
    char *p;
    if ((p = strstr(line, "st0="))) st0 = atoi(p + 4);
    if ((p = strstr(line, "head0="))) head0 = atoi(p + 6);
    if ((p = strstr(line, "abs="))) abs_ok = atoi(p + 4);
    if ((p = strstr(line, "rel="))) rel_ok = atoi(p + 4);
    int rel = !abs_ok && rel_ok;
    int wr[32][2], mv[32][2], nx[32][2], kind[32][2], nq = 0; /* kind: 0 move, 1 halt, 2 free */
    p = tb + 6;
    while (*p && *p != '\n' && *p != '\r') {
      for (int a = 0; a < 2; a++) {
        while (*p == ' ') p++;
        int unread = 0;
        if (*p == '?') { unread = 1; p++; }
        if (*p == '-') { kind[nq][a] = unread ? 2 : 1; p += 3; }
        else {
          kind[nq][a] = 0;
          wr[nq][a] = *p++ - '0';
          mv[nq][a] = (*p++ == 'R') ? 1 : -1;
          nx[nq][a] = (int)strtol(p, &p, 10);
        }
      }
      nq++;
      while (*p == ' ') p++;
      if (*p == '|') p++; else break;
    }
    memset(tp, 0, (size_t)(L + R + 1));
    ll head = 0, t = 0, minhead = 0, ones = 0;
    int st = 0, failed = 0;
    char verdict[512] = "";
    for (int n = 1; n <= nmax && !failed; n++) {
      ll blk = 1LL << n;
      ll limit = t + 8LL * (1LL << (2 * n)) + 100000;
      ll want = rel ? head0 + blk - 2 : head0;
      for (;;) {
        if (st == st0 && head == want && ones == blk) {
          int ok = 1;
          for (ll k = -1; k < blk - 1; k++) if (!T(k)) { ok = 0; break; }
          if (ok) break;
        }
        if (t >= limit) { sprintf(verdict, "stage %d not reached (limit)", n); failed = 1; break; }
        if (head <= -L + 2 || head >= R - 2) { sprintf(verdict, "stage %d: off tape", n); failed = 1; break; }
        int s = T(head), k = kind[st][s];
        if (k == 2) { sprintf(verdict, "free cell %c%d read before stage %d", 'A' + st, s, n); failed = 1; break; }
        if (k == 1) { sprintf(verdict, "halt before stage %d", n); failed = 1; break; }
        int w = wr[st][s];
        if (w != s) { T(head) = (unsigned char)w; if (head >= -2) ones += w ? 1 : -1; }
        head += mv[st][s];
        if (head < minhead) minhead = head;
        st = nx[st][s];
        t++;
      }
      if (failed) break;
      char *q = ctr[n]; *q = 0;
      int len = 0;
      for (ll k = -3; k >= minhead - 1 && len < 8000; k--)
        if (T(k)) len += sprintf(q + len, "%lld,", k);
      int c = n <= 2 ? 0 : (n <= 8 ? 1 : 2);
      int first = c == 0 ? 1 : (c == 1 ? 3 : 9);
      if (n != first && strcmp(ctr[n], ctr[first]))
        { sprintf(verdict, "counter changes at n=%d (%s -> %s)", n, ctr[first][0] ? ctr[first] : "none", ctr[n][0] ? ctr[n] : "none"); failed = 1; }
    }
    if (!failed) {
      if (nmax >= 9 && (!strcmp(ctr[1], ctr[3]) || !strcmp(ctr[3], ctr[9]) || !strcmp(ctr[1], ctr[9])))
        sprintf(verdict, "counter equal for different c");
      else sprintf(verdict, "CONSISTENT to n=%d (c2 counter %s)", nmax, ctr[9]);
    }
    if (strncmp(verdict, "CONSISTENT", 10) == 0) ngood++; else nbad++;
    printf("%d: %s\n", id, verdict);
    fflush(stdout);
  }
  fprintf(stderr, "machines %d consistent %d failed %d\n", id, ngood, nbad);
  return 0;
}
