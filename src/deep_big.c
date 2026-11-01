/* deep.c with a 4M-cell tape on each side (for machines that leave deep.c's 4096 cells);
   at the step limit it prints head, state and the extent of the ones.
   Deep stage check for e13 survivors: run from the blank tape, find stage
   starts (state st0, head at head0 or at head0 - 2 + 2^n (abs/rel in the
   survivor line), ones at cells >= -2 exactly -1 .. 2^n - 2),
   print the number of counter ones (left of -2) and the leftmost one at each
   start, stop at the first
   stage not reached within the step limit or at a halt.
   Usage: deep NMAX < survivors.txt   (lines with "st0=S ... table=...") */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LEFT (1 << 22)
static unsigned char *tape;

int main(int argc, char **argv) {
  int nmax = argc > 1 ? atoi(argv[1]) : 16;
  long long size = LEFT + (1LL << nmax) + (1 << 22) + 4096;
  tape = malloc(size);
  char line[4096];
  int id = 0;
  while (fgets(line, sizeof line, stdin)) {
    char *tp = strstr(line, "table="), *sp = strstr(line, "st0=");
    if (!tp || !sp) continue;
    id++;
    int st0 = atoi(sp + 4);
    /* stage starts: absolute (head at head0) or relative to the block end (head - 2^n ==
       head0 - 2), as the search reports with abs= / rel= */
    char *hp = strstr(line, "head0="), *ap = strstr(line, "abs="), *rp = strstr(line, "rel=");
    long long head0 = hp ? atoi(hp + 6) : -2;
    int rel = ap && rp && !atoi(ap + 4) && atoi(rp + 4);
    /* halt[q][a]: 1 = halt cell, 2 = free cell the search never read ('?---'), which may be
       any move, so reading it is reported apart from a real halt */
    int w[32][2], d[32][2], nx[32][2], halt[32][2], nq = 0;
    char *p = tp + 6;
    for (;;) {                          /* rows separated by '|', cells by ' ' */
      for (int a = 0; a < 2; a++) {
        while (*p == ' ') p++;
        int unread = *p == '?';
        if (unread) p++;
        if (p[0] == '-') { halt[nq][a] = unread ? 2 : 1; p += 3; }
        else { halt[nq][a] = 0; w[nq][a] = p[0] - '0'; d[nq][a] = p[1] == 'R' ? 1 : -1;
               nx[nq][a] = strtol(p + 2, &p, 10); }
      }
      nq++;
      while (*p == ' ') p++;
      if (*p != '|') break;
      p++;
    }
    memset(tape, 0, size);
    long long head = 0, t = 0;
    int st = 0;
    long long ones_ge = 0;                   /* ones at cells >= -2 */
    printf("survivor %d (nq=%d st0=%d%s):", id, nq, st0, rel ? " rel" : "");
    int failed = 0;
    for (int n = 1; n <= nmax && !failed; n++) {
      long long blk = 1LL << n, limit = t + 64LL * (1LL << (2 * n)) + 100000;
      for (;;) {
        if (st == st0 && head == (rel ? blk + head0 - 2 : head0) && ones_ge == blk) {
          int ok = 1;
          for (long long k = -1; k < blk - 1; k++) if (!tape[LEFT + k]) { ok = 0; break; }
          if (ok) break;
        }
        if (t >= limit) { { long long lo1 = 0, hi1 = 0, ones = 0; for (long long k = -LEFT + 1; k < size - LEFT; k++) if (tape[LEFT + k]) { if (!ones) lo1 = k; hi1 = k; ones++; } printf(" | stage %d NOT reached (limit, t=%lld) head=%lld state=%c ones=%lld in [%lld,%lld]", n, t, head, (char)(65 + st), ones, lo1, hi1); } failed = 1; break; }
        if (head < -LEFT + 2 || head > blk + (1 << 22)) { printf(" | stage %d: off tape at %lld", n, head); failed = 1; break; }
        int s = tape[LEFT + head];
        if (halt[st][s] == 2) { printf(" | FREE cell %c%d read before stage %d (t=%lld)", 'A' + st, s, n, t); failed = 1; break; }
        if (halt[st][s]) { printf(" | HALT before stage %d (t=%lld)", n, t); failed = 1; break; }
        int ww = w[st][s];
        if (ww != s) { tape[LEFT + head] = ww; if (head >= -2) ones_ge += ww ? 1 : -1; }
        head += d[st][s];
        st = nx[st][s];
        t++;
      }
      if (failed) break;
      int cnt = 0; long long lo = 0;
      for (long long k = -3; k >= -LEFT + 1; k--) if (tape[LEFT + k]) { cnt++; lo = k; }
      printf(" n%d:%d", n, cnt);
      if (cnt) printf("@%lld", lo);
      fflush(stdout);
    }
    printf("\n");
  }
  return 0;
}
