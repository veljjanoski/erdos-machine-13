/* Local search for an (N-1)-state Erdos powers-of-two machine around an
   N-state one (erdos_14, so 13 states).

   A candidate is the base machine with one state X deleted (states above X
   renumbered), every rule that pointed to X rewritten ("forced" cells), and
   at most FMAX cells rewritten in total. Cells are fixed lazily: a cell gets
   its value the first time the machine reads it (forced cells: any of the
   4 (N-1) moves; other cells: keep, or any other move while the budget lasts).
   Machines that differ only in cells never read behave the same, so each leaf
   of the search stands for all of them. A halt is never chosen at a read: it
   would halt too early; the halt goes into an unread cell later.

   Filter (same tape layout as erdos_14): starting on the blank tape, for each
   n = 1..NS in order, the cells >= -2 must at some step be exactly the block
   {-1, ..., 2^n - 2} of ones (the stage start for 2^n), by step
   BUDGET_MULT * T14S[n] + 64. A match counts as the stage start only if the
   state is the one at stage 1 and the head is either at the same cell at
   every stage or at the same offset from the block end; other matches are
   ignored. At a stage start the cells < -2 (the counter) must be the same as
   at earlier stages with the same number c of exceptions so far (c = 0 for
   n <= 2, 1 for 3..8, 2 for 9) and differ from those with another c, or the
   machine is rejected. Halting rejects. Survivors are printed for an
   independent check in Python.

   Modes:
     e13search base   machine.txt              (the base machine)
     e13search mutant machine.txt              (all single-cell changes)
     e13search delete machine.txt FMAX X...    (the search)
     e13search free machine.txt X CELL...      (delete X; the named cells
                                                of the base machine, like D0,
                                                and the rules into X take any
                                                move; nothing else changes)
     e13search freesym machine.txt X LO HI CELL...
                                               (free, with the base states
                                                LO..HI, X among them, taken
                                                as interchangeable: a move may
                                                go to one already in use or
                                                to the next unused one, so
                                                each machine is met once)
     e13search estimate machine.txt PROBES X FMAX LO HI [CELL...]
                                               (Knuth's random-probe estimate
                                                of the search size: free mode
                                                if cells are given, else
                                                delete with FMAX; LO = -1:
                                                no symmetry)

   The search runs serially down to E13_SPLIT branchings (env, default 2),
   then the nodes there go to the threads; each node keeps only its table, and
   the worker replays it from the blank tape. Progress lines go to stderr.  */
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NS 9
#define LO 64
#define HI 1300
#define TAPE (LO + HI)
#define BUDGET_MULT 4
#define MAXQ 16
#define DEAD (-1)
#define SURV (-2)
#define FORCED (-2)

static uint64_t KEY[TAPE];
static uint64_t TGT[NS + 1];
static const long T14S[NS + 1] = {0, 2, 22, 90, 268, 944, 3434, 13450, 52068, 208864};
static const int CVAL[NS + 1] = {0, 0, 0, 1, 1, 1, 1, 1, 1, 2};
static int NQ, base[2 * MAXQ]; /* (next << 2) | (dir << 1) | write, dir 1 = R; -1 = halt */
static int SYM_LO = -1, SYM_HI = -1; /* interchangeable states, renumbered; -1 = off */
static int SYM_N0 = 0;               /* the first SYM_N0 of them are named by fixed cells */
static int SPLIT = 2;                /* branchings done serially before the threads start */

typedef struct {
  int head, st, nt, stage, st0, head0, abs_ok, rel_ok, lo, hi, cs[3];
  long t, deadline;
  uint64_t hb, hc, ch[3];
  unsigned char tape[TAPE];
} Sim;

typedef struct {
  long t; /* steps to reach the node; run() from the blank tape stops at the same cell */
  int tab[2 * MAXQ], budget, nintro;
  unsigned char fixed[2 * MAXQ];
} Node;

typedef struct {
  int X, nq, collecting;
  long long hist[NS + 1], leaves, survivors, steps;
  Node *front;
  long nfront, capfront;
  Sim buf[2 * MAXQ + 2], kids[4 * MAXQ + 1];
} Ctx;

static uint64_t splitmix(uint64_t *s) {
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static int parse_cell(const char *s) {
  if (s[0] == '-') return -1;
  return (atoi(s + 2) << 2) | ((s[1] == 'R') << 1) | (s[0] - '0');
}

static void load(const char *path) {
  FILE *f = fopen(path, "r");
  char a[16], b[16];
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
  NQ = 0;
  while (NQ < MAXQ && fscanf(f, "%15s %15s", a, b) == 2) {
    base[2 * NQ] = parse_cell(a);
    base[2 * NQ + 1] = parse_cell(b);
    NQ++;
  }
  fclose(f);
}

static void sim_init(Sim *s) {
  memset(s, 0, sizeof *s);
  s->st0 = -1;
  s->abs_ok = s->rel_ok = 1;
  s->nt = 1;
  s->deadline = BUDGET_MULT * T14S[1] + 64;
}

/* Runs until the machine reads an unfixed cell (returns its index), dies
   (DEAD) or reaches the start of stage NS (SURV). Re-entry at the same
   configuration is harmless: an accepted stage start moves the target on. */
static int run(Sim *s, const int *tab, const unsigned char *fixed) {
  for (;;) {
    if (s->hb == TGT[s->nt]) {
      int c = CVAL[s->nt], ok = 1, abs_ok = s->abs_ok, rel_ok = s->rel_ok;
      if (s->nt > 1) {
        if (s->st != s->st0) ok = 0;
        if (s->head != s->head0) abs_ok = 0;
        if (s->head - (1 << s->nt) != s->head0 - 2) rel_ok = 0;
        if (!abs_ok && !rel_ok) ok = 0;
      }
      if (ok) {
        if (s->cs[c] && s->ch[c] != s->hc) return DEAD;
        for (int j = 0; j < 3; j++)
          if (j != c && s->cs[j] && s->ch[j] == s->hc) return DEAD;
        s->cs[c] = 1;
        s->ch[c] = s->hc;
        if (s->nt == 1) {
          s->st0 = s->st;
          s->head0 = s->head;
        }
        s->abs_ok = abs_ok;
        s->rel_ok = rel_ok;
        s->stage = s->nt;
        if (s->nt == NS) return SURV;
        s->nt++;
        s->deadline = BUDGET_MULT * T14S[s->nt] + 64;
      }
    }
    if (s->t > s->deadline) return DEAD;
    int p = s->head + LO, sym = s->tape[p], i = 2 * s->st + sym;
    if (!fixed[i]) return i;
    int code = tab[i];
    if (code < 0) return DEAD;
    int w = code & 1;
    if (w != sym) {
      s->tape[p] = (unsigned char)w;
      if (s->head >= -2) s->hb ^= KEY[p]; else s->hc ^= KEY[p];
    }
    s->head += (code & 2) ? 1 : -1;
    if (s->head <= -LO || s->head >= HI - 1) return DEAD;
    s->st = code >> 2;
    s->t++;
  }
}

static void print_cell(FILE *o, int code) {
  if (code < 0) fprintf(o, "---");
  else fprintf(o, "%d%c%d", code & 1, (code & 2) ? 'R' : 'L', code >> 2);
}

static void report(Ctx *x, const Sim *s, const int *tab, const unsigned char *fixed) {
#pragma omp critical
  {
    printf("SURVIVOR X=%d st0=%d head0=%d abs=%d rel=%d steps=%ld table=", x->X, s->st0, s->head0,
           s->abs_ok, s->rel_ok, s->t);
    for (int q = 0; q < x->nq; q++) {
      for (int a = 0; a < 2; a++) {
        int i = 2 * q + a;
        if (!fixed[i]) fputc('?', stdout);
        print_cell(stdout, tab[i] == FORCED ? -1 : tab[i]);
        if (a == 0) fputc(' ', stdout);
      }
      fputc(q + 1 < x->nq ? '|' : '\n', stdout);
    }
    fflush(stdout);
  }
}

/* Symmetry breaking: returns the number of interchangeable states in use
   after a cell gets move `code`, or -1 if the move goes past the next unused
   one. Fixed cells never point into the set (checked in setup), so the states
   come into use in the order the machine first names them. */
static int sym_next(int code, int nintro) {
  if (SYM_LO < 0 || code < 0) return nintro;
  int q = code >> 2;
  if (q < SYM_LO || q > SYM_HI) return nintro;
  if (q - SYM_LO < nintro) return nintro;
  return q - SYM_LO == nintro ? nintro + 1 : -1;
}

static void explore(Ctx *x, int depth, int splits, int *tab, unsigned char *fixed, int budget, int nintro) {
  Sim *s = &x->buf[depth];
  long t0 = s->t;
  int i = run(s, tab, fixed);
  x->steps += s->t - t0;
  if (i == DEAD) { x->hist[s->stage]++; x->leaves++; return; }
  if (i == SURV) { x->survivors++; x->leaves++; report(x, s, tab, fixed); return; }
  int orig = tab[i], multi = orig == FORCED || budget > 0;
  if (x->collecting && splits == SPLIT) {
    if (x->nfront == x->capfront) {
      x->capfront = x->capfront ? 2 * x->capfront : 1024;
      x->front = (Node *)realloc(x->front, sizeof(Node) * (size_t)x->capfront);
      if (!x->front) { fprintf(stderr, "out of memory at front=%ld\n", x->nfront); exit(1); }
    }
    Node *nd = &x->front[x->nfront++];
    nd->t = s->t;
    memcpy(nd->tab, tab, sizeof nd->tab);
    memcpy(nd->fixed, fixed, sizeof nd->fixed);
    nd->budget = budget;
    nd->nintro = nintro;
    return;
  }
  Sim *c = &x->buf[depth + 1];
  fixed[i] = 1;
  if (orig != FORCED) {
    *c = *s;
    explore(x, depth + 1, splits + multi, tab, fixed, budget, nintro);
  }
  if (multi) {
    int nb = orig == FORCED ? budget : budget - 1;
    for (int v = 0; v < 4 * x->nq; v++) {
      int ni = sym_next(v, nintro);
      if (v == orig || ni < 0) continue;
      tab[i] = v;
      *c = *s;
      explore(x, depth + 1, splits + 1, tab, fixed, nb, ni);
    }
  }
  tab[i] = orig;
  fixed[i] = 0;
}

/* One random root-to-leaf path; adds weighted counts to out[0..3] =
   leaves, nodes, simulated steps, survivors (Knuth 1975: each average is an
   unbiased estimate of the full search's total). Every child of a node on
   the path is run, so children that die at once are counted exactly and the
   path only picks among the live ones; plain Knuth probes were far off on
   the finished searches because most of the tree sits under rare choices. */
static void probe(Ctx *x, const int *tab0, int budget, uint64_t *rng, double *out) {
  int tab[2 * MAXQ], nintro = SYM_N0;
  unsigned char fixed[2 * MAXQ] = {0};
  double w = 1;
  Sim *s = &x->buf[0];
  memcpy(tab, tab0, sizeof tab);
  sim_init(s);
  int i = run(s, tab, fixed);
  out[1] += 1;
  out[2] += (double)s->t;
  if (i < 0) {
    out[0] += 1;
    if (i == SURV) out[3] += 1;
    return;
  }
  for (;;) {
    int orig = tab[i], vals[4 * MAXQ + 1], nbs[4 * MAXQ + 1], nis[4 * MAXQ + 1], k = 0;
    if (orig != FORCED) { vals[k] = orig; nbs[k] = budget; nis[k] = nintro; k++; }
    if (orig == FORCED || budget > 0) {
      for (int v = 0; v < 4 * x->nq; v++) {
        int ni = sym_next(v, nintro);
        if (v == orig || ni < 0) continue;
        vals[k] = v; nbs[k] = orig == FORCED ? budget : budget - 1; nis[k] = ni; k++;
      }
    }
    int live[4 * MAXQ + 1], nlive = 0;
    fixed[i] = 1;
    for (int c = 0; c < k; c++) {
      Sim *ch = &x->kids[c];
      *ch = *s;
      tab[i] = vals[c];
      long t0 = ch->t;
      int r = run(ch, tab, fixed);
      out[1] += w;
      out[2] += w * (double)(ch->t - t0);
      if (r < 0) {
        out[0] += w;
        if (r == SURV) out[3] += w;
      } else {
        live[nlive++] = c;
      }
    }
    if (!nlive) return;
    int c = live[splitmix(rng) % (uint64_t)nlive];
    w *= nlive;
    tab[i] = vals[c];
    budget = nbs[c];
    nintro = nis[c];
    *s = x->kids[c];
    i = run(s, tab, fixed); /* stops at once at the same unfixed cell */
  }
}

/* Builds the table with X deleted; returns the budget for further changes,
   or -1 if the forced cells alone exceed FMAX. */
static int setup(int X, int fmax, const unsigned char *freebase, int *tab) {
  int nF = 0;
  for (int q = 0, s = 0; q < NQ; q++) {
    if (q == X) continue;
    for (int a = 0; a < 2; a++) {
      int code = base[2 * q + a], idx = 2 * s + a;
      if (freebase && freebase[2 * q + a]) {
        code = FORCED;
      } else if (code >= 0 && (code >> 2) == X) {
        code = FORCED;
        nF++;
      } else if (code >= 0) {
        int n = code >> 2;
        code = ((n - (X >= 0 && n > X)) << 2) | (code & 3); /* X < 0: no deletion */
      }
      tab[idx] = code;
    }
    s++;
  }
  if (SYM_LO >= 0) {
    /* fixed cells may name only the first few states of the set (these keep
       their labels); the others must have both cells free */
    unsigned char named[MAXQ] = {0};
    int n0 = 0;
    for (int i = 0; i < 2 * (NQ - 1); i++)
      if (tab[i] >= 0 && (tab[i] >> 2) >= SYM_LO && (tab[i] >> 2) <= SYM_HI) named[(tab[i] >> 2) - SYM_LO] = 1;
    while (n0 <= SYM_HI - SYM_LO && named[n0]) n0++;
    for (int q = SYM_LO + n0; q <= SYM_HI; q++)
      if (named[q - SYM_LO] || tab[2 * q] != FORCED || tab[2 * q + 1] != FORCED) {
        fprintf(stderr, "symmetric state %d is named by a fixed cell or has one\n", q);
        exit(1);
      }
    if (SYM_LO == 0) { fprintf(stderr, "start state in the symmetric states\n"); exit(1); }
    SYM_N0 = n0;
  }
  if (freebase) fmax = nF; /* only the free and forced cells change */
  if (nF > fmax) { fprintf(stderr, "X=%d forced=%d > fmax\n", X, nF); return -1; }
  return fmax - nF;
}

static void run_delete(int X, int fmax, const unsigned char *freebase) {
  int tab[2 * MAXQ], nq = NQ - (X >= 0);
  unsigned char fixed[2 * MAXQ] = {0};
  int budget = setup(X, fmax, freebase, tab);
  if (budget < 0) return;
  double w0 = omp_get_wtime();
  Ctx *top = (Ctx *)calloc(1, sizeof(Ctx));
  top->X = X;
  top->nq = nq;
  top->collecting = 1;
  sim_init(&top->buf[0]);
  explore(top, 0, 0, tab, fixed, budget, SYM_N0);
  long long hist[NS + 1], leaves = top->leaves, survivors = top->survivors, steps = top->steps;
  memcpy(hist, top->hist, sizeof hist);
  long nfront = top->nfront, done = 0;
  long long steps_done = steps;
  fprintf(stderr, "X=%d split=%d front=%ld (%.0f MB) after %.1fs\n", X, SPLIT, nfront,
          sizeof(Node) * (double)nfront / 1e6, omp_get_wtime() - w0);
  fflush(stderr);
#pragma omp parallel
  {
    Ctx *x = (Ctx *)calloc(1, sizeof(Ctx));
    x->X = X;
    x->nq = nq;
    int k;
#pragma omp for schedule(dynamic, 1)
    for (k = 0; k < (int)nfront; k++) {
      Node *nd = &top->front[k];
      int t2[2 * MAXQ];
      unsigned char f2[2 * MAXQ];
      long long s0 = x->steps;
      memcpy(t2, nd->tab, sizeof t2);
      memcpy(f2, nd->fixed, sizeof f2);
      sim_init(&x->buf[0]);
      explore(x, 0, SPLIT, t2, f2, nd->budget, nd->nintro);
      x->steps -= nd->t; /* the replay to the node was counted when it was collected */
#pragma omp critical(progress)
      {
        done++;
        steps_done += x->steps - s0;
        if (100LL * done / nfront != 100LL * (done - 1) / nfront) {
          fprintf(stderr, "  %lld%% of nodes (%ld/%ld) steps=%.3e wall=%.0fs\n", 100LL * done / nfront, done,
                  nfront, (double)steps_done, omp_get_wtime() - w0);
          fflush(stderr);
        }
      }
    }
#pragma omp critical
    {
      for (int j = 0; j <= NS; j++) hist[j] += x->hist[j];
      leaves += x->leaves;
      survivors += x->survivors;
      steps += x->steps;
    }
    free(x);
  }
  fprintf(stderr, "X=%d budget=%d split=%d front=%ld leaves=%lld survivors=%lld steps=%lld wall=%.1fs dead-at-stage(0..%d):",
          X, budget, SPLIT, nfront, leaves, survivors, steps, omp_get_wtime() - w0, NS);
  for (int j = 0; j <= NS; j++) fprintf(stderr, " %lld", hist[j]);
  fprintf(stderr, "\n");
  fflush(stderr);
  free(top->front);
  free(top);
}

/* LO..HI are base states (X among them), stored renumbered. */
static void set_sym(int X, int lo, int hi) {
  SYM_LO = lo - (lo > X);
  SYM_HI = hi - (hi >= X);
}

static int run_full(const int *tab, Sim *s) {
  unsigned char fixed[2 * MAXQ];
  memset(fixed, 1, sizeof fixed);
  sim_init(s);
  return run(s, tab, fixed);
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: e13search base|mutant|delete machine.txt [FMAX X...]\n"); return 1; }
  uint64_t seed = 20261028;
  const char *split_env = getenv("E13_SPLIT");
  if (split_env) SPLIT = atoi(split_env);
  for (int i = 0; i < TAPE; i++) KEY[i] = splitmix(&seed);
  for (int n = 1; n <= NS; n++) {
    TGT[n] = 0;
    for (int c = -1; c <= (1 << n) - 2; c++) TGT[n] ^= KEY[c + LO];
  }
  load(argv[2]);
  Sim *s = (Sim *)malloc(sizeof(Sim));
  if (!strcmp(argv[1], "base")) {
    int r = run_full(base, s);
    printf("base states=%d result=%s stage=%d st0=%d head0=%d abs=%d rel=%d steps=%ld\n", NQ,
           r == SURV ? "SURV" : "DEAD", s->stage, s->st0, s->head0, s->abs_ok, s->rel_ok, s->t);
  } else if (!strcmp(argv[1], "mutant")) {
    int tab[2 * MAXQ];
    for (int i = 0; i < 2 * NQ; i++)
      for (int v = -1; v < 4 * NQ; v++) {
        if (v == base[i]) continue;
        memcpy(tab, base, sizeof tab);
        tab[i] = v;
        int r = run_full(tab, s);
        printf("%d %d %d %d\n", i, v, s->stage, r == SURV);
      }
  } else if (!strcmp(argv[1], "delete")) {
    int fmax = atoi(argv[3]);
    for (int i = 4; i < argc; i++) run_delete(atoi(argv[i]), fmax, NULL);
  } else if (!strcmp(argv[1], "free") || !strcmp(argv[1], "freesym")) {
    unsigned char freebase[2 * MAXQ] = {0};
    int X = atoi(argv[3]), a0 = 4;
    if (!strcmp(argv[1], "freesym")) {
      set_sym(X, atoi(argv[4]), atoi(argv[5]));
      a0 = 6;
    }
    for (int i = a0; i < argc; i++) freebase[2 * (argv[i][0] - 'A') + (argv[i][1] - '0')] = 1;
    run_delete(X, 0, freebase);
  } else if (!strcmp(argv[1], "estimate")) {
    long probes = atol(argv[3]);
    int X = atoi(argv[4]), fmax = atoi(argv[5]), tab[2 * MAXQ];
    unsigned char freebase[2 * MAXQ] = {0};
    if (atoi(argv[6]) >= 0) set_sym(X, atoi(argv[6]), atoi(argv[7]));
    for (int i = 8; i < argc; i++) freebase[2 * (argv[i][0] - 'A') + (argv[i][1] - '0')] = 1;
    int budget = setup(X, fmax, argc > 8 ? freebase : NULL, tab);
    if (budget < 0) return 1;
    double sum[4] = {0}, sq[4] = {0}, w0 = omp_get_wtime();
#pragma omp parallel
    {
      Ctx *x = (Ctx *)calloc(1, sizeof(Ctx));
      double ls[4] = {0}, lq[4] = {0};
      /* splitmix streams that start a step apart repeat each other, so each
         thread starts at a hashed point of the sequence */
      const char *seed_env = getenv("E13_SEED");
      uint64_t rng = 0x5EED0000ULL + (uint64_t)omp_get_thread_num() +
                     ((seed_env ? strtoull(seed_env, NULL, 10) : 0) << 20);
      rng = splitmix(&rng);
      x->nq = NQ - (X >= 0);
      int p;
#pragma omp for schedule(dynamic, 256)
      for (p = 0; p < (int)probes; p++) {
        double o[4] = {0};
        probe(x, tab, budget, &rng, o);
        for (int j = 0; j < 4; j++) { ls[j] += o[j]; lq[j] += o[j] * o[j]; }
      }
#pragma omp critical
      for (int j = 0; j < 4; j++) { sum[j] += ls[j]; sq[j] += lq[j]; }
      free(x);
    }
    const char *name[4] = {"leaves", "nodes", "steps", "survivors"};
    printf("estimate X=%d budget=%d probes=%ld wall=%.1fs\n", X, budget, probes, omp_get_wtime() - w0);
    for (int j = 0; j < 4; j++) {
      double m = sum[j] / probes, var = sq[j] / probes - m * m;
      printf("  %-9s %.3e  (+- %.1e standard error)\n", name[j], m, var > 0 ? sqrt(var / probes) : 0.0);
    }
  }
  free(s);
  return 0;
}
