"""GPU version of `e13search free/freesym` (CuPy RawKernel, NVRTC).

Same search as the C program: the base machine with state X deleted, the
named cells (and the rules into X) free, the rest fixed; a free cell gets its
move the first time the machine reads it; same filter as run() in
e13search.c (Zobrist hashes, stage targets, counter rule, deadlines) and the
same symmetry rule. The counts (leaves, steps, deaths per stage, survivors)
must equal the C program's exactly.

Tree: a live node is a machine stopped at the first read of an unassigned
free cell (or at the step cap, to keep each GPU call short). One GPU thread
takes one (node, move) pair, copies the node, sets the cell and runs on. Live
children go to an output buffer and then onto a depth-first stack in GPU
memory.

Usage: python e13gpu.py machine.txt X LO HI CELL...   (LO = -1: no symmetry)
Survivors go to stdout in the C format, the summary and progress to stderr.

Checkpoints: with env E13_CKPT=file, the stack, counters and survivors are
saved every E13_CKPT_SEC seconds (default 1800), between GPU calls, so the
saved state is exact; a run started with an existing file resumes from it.
Clean stop: create the file <E13_CKPT>.stop; within a second the run saves a
checkpoint, removes the stop file and exits with code 3 (nothing lost).
"""
import os
import sys
import time

import cupy as cp
import numpy as np

NS, LO, HI, BUDGET_MULT = 9, 64, 1300, 4
TAPE = LO + HI
TW = (TAPE + 31) // 32
T14S = [0, 2, 22, 90, 268, 944, 3434, 13450, 52068, 208864]
CVAL = [0, 0, 0, 1, 1, 1, 1, 1, 1, 2]
FORCED, HALT_CODE, FREE_CODE = -2, 62, 63
M64 = (1 << 64) - 1


def splitmix(state):
    state = (state + 0x9E3779B97F4A7C15) & M64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return state, z ^ (z >> 31)


def parse_cell(s):
    if s[0] == '-':
        return -1
    return (int(s[2:]) << 2) | ((s[1] == 'R') << 1) | int(s[0])


def load(path):
    base = []
    for line in open(path):
        a = line.split()
        if len(a) >= 2:
            base += [parse_cell(a[0]), parse_cell(a[1])]
    return base


def setup(base, X, freebase, lo, hi):
    """Mirror of set_sym() + setup() in e13search.c (free mode, budget 0)."""
    nq0 = len(base) // 2
    tab = []
    for q in range(nq0):
        if q == X:
            continue
        for a in (0, 1):
            code = base[2 * q + a]
            if freebase[2 * q + a]:
                code = FORCED
            elif code >= 0 and (code >> 2) == X:
                code = FORCED
            elif code >= 0:
                n = code >> 2
                code = ((n - (n > X)) << 2) | (code & 3)
            tab.append(code)
    sym_lo = sym_hi = -1
    n0 = 0
    if lo >= 0:
        sym_lo = lo - (lo > X)
        sym_hi = hi - (hi >= X)
        named = [0] * (sym_hi - sym_lo + 1)
        for c in tab:
            if c >= 0 and sym_lo <= (c >> 2) <= sym_hi:
                named[(c >> 2) - sym_lo] = 1
        while n0 <= sym_hi - sym_lo and named[n0]:
            n0 += 1
        for q in range(sym_lo + n0, sym_hi + 1):
            if named[q - sym_lo] or tab[2 * q] != FORCED or tab[2 * q + 1] != FORCED:
                sys.exit(f'symmetric state {q} is named by a fixed cell or has one')
        if sym_lo == 0:
            sys.exit('start state in the symmetric states')
    return tab, sym_lo, sym_hi, n0


KERNEL = r'''
typedef unsigned long long u64;
#define NQ %(NQ)d
#define NV (4 * NQ)
#define NS_ %(NS)d
#define LO_ %(LO)d
#define HI_ %(HI)d
#define TAPE_ %(TAPE)d
#define TW_ %(TW)d
#define SYM_LO %(SYM_LO)d
#define SYM_HI %(SYM_HI)d
#define BUDGET_MULT_ %(BUDGET_MULT)d
__constant__ int T14S_[NS_ + 1] = {%(T14S)s};
__constant__ int CVAL_[NS_ + 1] = {%(CVAL)s};

/* flags: bit 0 abs_ok, bit 1 rel_ok, bits 2..4 cs[0..2].
   tab: 6-bit codes, 10 cells per word; 62 = halt, 63 = unassigned free cell. */
struct Node {
  u64 hb, hc, ch0, ch1, ch2;
  int head, st, nt, stage, st0, head0;
  int flags, t, deadline, nintro, pend, mask;
  u64 tab[3];
  unsigned int tape[TW_];
};

extern "C" __global__ void node_size(int *out) { out[0] = sizeof(Node); }

__device__ __forceinline__ int get_code(u64 a, u64 b, u64 c, int i) {
  u64 w = i < 10 ? a : (i < 20 ? b : c);
  int j = i < 10 ? i : (i < 20 ? i - 10 : i - 20);
  return (int)((w >> (6 * j)) & 63);
}

__device__ __forceinline__ void set_code(u64 &a, u64 &b, u64 &c, int i, int v) {
  int j = i < 10 ? i : (i < 20 ? i - 10 : i - 20);
  u64 m = ~(63ULL << (6 * j)), x = (u64)v << (6 * j);
  if (i < 10) a = (a & m) | x;
  else if (i < 20) b = (b & m) | x;
  else c = (c & m) | x;
}

__device__ __forceinline__ int sym_next(int code, int nintro) {
  if (SYM_LO < 0) return nintro;
  int q = code >> 2;
  if (q < SYM_LO || q > SYM_HI) return nintro;
  if (q - SYM_LO < nintro) return nintro;
  return q - SYM_LO == nintro ? nintro + 1 : -1;
}

/* ctri: [0] next work item, [1] children written, [2] survivors written.
   ctr: [0] leaves, [1] steps, [2] survivors, [3 + s] deaths at stage s. */
extern "C" __global__ void expand(const Node *__restrict__ in, int nin, Node *__restrict__ out,
                                  int *ctri, u64 *ctr, int *surv, int survcap,
                                  const u64 *__restrict__ keyg, const u64 *__restrict__ tgtg, int cap,
                                  int slice) {
  __shared__ u64 key[TAPE_];
  __shared__ u64 tgt_sh[NS_ + 1];
  extern __shared__ unsigned int tsh[];
  for (int i = threadIdx.x; i < TAPE_; i += blockDim.x) key[i] = keyg[i];
  if (threadIdx.x <= NS_) tgt_sh[threadIdx.x] = tgtg[threadIdx.x];
  __syncthreads();
  const int BD = blockDim.x;
  unsigned int *tp = tsh + threadIdx.x; /* word w of this thread's tape at tp[w * BD] */
  u64 l_leaves = 0, l_steps = 0, l_surv = 0;
  unsigned int l_hist[NS_ + 1];
  for (int j = 0; j <= NS_; j++) l_hist[j] = 0;
  const int nitems = nin * NV;
  /* this thread's current machine; it runs `slice` steps per round, so one long run does
     not hold up the other threads of its warp (they wait only for the round) */
  u64 hb = 0, hc = 0, ch0 = 0, ch1 = 0, ch2 = 0, ta = 0, tb = 0, tc = 0, tgt = 0;
  int head = 0, st = 0, nt = 1, stage = 0, st0 = 0, head0 = 0, flags = 0, t = 0, deadline = 0, mask = 0;
  int nintro = 0, t0 = 0, seg = 0, active = 0;
  for (;;) {
    if (!active) {
      int k = atomicAdd(&ctri[0], 1);
      if (k >= nitems) break;
      int ni = k / NV, v = k - ni * NV;
      const Node *nd = in + ni;
      int pend = nd->pend;
      nintro = nd->nintro;
      if (pend < 0) {
        if (v) continue; /* a capped node just goes on */
      } else {
        int n2 = sym_next(v, nintro);
        if (n2 < 0) continue;
        nintro = n2;
      }
      hb = nd->hb; hc = nd->hc; ch0 = nd->ch0; ch1 = nd->ch1; ch2 = nd->ch2;
      head = nd->head; st = nd->st; nt = nd->nt; stage = nd->stage; st0 = nd->st0; head0 = nd->head0;
      flags = nd->flags; t = nd->t; deadline = nd->deadline; mask = nd->mask;
      ta = nd->tab[0]; tb = nd->tab[1]; tc = nd->tab[2];
      for (int w = 0; w < TW_; w++) tp[w * BD] = nd->tape[w];
      if (pend >= 0) set_code(ta, tb, tc, pend, v);
      t0 = t;
      seg = 0;
      tgt = tgt_sh[nt];
      active = 1;
    }
    int res = 4, pi = -1; /* res: 0 dead, 1 survivor, 2 free-cell read, 3 step cap, 4 round over */
    for (int r = 0; r < slice; r++) {
      if (seg >= cap) { res = 3; break; } /* checked first: resuming re-enters here */
      if (hb == tgt) {
        int c = CVAL_[nt], ok = 1, abs_ok = flags & 1, rel_ok = (flags >> 1) & 1;
        if (nt > 1) {
          if (st != st0) ok = 0;
          if (head != head0) abs_ok = 0;
          if (head - (1 << nt) != head0 - 2) rel_ok = 0;
          if (!abs_ok && !rel_ok) ok = 0;
        }
        if (ok) {
          u64 chc = c == 0 ? ch0 : (c == 1 ? ch1 : ch2);
          if (((flags >> (2 + c)) & 1) && chc != hc) { res = 0; break; }
          if ((c != 0 && (flags & 4) && ch0 == hc) || (c != 1 && (flags & 8) && ch1 == hc) ||
              (c != 2 && (flags & 16) && ch2 == hc)) { res = 0; break; }
          flags |= 4 << c;
          if (c == 0) ch0 = hc; else if (c == 1) ch1 = hc; else ch2 = hc;
          if (nt == 1) { st0 = st; head0 = head; }
          flags = (flags & ~3) | abs_ok | (rel_ok << 1);
          stage = nt;
          if (nt == NS_) { res = 1; break; }
          nt++;
          deadline = BUDGET_MULT_ * T14S_[nt] + 64;
          tgt = tgt_sh[nt];
        }
      }
      if (t > deadline) { res = 0; break; }
      int p = head + LO_;
      unsigned int word = tp[(p >> 5) * BD];
      int sym = (word >> (p & 31)) & 1, i = 2 * st + sym;
      mask |= 1 << i;
      int code = get_code(ta, tb, tc, i);
      if (code == 63) { res = 2; pi = i; break; }
      if (code == 62) { res = 0; break; }
      if ((code & 1) != sym) {
        tp[(p >> 5) * BD] = word ^ (1u << (p & 31));
        if (head >= -2) hb ^= key[p]; else hc ^= key[p];
      }
      head += (code & 2) ? 1 : -1;
      if (head <= -LO_ || head >= HI_ - 1) { res = 0; break; }
      st = code >> 2;
      t++;
      seg++;
    }
    if (res == 4) continue; /* same machine, next round */
    active = 0;
    l_steps += (u64)(t - t0);
    if (res == 0) {
      l_leaves++;
      l_hist[stage]++;
    } else if (res == 1) {
      l_leaves++;
      l_surv++;
      int s = atomicAdd(&ctri[2], 1);
      if (s < survcap) {
        int *b = surv + 32 * s;
        for (int j = 0; j < 2 * NQ; j++) b[j] = get_code(ta, tb, tc, j);
        b[26] = mask; b[27] = st0; b[28] = head0; b[29] = flags & 1; b[30] = (flags >> 1) & 1; b[31] = t;
      }
    } else {
      Node *o = out + atomicAdd(&ctri[1], 1);
      o->hb = hb; o->hc = hc; o->ch0 = ch0; o->ch1 = ch1; o->ch2 = ch2;
      o->head = head; o->st = st; o->nt = nt; o->stage = stage; o->st0 = st0; o->head0 = head0;
      o->flags = flags; o->t = t; o->deadline = deadline; o->nintro = nintro;
      o->pend = res == 2 ? pi : -1; o->mask = mask;
      o->tab[0] = ta; o->tab[1] = tb; o->tab[2] = tc;
      for (int w = 0; w < TW_; w++) o->tape[w] = tp[w * BD];
    }
  }
  atomicAdd(&ctr[0], l_leaves);
  atomicAdd(&ctr[1], l_steps);
  atomicAdd(&ctr[2], l_surv);
  for (int j = 0; j <= NS_; j++)
    if (l_hist[j]) atomicAdd(&ctr[3 + j], (u64)l_hist[j]);
}

/* the root: blank tape, as sim_init() */
extern "C" __global__ void init_root(Node *o, u64 ta, u64 tb, u64 tc, int nintro) {
  o->hb = o->hc = o->ch0 = o->ch1 = o->ch2 = 0;
  o->head = 0; o->st = 0; o->nt = 1; o->stage = 0; o->st0 = -1; o->head0 = 0;
  o->flags = 3; o->t = 0; o->deadline = BUDGET_MULT_ * T14S_[1] + 64; o->nintro = nintro;
  o->pend = -1; o->mask = 0;
  o->tab[0] = ta; o->tab[1] = tb; o->tab[2] = tc;
  for (int w = 0; w < TW_; w++) o->tape[w] = 0;
}
'''


def cell_str(code):
    if code < 0:
        return '---'
    return f"{code & 1}{'R' if code & 2 else 'L'}{code >> 2}"


def main():
    t_start = time.time()
    path, X, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    base = load(path)
    freebase = [0] * len(base)
    for c in sys.argv[5:]:
        freebase[2 * (ord(c[0]) - ord('A')) + int(c[1])] = 1
    tab, sym_lo, sym_hi, n0 = setup(base, X, freebase, lo, hi)
    nq = len(tab) // 2
    assert nq <= 15 and 4 * nq <= 61  # 26 cells fit 3 words; codes below 62

    seed, key = 20261028, []
    for _ in range(TAPE):
        seed, k = splitmix(seed)
        key.append(k)
    tgt = [0] * (NS + 1)
    for n in range(1, NS + 1):
        for c in range(-1, (1 << n) - 1):
            tgt[n] ^= key[c + LO]

    words = [0, 0, 0]
    for i, c in enumerate(tab):
        code = FREE_CODE if c == FORCED else (HALT_CODE if c < 0 else c)
        words[i // 10] |= code << (6 * (i % 10))

    src = KERNEL % dict(NQ=nq, NS=NS, LO=LO, HI=HI, TAPE=TAPE, TW=TW, SYM_LO=sym_lo, SYM_HI=sym_hi,
                        BUDGET_MULT=BUDGET_MULT, T14S=', '.join(map(str, T14S)),
                        CVAL=', '.join(map(str, CVAL)))
    mod = cp.RawModule(code=src, options=('-std=c++14',))
    k_size, k_expand, k_root = (mod.get_function(n) for n in ('node_size', 'expand', 'init_root'))
    tmp = cp.zeros(1, dtype=cp.int32)
    k_size((1,), (1,), (tmp,))
    NB = int(tmp[0])
    nv = 4 * nq

    stack_cap = 40_000_000  # nodes (11.5 GB; V3 peaked at 3.6M)
    bmax = 1 << 17
    stack = cp.empty(stack_cap * NB, dtype=cp.uint8)
    out = cp.empty(bmax * nv * NB, dtype=cp.uint8)
    survcap = 200_000
    surv = cp.zeros(survcap * 32, dtype=cp.int32)
    ctri = cp.zeros(3, dtype=cp.int32)
    ctr = cp.zeros(3 + NS + 1, dtype=cp.uint64)
    keyg = cp.asarray(np.array(key, dtype=np.uint64))
    tgtg = cp.asarray(np.array(tgt, dtype=np.uint64))
    k_root((1,), (1,), (stack, np.uint64(words[0]), np.uint64(words[1]), np.uint64(words[2]), np.int32(n0)))

    dev = cp.cuda.Device()
    nsm = dev.attributes['MultiProcessorCount']
    bd = 128
    shmem = TW * bd * 4
    grid = nsm * 3
    cap = 20000
    slice_steps = int(os.environ.get('E13_SLICE', '32'))  # steps per round (env, for tuning)
    b = 256
    top, peak, launches, prior = 1, 1, 0, 0.0
    ckpt = os.environ.get('E13_CKPT')
    ckpt_sec = float(os.environ.get('E13_CKPT_SEC', '1800'))
    config = ' '.join(sys.argv[2:]) + f' {os.path.basename(path)} NB={NB}'
    slot = 'a'
    chunk = 1 << 28
    if ckpt and os.path.exists(ckpt):
        # the pointer file names the last complete slot (a or b)
        with open(ckpt) as f:
            old = f.read().strip()
        with np.load(f'{ckpt}.{old}.npz') as npz:  # closed after, so the file can be removed
            z = {name: npz[name] for name in npz.files}
        if str(z['config']) != config:
            sys.exit(f'checkpoint is for another run: {z["config"]}')
        top, peak, launches, prior = int(z['top']), int(z['peak']), int(z['launches']), float(z['wall'])
        with open(f'{ckpt}.{old}.bin', 'rb') as f:
            for o in range(0, top * NB, chunk):
                part = np.frombuffer(f.read(min(chunk, top * NB - o)), dtype=np.uint8)
                stack[o:o + part.size] = cp.asarray(part)
        if os.path.getsize(f'{ckpt}.{old}.bin') != top * NB:
            sys.exit('checkpoint stack file has the wrong size')
        ctr[:] = cp.asarray(z['ctr'])
        ctri[2] = int(z['nsurv'])
        surv[:z['surv'].size] = cp.asarray(z['surv'])
        slot = 'b' if old == 'a' else 'a'
        print(f'resumed from {ckpt}.{old}: stack={top} wall so far {prior:.0f}s', file=sys.stderr, flush=True)
    t_start -= prior
    last = last_ckpt = time.time()
    print(f'X={X} nq={nq} free={len(sys.argv) - 5} sym={sym_lo}..{sym_hi} n0={n0} node={NB}B '
          f'grid={grid}x{bd}', file=sys.stderr, flush=True)
    if ckpt and os.path.exists(ckpt + '.stop'):
        os.remove(ckpt + '.stop')  # left over from the last clean stop
    last_stop_check = 0.0
    while top > 0:
        stop = False
        if ckpt and time.time() - last_stop_check > 1:
            last_stop_check = time.time()
            stop = os.path.exists(ckpt + '.stop')
        if ckpt and (stop or time.time() - last_ckpt > ckpt_sec):
            # write the other slot in full, then point to it: a crash mid-save keeps the old one
            ns = int(ctri[2].get())
            with open(f'{ckpt}.{slot}.bin', 'wb') as f:
                for o in range(0, top * NB, chunk):
                    f.write(stack[o:min(top * NB, o + chunk)].get())
                f.flush()
                os.fsync(f.fileno())
            np.savez(f'{ckpt}.{slot}.npz', config=config, top=top, peak=peak, launches=launches,
                     wall=time.time() - t_start, ctr=ctr.get(), nsurv=ns, surv=surv[:min(ns, survcap) * 32].get())
            with open(ckpt + '.tmp', 'w') as f:
                f.write(slot)
                f.flush()
                os.fsync(f.fileno())
            os.replace(ckpt + '.tmp', ckpt)
            slot = 'b' if slot == 'a' else 'a'
            last_ckpt = time.time()
            if stop:
                os.remove(ckpt + '.stop')
                c = ctr.get()
                print(f'STOPPED cleanly at wall={time.time() - t_start:.0f}s steps={int(c[1]):.3e} '
                      f'survivors={int(c[2])} stack={top}; run the same command to resume',
                      file=sys.stderr, flush=True)
                sys.exit(3)
        nin = min(b, top)
        base_i = top - nin
        ctri[0] = 0
        ctri[1] = 0
        t0 = time.time()
        k_expand((grid,), (bd,), (stack[base_i * NB:], np.int32(nin), out, ctri, ctr, surv, np.int32(survcap),
                                  keyg, tgtg, np.int32(cap), np.int32(slice_steps)), shared_mem=shmem)
        nout = int(ctri[1].get())  # waits for the kernel
        dt = time.time() - t0
        launches += 1
        if base_i + nout > stack_cap:
            sys.exit(f'stack overflow: {base_i + nout} nodes')
        if nout:
            stack[base_i * NB:(base_i + nout) * NB] = out[:nout * NB]
        top = base_i + nout
        peak = max(peak, top)
        # keep each call near 0.1 s (Windows stops GPU calls over ~2 s)
        if dt < 0.05 and b < bmax:
            b *= 2
        elif dt > 0.2 and b > 16:
            b //= 2
        if time.time() - last > 30:
            last = time.time()
            c = ctr.get()
            el = last - t_start
            print(f'  wall={el:.0f}s leaves={int(c[0])} steps={int(c[1]):.3e} ({int(c[1]) / el:.2e}/s) '
                  f'survivors={int(c[2])} stack={top} peak={peak} batch={b}', file=sys.stderr, flush=True)

    cp.cuda.Device().synchronize()
    c = [int(x) for x in ctr.get()]
    ns = int(ctri[2].get())
    sv = surv[:min(ns, survcap) * 32].get().reshape(-1, 32)
    for r in sv:
        cells = []
        for q in range(nq):
            row = []
            for a in (0, 1):
                i = 2 * q + a
                code = int(r[i])
                code = -1 if code >= 62 else code
                row.append(('' if (int(r[26]) >> i) & 1 else '?') + cell_str(code))
            cells.append(' '.join(row))
        print(f'SURVIVOR X={X} st0={r[27]} head0={r[28]} abs={r[29]} rel={r[30]} steps={r[31]} '
              f'table=' + '|'.join(cells))
    sys.stdout.flush()
    wall = time.time() - t_start
    print(f'X={X} budget=0 gpu leaves={c[0]} survivors={c[2]} steps={c[1]} wall={wall:.1f}s '
          f'dead-at-stage(0..{NS}): ' + ' '.join(str(x) for x in c[3:]), file=sys.stderr)
    print(f'  launches={launches} peak stack={peak} nodes ({peak * NB / 1e9:.2f} GB) '
          f'{c[1] / wall:.2e} steps/s', file=sys.stderr)
    if ns > survcap:
        print(f'WARNING: {ns} survivors, only {survcap} printed', file=sys.stderr)
    if ckpt:  # finished: nothing left to resume
        for name in (ckpt, f'{ckpt}.a.bin', f'{ckpt}.a.npz', f'{ckpt}.b.bin', f'{ckpt}.b.npz'):
            if os.path.exists(name):
                os.remove(name)


if __name__ == '__main__':
    main()
