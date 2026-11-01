# Searching for a 13-state Erdős powers-of-two machine: a negative result

**Goal.** Find a 2-symbol Turing machine with 13 states that halts if and only if some 2^n with n > 8
has no digit 2 in base 3. Erdős conjectured that no such n exists. The smallest known 2-symbol machine
is Tristan Stérin's 14-state `erdos_14` (October 2026), down from the 15-state machine of Stérin and
Woods (arXiv:2107.12475). Our Lean proof that `erdos_14` is correct is in
[google-deepmind/cryptids#19](https://github.com/google-deepmind/cryptids/pull/19), which is open and
awaiting review.

**Result.** No 13-state machine was found. Within the search rules below, none is "one module away"
from `erdos_14`:

- with the other two modules kept as in `erdos_14`, the counter and end of stage cannot be done in
  3 states with a cell left for the halt, the doubling in 2, or the division in 6 (for any value of
  the dispatch cell `C0`);
- deleting any one state and changing up to 4 cells (the rewritten rules into it included) gives
  nothing;
- changing up to 3 cells around any of the 61 three-state counters gives nothing, and neither does
  re-searching the counter together with one of the cells `C1`, `D1`, `C0` or `G1`.

This is **not** a lower bound. A 13-state machine with a different tape layout, a different stage
structure, or a much slower division would not be seen by these searches. The total work was about
1.4e16 simulated steps, mostly on GPUs.

## How `erdos_14` works

```
1LB1RA_1LC1LB_0RD0RA_1LM0RE_1LF1RE_0LJ1LG_1RK0LH_1RI1LH_0LJ0RE_0RD1LJ_0LL1RK_1RN0LC_0LN---_0RK0LM
```

The tape holds a separator 0 at cell −2, the number 2^n in unary as a block of ones at cells
−1 … 2^n − 2, and a small counter left of the separator. **Stage n** starts in state `C` at cell −2
with the block for 2^n and checks that power of two:

| Module | States | What it does |
|---|---|---|
| doubling + boot | `A B C` | `A`, `B` write the first block; later stages double the block in place. `C0` (at the separator) dispatches to the division. |
| division | `D`–`J` | Reads the base-3 digits of the block length by repeated division by 3, using two holes in the block. A digit 2 ends the stage (`G0`), quotient 0 means "no digit 2" (`D0`). |
| counter + end of stage | `K L M N` | On an exception, advances the counter (cells −3, −4, −5 = `000`, `100`, `111` after 0, 1, 2 exceptions); the third exception halts (`M1`, the only halt). Then removes the extra cell left by the division and hands over to `C`. |

An **exception** is a power of two with no digit 2 in base 3. The known ones (for n ≥ 1) are
2^2 = 11 and 2^8 = 100111, so the counter is at 0 for stages 1 and 2, at 1 for stages 3 … 8 and at
2 from stage 9 on. `machines/erdos_14.txt` has the table with states `A`–`N`
numbered 0–13.

## The search

**Candidates.** Delete one state of `erdos_14` (rules into it become free), then either let a chosen
set of cells take any move, or allow up to FMAX changed cells anywhere (the rewritten rules into the
deleted state count towards FMAX); everything else stays as in `erdos_14`. A cell is one state
reading one symbol, e.g. `D0`. A free cell gets its value when the
machine first reads it, so machines that differ only in cells never read count as one **leaf** (the
usual tree-normal-form enumeration). A halt is never chosen at a read, since that would halt too
early; it goes into a cell that is still unread at the end.

**Filter.** Run from the blank tape. For n = 1 … 9 in order, the machine must reach a **stage start**:
a step at which the ones at cells ≥ −2 are exactly the block −1 … 2^n − 2 (the first such step with
the right state and head counts). Then:

- the state must be the same at every stage start, and the head either at the same cell every time
  or at the same distance from the block end;
- stage n must be reached within 4 × the step at which `erdos_14` reaches it, plus 64 (counted from
  the blank tape; 208,864 steps for stage 9 in `erdos_14`);
- at each stage start the cells left of −2 (the counter) must be the same for all stages with the
  same number of exceptions so far, and different for different numbers;
- the machine must stay within 64 cells left and 1300 right of the start, and must not halt.

A **survivor** reaches the start of stage 9. The block and counter regions are compared by 64-bit
hashes. A hash collision could drop a machine or keep a wrong one; over about 1.4e16 steps in all,
the expected number of collisions is below 0.001. Every survivor that could still hold a halt was
re-checked by plain simulation, without hashing.

**Symmetry.** When a set of states is fully free (for the division: `D`–`I`, with `J` deleted), the
states are interchangeable, and a move may only name a state already in use or the next unused one.
Each machine is then met once up to relabelling of those states.

**Survivor checks.** Survivors go through up to three independent checkers:

- `check.py` (no hashing) re-runs stages 1 … 11, then starts the machine from crafted stage starts
  (counter at 0, 1, 2; block length m = 2 … 40). It must halt exactly when the counter is at 2 and m
  has no digit 2, and otherwise reach the next stage.
- `deep.c` runs stages 1 … 15 or 16 and prints the number and leftmost position of the counter ones
  at each stage start.
- `recheck.c` (written separately, also without hashing) does the same and compares the exact
  counter contents.

## Results

| Search | Free part | Leaves | Time | Survivors | Usable |
|---|---|---|---|---|---|
| **1. Local edits:** delete any one state, change ≤ 2 / ≤ 3 / ≤ 4 cells | any cell | 1.5e5 / 3.6e7 / 6.0e9 | 0.1 s / 18 s / 41 min (CPU) | 1 / 11 / 57 | 0 |
| **2. Counter + end in 3 states** (delete `N`) | `D0 G0 K0 K1 L1 M1` and the rules into `N` | 1.85e10 (3.08e9 with symmetry) | 52 min (CPU) | 366 (61 up to relabelling) | 0 |
| **3. Doubling + boot in 2 states** (delete `A`, `B` or `C`) | the other two rows and the rules into the deleted state (e.g. `L1` for `C`) | 1.6e5 / 1.6e5 / 1.0e6 | 0.1 s (CPU) | 0 | 0 |
| **4. Division in 6 states** (delete `J`, `D`–`I` interchangeable) | | | | | |
| V4: `C0` and the `D`, `E` rows kept | 8 cells | 8.3e8 | 48 s (CPU), 4.7 s (GPU) | 0 | 0 |
| V3: `C0` and the `D` row kept | 10 cells | 5.3e10 | 79 min (CPU), 139 s (GPU) | 0 | 0 |
| V2: `C0` kept | 12 cells | 3.13e13 | 20.1 h (RTX 3090) | 12 | 0 |
| V1: `C0` free too (branches `C0` = `1LD`, `0LD`, `1RD`) | 13 cells | 3.71e13 + 4.14e13 + 4.23e13 | 9.6 + 11.4 + 10.5 h (RTX 5090) | 0 + 0 + 2868 | 0 |
| **5a. Spare division cell:** ≤ 2 / ≤ 3 changed cells around each of the 61 three-state counters | any cell | 1.85e7 / 2.85e9 | 34 s / 33 min (CPU) | 1074 / 2115 | 0 |
| **5b. Counter + one more cell** (`C1`, `D1`, `C0` or `G1` free as well) | 9 cells | 9.1e10 – 2.7e11 each | 23 min in all (RTX 3090) | 61 / 70 / 185 / 61 | 0 |

The steps for the division searches are 9.7e10 (V4), 5.8e12 (V3), 3.30e15 (V2) and 1.10e16 (V1).

**Why the survivors are not usable.**

- *Searches 1, 2 and 5:* almost all survivors read every cell before stage 9, so no cell is left for
  the halt. The rest leave a cell unread: 5 in search 1, 30 and 115 in 5a (≤ 2 and ≤ 3 changes), 1
  in 5b. With the halt put in that cell (in each unread cell in turn, where there are two), every one
  of them counts correctly to stage 11, but none halts from `check.py`'s crafted starts with the
  counter at 2 and a block with no digit 2: it runs on (no halt within 40m² + 5000 steps) or starts
  the next stage. A usable machine must halt at the stage of any 2^n without a digit 2. No such 2^n
  with n > 8 is known, so crafted starts are the only test of that branch, and we count a failure
  there as unusable. So 0 pass.
- *Search 2:* all 61 counters use all six cells of their three states. This is why search 5 looks
  for a division that leaves one of its own cells free for the halt.
- *Division (V2 and V1):* all 2880 survivors (12 from V2, 2868 from V1) fail within a few stages:
  2877 break the counter rule and 3 run away. All 12 of V2 and 2856 of V1 change their counter at
  the start of stage 15, as if 2^14 = 211110211 (base 3) had no digit 2; 9 of V1 do it at stage 13
  (2^12 = 12121201). The other 3 move off to the left and do not reach stage 12 within 1.08e9 steps
  (80 × `erdos_14`, checked with `deep_big.c`). Their rows `A`–`I` never name `K`–`N`, so they can
  never reach the halt. No survivor reads a never-read free cell before it fails, so the verdicts
  hold for every way of filling those cells. The 12 V2 survivors also fail 107 of `check.py`'s 117 crafted stage tests.
- *The stage-15 pattern:* the 2868 machines that change at stage 15 (12 from V2, 2856 from V1) count
  exceptions at 2^2, 2^8 and 2^14. That is the pattern
  of a test for 2^k ≡ 4 (mod 9), i.e. last two base-3 digits `11`, which holds exactly for
  k ≡ 2 (mod 6). This is our reading of them; we did not prove it. A 9-stage filter cannot tell such
  a machine from a real one, which is why the deep checks exist.

**`C0` is covered completely.** `C0` is the first free cell V1 reads (at step 2, so it cannot be the
halt). Its 32 values are the 4 moves to each of 8 targets: `A`, `B`, `C`, `K`–`N`, and `D`, which
stands for any of the interchangeable `D`–`I`. The 28 that do not name `D` never reach another free
cell, since the fixed rows never name `D`–`I`; each is a single leaf that dies within 153 steps
(`v1c0_est.txt`). The four that name `D` are `0RD` (= V2) and the three V1 branches.

## What this does and does not show

Read together, the one-module searches suggest that in this layout a 13-state machine would need 3 states
for the doubling, 7 for the division and 3 for the counter: 26 cells, one of which must be the halt.
The three-state counters use all their cells, and in `erdos_14` all six doubling cells are in use, so
the halt would most likely have to sit in a division cell that the division does not need. Search 5a
looked for that within 3 changed cells of the 61 counters and found nothing. This is not a proof, since each
module was shrunk with the other two fixed.

Not covered:

- **Two modules changed at once.** The doubling and counter together are estimated at 5e17 steps
  in the cheapest case (100+ days on an RTX 3090) and up to about 1e19; freeing only the `K` and `L` rows on top
  of V2 is estimated at 2.7e22 steps. These are out of reach for this search.
- **Other layouts.** The filter assumes `erdos_14`'s tape layout: block at −1 … 2^n − 2, separator at
  −2, counter on the left and a function of the exception count.
- **Slower or wider machines.** Anything that needs more than 4 × `erdos_14`'s steps (+64), counted
  from the blank tape, to reach some stage n ≤ 9, or that leaves the window of 64 cells left and 1300
  right before stage 9.
- **Early stage-start configurations.** The first step with the right state, head and exact block
  is taken as the stage start, and a wrong counter there rejects the machine. A machine that passes
  through that configuration before its counter is final is not covered.
- **Free boot.** When cells of `A`, `B` or `C` can change (searches 1, 3 and 5a), the filter fixes the
  stage-start state and head at the first moment the 2-block appears. A boot that passes the 2-block
  in another state and comes back later is not covered. V1 and 5b are not affected: their free cells
  `C0` and `C1` are first read after that moment.
- **Larger counter extensions.** The counter re-searched together with a whole row (`A`, `B`, `E`,
  `F`, `H`, `I` or `J`), or with `C0 C1` or `F0 I0`, was not run (about 2 GPU days in all).

Things tried that did not help: SAT synthesis of the free cells (CaDiCaL on the unrolled run) was far
slower than the search here, giving no answer in 10 minutes on cases the search finishes in 48 s.
Several hand-designed layouts (quotient on the far side of the gap, a unary counter behind a fixed 0,
a sparse block, a mirrored division, shared scan states) each fail on paper or in simulation.

## Validation

- `erdos_14` passes the filter (stage 9 at step 208,864) and `check.py`
  (`results/validation/base_erdos14.txt`, `check_erdos14.txt`).
- All 1568 single-cell changes of `erdos_14` were run (`results/validation/mutants.txt`). 61 pass:
  56 only change the halt cell `M1`, which is never read before stage 9; the other 5 are `C0` = `0RI`,
  `0RJ` and `N1` = `0LK`, `0LL`, `0LN`.
- Split tests: a 15-state copy of `erdos_14` with a duplicate of `D` (`machines/split_15.txt`).
  Deleting the duplicate and allowing 1 change gives back `erdos_14` plus the two `C0` variants
  (`split_test.txt`). So does deleting it and freeing `K0 K1 L1` in free mode, the code path of V1
  and V2 (`split_test_free.txt`).
- The GPU search gives exactly the CPU counts (leaves, steps, deaths per stage) on V4, V3 and the
  symmetric counter search (`log_v4.txt` / `g_v4.log`, `log_v3.txt` / `g_v3.log`, `log_modsym.txt`
  / `g_modsym.log`). A stop and resume from a checkpoint also gave exact counts in development tests
  (those logs were not kept); the V2 run itself was stopped and resumed once.
- CPU logs without a `split=` field come from earlier builds of the same program (the output format
  changed over time, and the X = −1 option was added for search 5a). The shipped `src/e13search.c`
  reproduces V4, the base check and both split tests exactly (re-run on 1 November 2026,
  `results/validation/v4_rerun.txt` and the files above).
- With symmetry the counter search gives 61 = 366 / 3! survivors; the 366 reduce to exactly these 61
  under relabelling of `K`, `L`, `M`.
- `check.py` caught a bug in an early version of the search (a wrong counter at the first stage
  match was not rejected). `deep.c` had a 32-bit overflow under MSVC and first ignored stage starts
  relative to the block end; both were fixed and re-run before the results here.
- `recheck.c` reproduces all 2868 V1 verdicts, including the split into the classes above.

## Reproducing

The CPU code was built with MSVC (`cl /O2 /openmp src\e13search.c`, `cl /O2 src\deep.c`,
`cl /O2 src\recheck.c`) and run on 4 to 20 threads of one desktop PC. `deep.c` was also built with `gcc -O2` on Linux;
`e13search.c` was not tried with gcc. The GPU search needs Python 3, NumPy and CuPy; it ran on an
RTX 3090 (Windows) and a rented RTX 5090 (Linux, CUDA 13.0, CuPy 14.2.0). States are numbered from 0
(`A` = 0 … `N` = 13).

```
# 1. local edits (FMAX = 2, 3, 4)
e13search delete machines/erdos_14.txt FMAX 0 1 2 3 4 5 6 7 8 9 10 11 12 13
# 2. counter in 3 states, without and with symmetry
e13search free    machines/erdos_14.txt 13 D0 G0 K0 K1 L1 M1
e13search freesym machines/erdos_14.txt 13 10 13 D0 G0 K0 K1 L1 M1
# 3. doubling in 2 states
e13search free machines/erdos_14.txt 2 A0 A1 B0 B1
e13search free machines/erdos_14.txt 1 A0 A1 C0 C1
e13search free machines/erdos_14.txt 0 B0 B1 C0 C1
# 4. division in 6 states: V4, V3 (CPU), V2 and V1 (GPU)
e13search freesym machines/erdos_14.txt 9 3 9 F0 F1 G0 G1 H0 H1 I0 I1
e13search freesym machines/erdos_14.txt 9 3 9 E0 E1 F0 F1 G0 G1 H0 H1 I0 I1
python src/e13gpu.py machines/erdos_14.txt 9 3 9 D0 D1 E0 E1 F0 F1 G0 G1 H0 H1 I0 I1
python src/e13gpu.py machines/c0_1RD.txt   9 3 9 D0 D1 E0 E1 F0 F1 G0 G1 H0 H1 I0 I1   # also c0_1LD, c0_0LD
# 5a. around each three-state counter (X = -1: no state deleted; FMAX = 2, 3)
e13search delete machines/counters_3state/base_01.txt FMAX -1
# 5b. counter + one cell (also D1, C0, G1)
python src/e13gpu.py machines/erdos_14.txt 13 10 13 K0 K1 L0 L1 M0 M1 D0 G0 C1
# survivor checks
python src/check.py survivors.txt
deep 15 < survivors.txt > deep.txt && python src/classify_deep.py survivors.txt deep.txt
recheck 15 < survivors.txt
```

Survivors go to stdout and the summary to stderr. With `E13_CKPT=file` the GPU search saves a
checkpoint every `E13_CKPT_SEC` seconds (default 1800; 600 in the V1 run) and resumes from it.
X = −1 (no state deleted) works only in the CPU `delete` mode. `src/run_v1.sh` is the script that
ran the three V1 branches on the rented machine (its paths are that machine's).

## Files

- `src/e13search.c`: CPU search (modes `base`, `mutant`, `delete`, `free`, `freesym`, `estimate`).
- `src/e13gpu.py`: GPU version of `free` / `freesym`.
- `src/check.py`, `src/deep.c`, `src/deep_big.c`, `src/recheck.c`, `src/classify_deep.py`: survivor
  checks (`deep_big.c` is `deep.c` with a 4M-cell tape).
- `machines/`: `erdos_14`, the V1 branch tables, the split-test machine and the 61 three-state
  counters (as 13-state machines with `N` deleted).
- `results/1_local_edits` … `results/5_counter_plus`: logs and survivor lines. The totals and deaths
  per stage are on the `leaves=` line (the last line of CPU logs, the one before last of GPU logs;
  multi-state logs have one such line per deleted state, which the table adds up). Older CPU logs
  use an older output format (see Validation).
- **Survivor lines:** `st0` and `head0` are the stage-start state and head, `abs` / `rel` say which
  head rule held (same cell / same distance from the block end), and `table=` lists the rows,
  separated by `|`. States after the deleted one are renumbered (in V1, `K`–`N` are 9–12), and
  interchangeable states are numbered in the order the search first names them. A `?` marks a cell
  the search never read. `?---` is either a free cell (it may hold any move or the halt) or the
  fixed halt `M1`. `?` followed by a move is, in `free` / `freesym` output, a fixed cell that keeps
  that move; in `delete` output it may still be changed within the remaining budget.
- `results/4_division_6/v1`: V1 logs, the 2868 survivors, the deep check (`deep15_1RD.txt`), its
  classification (`classify_1RD.txt`), the 3 runaways (`offtape3.txt`, `runaways_deep_big.txt`) and
  the per-`C0` size estimates (`v1c0_est.txt`).
- `results/validation`: the mutant run (`cell value stage passed`, cell = 2 · state + symbol, value =
  4 · next + 2 · right + write, −1 = halt), the deep check of `erdos_14`, the `recheck.c` output and
  the re-runs of 1 November 2026.

## Credits

`erdos_14` is due to Tristan Stérin; it is listed on the
[Cryptids](https://wiki.bbchallenge.org/wiki/Cryptids) page of the bbchallenge wiki. The 15-state
machine is from Tristan Stérin and Damien Woods ([arXiv:2107.12475](https://arxiv.org/abs/2107.12475)).

Claude (Anthropic) was used as assistant.
