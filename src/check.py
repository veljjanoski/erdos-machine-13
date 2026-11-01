"""Independent check of e13search survivors (no hashing, no shared code).

For a machine table (cells may carry a '?' = never read by the search run):
  1. Run from the blank tape. Stage n starts at the first step at which the
     ones at cells >= -2 are exactly the block -1 .. 2^n - 2 and the state and
     head match stage 1 (head at the same cell, or at the same offset from the
     block end). Check n = 1 .. NMAX in order, no halt, and record the ones at
     cells < -2 (the counter) at each start. The counter must depend only on
     c(n) = #{1 <= k < n : 2^k has no digit 2 in base 3} and differ between
     values of c.
  2. From crafted stage starts (counter of c = 0, 1, 2 as recorded, block of
     m ones), the machine must halt iff c = 2 and m has no digit 2, and
     otherwise reach the next stage start with block 2m and counter of c or
     c + 1. Run for m = 2 .. MMAX.
The halt goes into a cell marked '?---', or into each '?' cell in turn.
"""
import sys

NMAX = 11
MMAX = 40


def parse(table):
    rows = []
    unread = []
    for q, row in enumerate(table.split('|')):
        cells = []
        for a, cell in enumerate(row.split()):
            if cell.startswith('?'):
                unread.append(2 * q + a)
                cell = cell[1:]
            cells.append(None if cell == '---' else (int(cell[0]), 1 if cell[1] == 'R' else -1, int(cell[2:])))
        rows.append(cells)
    return rows, unread


def no_two(m):
    while m:
        if m % 3 == 2:
            return False
        m //= 3
    return True


def is_block(tape, size, ones_ge):
    if ones_ge != size:
        return False
    return all(k in tape for k in range(-1, size - 1))


def counter(tape):
    return tuple(sorted(k for k in tape if k < -2))


def run_tracked(rows, tape, head, st, stop, limit):
    """Simulate; tape is a set of ones. stop(tape, head, st, ones_ge) is
    tested before each step. Returns (halt|stop|limit, steps, head, state)."""
    ones_ge = sum(1 for k in tape if k >= -2)
    for t in range(limit):
        if stop(tape, head, st, ones_ge):
            return 'stop', t, head, st
        sym = 1 if head in tape else 0
        r = rows[st][sym]
        if r is None:
            return 'halt', t, head, st
        w, d, nxt = r
        if w and not sym:
            tape.add(head)
            if head >= -2:
                ones_ge += 1
        elif sym and not w:
            tape.discard(head)
            if head >= -2:
                ones_ge -= 1
        head += d
        st = nxt
    return 'limit', limit, head, st


def stages(rows):
    tape, head, st = set(), 0, 0
    starts, total = [], 0
    for n in range(1, NMAX + 1):
        size = 1 << n

        def stop(tp, h, s, og, size=size):
            if not is_block(tp, size, og):
                return False
            if not starts:
                return True
            s0, h0 = starts[0][1], starts[0][2]
            return s == s0 and (h == h0 or h - size == h0 - 2)

        res, t, head, st = run_tracked(rows, tape, head, st, stop, 40 * (4 ** n) + 1000)
        total += t
        if res != 'stop':
            return starts, f'stage {n}: {res} at step {total}'
        starts.append((n, st, head, counter(tape), total))
    return starts, 'ok'


def check(table, verbose=True):
    rows0, unread = parse(table)
    halts = [i for i in range(2 * len(rows0)) if rows0[i // 2][i % 2] is None]
    placements = [None] if halts else [i for i in unread]
    out = []
    for place in placements:
        rows = [list(r) for r in rows0]
        if place is not None:
            rows[place // 2][place % 2] = None
        starts, msg = stages(rows)
        if msg != 'ok':
            out.append((place, 'FAIL', msg))
            continue
        # counter consistency
        cnt = {}
        bad = None
        for (n, s, h, ctr, t) in starts:
            c = sum(1 for k in range(1, n) if no_two(2 ** k))
            if c in cnt and cnt[c] != ctr:
                bad = f'counter differs within c={c} at n={n}'
            cnt.setdefault(c, ctr)
        if len(set(cnt.values())) != len(cnt):
            bad = 'counter equal for different c'
        if bad:
            out.append((place, 'FAIL', bad))
            continue
        s0, h0 = starts[0][1], starts[0][2]
        absolute = all(h == h0 for (_, _, h, _, _) in starts)
        # crafted stage starts
        fails = []
        for c in (0, 1, 2):
            if c not in cnt:
                continue
            for m in range(2, MMAX + 1):
                tape = set(cnt[c]) | set(range(-1, m - 1))
                head = h0 if absolute else h0 - 2 + m
                nxt = 2 * m
                want_c = c + 1 if no_two(m) else c
                want_halt = no_two(m) and c == 2
                target = cnt.get(want_c)

                def stop(tp, h, s, og, nxt=nxt, target=target):
                    if not is_block(tp, nxt, og):
                        return False
                    hh = h0 if absolute else h0 - 2 + nxt
                    return s == s0 and h == hh

                res, t, hd, st = run_tracked(rows, tape, head, s0, stop, 40 * m * m + 5000)
                if want_halt:
                    if res != 'halt':
                        fails.append((c, m, f'expected halt, got {res}'))
                elif res != 'stop':
                    fails.append((c, m, f'expected stage 2m, got {res} at {t}'))
                elif want_c in cnt and counter(tape) != target:
                    fails.append((c, m, f'counter {counter(tape)} != c={want_c} {target}'))
                elif want_c not in cnt and counter(tape) in cnt.values():
                    fails.append((c, m, 'counter for c=3 equals an earlier one'))
        status = 'PASS' if not fails else 'FAIL'
        out.append((place, status, f'stages ok to n={NMAX} (step {starts[-1][4]}); counters {cnt}; '
                    f'crafted fails {len(fails)}: {fails[:4]}'))
    return out


if __name__ == '__main__':
    for line in open(sys.argv[1]):
        if 'table=' not in line:
            continue
        table = line.split('table=')[1].strip()
        print(line.strip())
        for place, status, msg in check(table):
            print(f'   halt@{place}: {status} {msg}')
        sys.stdout.flush()
