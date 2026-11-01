"""Classify deep.c output lines. A real Erdos machine must keep the counter of c(n) at every
stage start: c = 0 for n = 1, 2; 1 for n = 3..8; 2 for n = 9..NMAX (2^9 .. 2^NMAX-1 all have a
ternary digit 2). Usage: python classify_deep.py SURVIVORS DEEPOUT [DEEPOUT ...]
(deep outputs in the same order as the survivor lines)."""
import re
import sys
from collections import Counter


def c_of(n):
    return 0 if n <= 2 else (1 if n <= 8 else 2)


surv = [l.rstrip('\n') for l in open(sys.argv[1]) if 'table=' in l]
deep = [l.rstrip('\n') for f in sys.argv[2:] for l in open(f) if l.startswith('survivor')]
assert len(surv) == len(deep), (len(surv), len(deep))
kinds = Counter()
examples = {}
for s, d in zip(surv, deep):
    body = d.split('):', 1)[1]
    stages, _, tail = body.partition('|')
    ctr = {int(m.group(1)): m.group(2) for m in re.finditer(r'n(\d+):(\S+)', stages)}
    first = {}
    kind = None
    for n in sorted(ctr):
        c = c_of(n)
        if c in first and ctr[first[c]] != ctr[n]:
            kind = f'counter changes at n={n} ({ctr[first[c]]} -> {ctr[n]})'
            break
        first.setdefault(c, n)
    if kind is None and len(set(ctr[first[c]] for c in first)) != len(first):
        kind = 'counter equal for different c'
    if kind is None and tail.strip():
        kind = re.sub(r'\(t=\d+\)|\(limit, t=\d+\)|at -?\d+', '', tail.strip())
    if kind is None:
        kind = f'CANDIDATE (consistent to n={max(ctr)})'
    rel = 'rel' if ' rel)' in d else 'abs'
    kinds[(kind, rel)] += 1
    examples.setdefault((kind, rel), s)
for (kind, rel), k in sorted(kinds.items(), key=lambda x: -x[1]):
    print(f'{k:6d}  {rel}  {kind}')
for key, s in examples.items():
    if key[0].startswith('CANDIDATE') or 'FREE' in key[0]:
        print('example', key, s)
