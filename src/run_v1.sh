#!/bin/bash
# V1 on the rented RTX 5090: the three fixed-C0 branches (C0 = 1LD, 0LD, 1RD; 0RD = V2 is done),
# one after another, each with its own checkpoint. Rerunning the script resumes: finished branches
# are skipped (done_<b>), an interrupted one continues from ck_<b> (its old log is kept as a part).
# After each branch: deep stage check to n = 16 (kills mod-9 fakes) and check.py on the survivors.
cd /root/e13
PY=/venv/main/bin/python
export PYTHONIOENCODING=utf-8 E13_CKPT_SEC=600
CELLS="D0 D1 E0 E1 F0 F1 G0 G1 H0 H1 I0 I1"
[ -x deep ] || gcc -O2 -o deep deep.c || exit 1
for b in 1LD 0LD 1RD; do
  [ -e done_$b ] && continue
  if [ -e ck_$b ] && [ -e log_$b.txt ]; then mv log_$b.txt log_${b}_part$(date +%s).txt; fi
  echo "$(date '+%F %T') branch $b start" >> run_v1.log
  E13_CKPT=/root/e13/ck_$b $PY e13gpu.py v1/c0_$b.txt 9 3 9 $CELLS > surv_$b.txt 2> log_$b.txt
  code=$?
  echo "$(date '+%F %T') branch $b exit $code (0 = finished, 3 = clean stop)" >> run_v1.log
  [ $code -eq 0 ] || exit $code
  ./deep 16 < surv_$b.txt > deep_$b.txt
  $PY check.py surv_$b.txt > check_$b.txt 2>&1
  touch done_$b
done
echo "$(date '+%F %T') ALL DONE" >> run_v1.log
