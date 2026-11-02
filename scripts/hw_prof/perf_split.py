#!/usr/bin/env python3
"""Split a multi-event `perf report --stdio` into one symbol table per event.

Usage: perf_split.py FULL_REPORT OUT_DIR
Writes OUT_DIR/by_symbol_<event>.txt (the report's percentage lines for that event) and
OUT_DIR/event_totals.tsv (event, approximate event count). (`perf report --event X` alone gives
empty tables, hence this.)
"""
import os
import re
import sys

full, out = sys.argv[1], sys.argv[2]
tables, totals, ev = {}, {}, None
for line in open(full):
    m = re.match(r"^# Samples: .* of event '([^']+)'", line)
    if m:
        ev = m.group(1).split(":")[0]
        tables[ev] = []
        continue
    m = re.match(r"^# Event count \(approx\.\): (\d+)", line)
    if m and ev:
        totals[ev] = int(m.group(1))
        continue
    if ev and re.match(r"^ +[0-9.]+%\s+\[\.\]", line):   # user-space symbols
        tables[ev].append(line)
os.makedirs(out, exist_ok=True)
for ev, rows in tables.items():
    with open(os.path.join(out, f"by_symbol_{ev}.txt"), "w") as f:
        f.writelines(rows)
with open(os.path.join(out, "event_totals.tsv"), "w") as f:
    for ev, n in totals.items():
        f.write(f"{ev}\t{n}\n")
print("events:", ", ".join(f"{e} ({len(tables[e])} symbols)" for e in tables))
