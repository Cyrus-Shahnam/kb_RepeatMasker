#!/usr/bin/env python3
"""
Name-chain check (lessons §2, §17). Run from the module root after every `make`.
Exit 0 = consistent. Also flags init leftovers and unfilled placeholders.
"""
import glob
import json
import os
import re
import sys

import yaml

errors, warnings = [], []
k = yaml.safe_load(open("kbase.yml"))["module-name"].strip()
spec = open(k + ".spec").read()
impl = open("lib/%s/%sImpl.py" % (k, k)).read()
names = {"kbase.yml": k,
         "KIDL": re.search(r"^\s*module\s+(\w+)", spec, re.M).group(1),
         "Impl class": re.search(r"^class\s+(\w+)", impl, re.M).group(1)}
mk = open("Makefile").read()
m = re.search(r"^SERVICE\s*=\s*(\S+)", mk, re.M)
if m:
    names["Makefile SERVICE"] = m.group(1)
if re.search(r"--java", mk):
    errors.append("Makefile still has --java flags (strip them, lessons §3)")
funcdefs = re.findall(r"^\s*funcdef\s+(\w+)", spec, re.M)
mapped = set()
for d in sorted(glob.glob("ui/narrative/methods/*/")):
    sj = json.load(open(os.path.join(d, "spec.json")))
    smap = sj["behavior"]["service-mapping"]
    names["service-mapping " + d] = smap["name"]
    meth = smap["method"]
    mapped.add(meth)
    if meth not in funcdefs:
        errors.append("%s: method %s has no funcdef (init leftover?)" % (d, meth))
    if "#BEGIN " + meth not in impl or "#END " + meth not in impl:
        errors.append("Impl lacks #BEGIN/#END markers for " + meth)
    if not re.search(r"def %s\(" % meth, impl):
        errors.append("Impl lacks def " + meth)
    sjp = {p["id"] for p in sj["parameters"]}
    dy = set(yaml.safe_load(open(os.path.join(d, "display.yaml")))["parameters"])
    im = {x["input_parameter"] for x in smap["input_mapping"] if "input_parameter" in x}
    if sjp != dy:
        errors.append("%s display.yaml mismatch %s" % (d, sjp ^ dy))
    if not sjp <= im:
        errors.append("%s unwired params %s" % (d, sjp - im))
for f in funcdefs:
    if f not in mapped:
        warnings.append("funcdef %s has no UI method directory" % f)
if len(set(names.values())) != 1:
    errors.append("names differ: %s" % names)
for path in ("kbase.yml", "Dockerfile") + tuple(glob.glob("ui/narrative/methods/*/spec.json")):
    if "REPLACE_WITH" in open(path).read():
        warnings.append("placeholder not filled in: %s" % path)
for w in warnings:
    print("WARN  " + w)
print("\n".join("ERROR " + e for e in errors) or "all consistent: %s %s" % (k, funcdefs))
sys.exit(1 if errors else 0)
