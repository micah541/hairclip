"""Shared machinery for the kit v5 builders (kit_step3_v5_build/build_v5.py, kit_step4b_v5_build/build_v5.py).

v5 = the v4 kit TREE (byte copy; the v4 tree is asserted equal to its zip first) + string patches (each `old` must occur exactly once in the
file) + the changed lines of inputs/MANIFEST.sha256 rewritten + CHANGES_v5.md + a deterministic zip.  The v4 tree, zip and build dir are only read.
Runs no certified computation (copy, string replace, sha256, ast.parse, zip)."""
import os, sys, ast, shutil, hashlib, zipfile

S = '/home/claude/s1spec'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''): h.update(c)
    return h.hexdigest()


def sha_b(b): return hashlib.sha256(b).hexdigest()


class Build:
    def __init__(self, v4_tree, v4_zip, v4_zip_sha, zip_prefix_v4, v5_tree, v5_zip, zip_prefix_v5, date_time):
        self.v4, self.v5, self.zip4, self.zip5 = v4_tree, v5_tree, v4_zip, v5_zip
        self.pre4, self.pre5, self.dt = zip_prefix_v4, zip_prefix_v5, date_time
        assert sha(v4_zip) == v4_zip_sha, ('v4 zip changed', v4_zip)
        self.changes = []            # dict(file, before, after, why, patches=[tags])
        self.pins = []               # dict(file, manifest, before, after)
        self.added = []

    def check_v4_tree_equals_zip(self):
        Z = zipfile.ZipFile(self.zip4); n = 0
        tree = set(os.path.relpath(os.path.join(dp, f), self.v4).replace(os.sep, '/') for dp, _, fs in os.walk(self.v4) for f in fs)
        zs = set()
        for nm in Z.namelist():
            r = nm[len(self.pre4):]; zs.add(r); n += 1
            assert sha_b(Z.read(nm)) == sha(os.path.join(self.v4, r)), ('v4 tree differs from its zip', r)
        assert tree == zs, ('v4 tree and zip have different file sets', sorted(tree ^ zs)[:5])
        print(f'v4 tree == v4 zip ({n} files)')

    def copy(self):
        assert not os.path.exists(self.v5), ('v5 tree exists; remove it first (it is a build product)', self.v5)
        shutil.copytree(self.v4, self.v5, symlinks=False)
        print('copied', self.v4, '->', self.v5)

    def patch(self, rel, pairs, why):
        """pairs: [(old, new, tag)]; every old must occur exactly once in the CURRENT file text"""
        p = os.path.join(self.v5, rel); raw = open(p, 'rb').read(); before = sha_b(raw)
        assert b'\r\n' not in raw, ('CRLF in a kit source', rel)
        t = raw.decode('utf-8'); tags = []
        for old, new, tag in pairs:
            c = t.count(old)
            assert c == 1, (rel, tag, 'old string occurs %d times' % c)
            t = t.replace(old, new); tags.append(tag)
        b = t.encode('utf-8'); open(p, 'wb').write(b)
        if rel.endswith('.py'): ast.parse(t, filename=rel)
        prev = [c for c in self.changes if c['file'] == rel]
        if prev: prev[0]['after'] = sha_b(b); prev[0]['patches'] += tags; prev[0]['why'] += ' ' + why
        else: self.changes.append(dict(file=rel, before=before, after=sha_b(b), why=why, patches=tags))

    def add_file(self, rel, text, why):
        p = os.path.join(self.v5, rel); assert not os.path.exists(p), ('exists', rel)
        open(p, 'wb').write(text.encode('utf-8')); self.added.append(dict(file=rel, sha256=sha(p), why=why))

    def regen_manifest_pins(self, man_rel):
        """rewrite the line of every CHANGED file in a 'sha256  relpath' manifest; record each"""
        p = os.path.join(self.v5, man_rel); lines = open(p, 'rb').read().decode('utf-8').split('\n'); out = []
        changed = {c['file']: c for c in self.changes}
        for l in lines:
            if l and not l.startswith('#') and '  ' in l:
                h, r = l.split('  ', 1)
                if r in changed:
                    assert h == changed[r]['before'], ('manifest pin was not the v4 file hash', r)
                    self.pins.append(dict(file=r, manifest=man_rel, before=h, after=changed[r]['after'])); l = changed[r]['after'] + '  ' + r
            out.append(l)
        before = sha(p); open(p, 'wb').write('\n'.join(out).encode('utf-8'))
        self.changes.append(dict(file=man_rel, before=before, after=sha(p), why='pins of the changed kit files regenerated (kit v5)', patches=['pins']))
        # every manifest line verifies against the v5 tree
        bad = []
        for l in out:
            if l and not l.startswith('#') and '  ' in l:
                h, r = l.split('  ', 1); q = os.path.join(self.v5, *r.split('/'))
                if not os.path.exists(q) or sha(q) != h: bad.append(r)
        assert not bad, ('manifest does not verify after regeneration', bad[:5])
        print(f'{man_rel}: {len(self.pins)} pins regenerated; all {sum(1 for l in out if l and not l.startswith("#"))} lines verify')

    def verify_unchanged_rest(self):
        """every file of the v5 tree that was not patched/added is byte-identical to v4"""
        touched = {c['file'] for c in self.changes} | {a['file'] for a in self.added}
        n = 0
        for dp, _, fs in os.walk(self.v5):
            for f in fs:
                r = os.path.relpath(os.path.join(dp, f), self.v5).replace(os.sep, '/')
                if r in touched: continue
                assert sha(os.path.join(dp, f)) == sha(os.path.join(self.v4, r)), ('unexpected change', r); n += 1
        print(f'{n} untouched files byte-identical to v4; {len(touched)} changed/added')

    def all_py_parse(self):
        bad = []
        for dp, _, fs in os.walk(self.v5):
            for f in fs:
                if f.endswith('.py'):
                    try: ast.parse(open(os.path.join(dp, f), encoding='utf-8').read())
                    except SyntaxError as e: bad.append((f, str(e)))
        assert not bad, bad
        print('every .py parses')

    def zip(self):
        names = sorted(os.path.relpath(os.path.join(dp, f), self.v5).replace(os.sep, '/') for dp, _, fs in os.walk(self.v5) for f in fs)
        if os.path.exists(self.zip5): os.remove(self.zip5)
        with zipfile.ZipFile(self.zip5, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            for r in names:
                zi = zipfile.ZipInfo(self.pre5 + r, date_time=self.dt); zi.compress_type = zipfile.ZIP_DEFLATED; zi.external_attr = 0o644 << 16
                z.writestr(zi, open(os.path.join(self.v5, r), 'rb').read(), compresslevel=9)
        self.nfiles = len(names); self.zsha = sha(self.zip5); self.zsize = os.path.getsize(self.zip5)
        # the zip unpacks to exactly the tree
        Z = zipfile.ZipFile(self.zip5)
        for nm in Z.namelist(): assert sha_b(Z.read(nm)) == sha(os.path.join(self.v5, nm[len(self.pre5):]))
        print('zip', self.zip5, self.zsize, self.zsha, self.nfiles, 'files')

    def changes_md(self, title, intro, audit_md):
        L = [f'# {title}', '', intro, '', '## Changed kit files (sha256 before = v4, after = v5)', '',
             '| file | v4 sha256 | v5 sha256 | patches | why |', '|---|---|---|---|---|']
        for c in self.changes:
            L.append(f"| `{c['file']}` | `{c['before']}` | `{c['after']}` | {', '.join(c['patches'])} | {c['why']} |")
        L += ['', '## Pins regenerated in `inputs/MANIFEST.sha256`', '']
        if self.pins:
            L += ['| pinned file | old pin (v4) | new pin (v5) | why |', '|---|---|---|---|']
            for q in self.pins:
                why = next(c['why'] for c in self.changes if c['file'] == q['file'])
                L.append(f"| `{q['file']}` | `{q['before']}` | `{q['after']}` | {why} |")
        else: L.append('none')
        L += ['', '## Added files', ''] + [f"- `{a['file']}`: {a['why']}" for a in self.added] + ['- `CHANGES_v5.md` (this file)', '', audit_md, '']
        return '\n'.join(L)


# ------------------------------------------------------------------------------------------------ patches common to both kits
GUARD_SPLIT_DEF = '''def _win_split(cmd):   # [kit patch v5, winfix 1] tested on Windows 11 / Python 3.13.7 (kit3v4 + kit4bv4 Windows run copies, Nov 24 2026)
    import ctypes
    from ctypes import wintypes
    f = ctypes.windll.shell32.CommandLineToArgvW; f.restype = ctypes.POINTER(wintypes.LPWSTR)
    n = ctypes.c_int(0); r = f(cmd, ctypes.byref(n))
    try: return [r[i] for i in range(n.value)]
    finally: ctypes.windll.kernel32.LocalFree(r)
'''


def patch_guard(B):
    B.patch('kit_guard.py', [
        ('"""kit_guard.py -- [kit v4] run a Python program with every WRITE confined to the kit.\n',
         '"""kit_guard.py -- [kit v4] run a Python program with every WRITE confined to the kit.   [kit patch v5: Windows command-line split, see _win_split]\n', 'docstring'),
        ('def _hook(ev, args):\n', GUARD_SPLIT_DEF + 'def _hook(ev, args):\n', 'winfix1-def'),
        ("        a = list(argv) if isinstance(argv, (list, tuple)) else [argv]\n",
         "        a = list(argv) if isinstance(argv, (list, tuple)) else [argv]\n"
         "        if os.name == 'nt' and len(a) == 1 and isinstance(a[0], str):   # [kit patch v5, winfix 1] Windows audits the joined command line\n"
         "            a = _win_split(a[0])\n", 'winfix1-split'),
    ], 'winfix 1: on Windows the subprocess.Popen audit event carries ONE joined command line, so the guard denied the kit\'s own guarded Python '
       'children; split it with CommandLineToArgvW (code of the tested Windows run copies). No effect on POSIX.')


CHART4 = 'writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py'
RSS_FALLBACK = '''try:
    import resource   # Unix only; used for the RSS report (timing sidecar and log lines, not compared, not part of the proof)
except ImportError:   # [kit patch v5, winfix 2] Windows: peak working set via psapi (as the tested kit4bv4 Windows env shim), else -1
    resource = None
'''
RSS_DEF_OLD = "def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024\n"
RSS_DEF_NEW = '''def rss():
    if resource is not None: return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024
    try:   # [kit patch v5, winfix 2] Windows: PeakWorkingSetSize in MB; -1 if unavailable
        import ctypes, ctypes.wintypes as _w
        class _PMC(ctypes.Structure):
            _fields_ = [("cb", _w.DWORD), ("PageFaultCount", _w.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        c = _PMC(); c.cb = ctypes.sizeof(c); k = ctypes.windll.kernel32; k.GetCurrentProcess.restype = ctypes.c_void_p
        f = ctypes.windll.psapi.GetProcessMemoryInfo; f.argtypes = [ctypes.c_void_p, ctypes.c_void_p, _w.DWORD]
        return c.PeakWorkingSetSize // (1 << 20) if f(k.GetCurrentProcess(), ctypes.byref(c), c.cb) else -1
    except Exception:
        return -1
'''


def patch_chart4(B):
    B.patch(CHART4, [('import resource\n', RSS_FALLBACK, 'winfix2-import'), (RSS_DEF_OLD, RSS_DEF_NEW, 'winfix2-rss')],
            'winfix 2: `resource` is Unix-only (RSS report only, written to the timing sidecar / log, never compared): import guarded, Windows peak working '
            'set via psapi, else -1. Replaces the env shim of the kit4bv4 Windows run. Certified computation unchanged.')
