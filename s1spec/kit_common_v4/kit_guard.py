"""kit_guard.py -- [kit v4] run a Python program with every WRITE confined to the kit.

usage:  python [interpreter flags] kit_guard.py <script.py> [args...]      (or  kit_guard.py -c <code> / -m <module> [args...])

Installs a sys.addaudithook() hook (PEP 578) before the program starts and then runs it as __main__ (runpy, with sys.path[0] = the script's
directory unless the interpreter runs isolated / safe-path, exactly as `python script.py` would).  The hook RAISES PermissionError when the
program (or a library it calls)
  * opens a file for writing / appending / creating, or renames, replaces, removes, truncates, links, chmods, creates or deletes a file or
    directory, at a path outside the kit root (the directory holding KIT_ROOT) and outside the system temporary directory;
  * changes its working directory to a place outside the kit;
  * starts a child process whose working directory lies outside the kit, or a Python child that does not go through this guard.
Python children started with subprocess are routed through this guard automatically (the interpreter command gets kit_guard.py inserted
after the interpreter flags), so the rule holds for grandchildren too.  Reads are not blocked; reads of files outside the kit and outside the
Python installation / OS directories are LOGGED to kit_out/guard/reads_<pid>.txt (check.py lists them as WARN).
Background: on Nov 15 2026 a kit stage resolved its working directory outside the kit and overwrote a file of the proof tree; this guard
makes that impossible for every process the kit starts."""
import os, sys, runpy, tempfile, subprocess

def _root(d=os.path.dirname(os.path.abspath(__file__))):
    while not os.path.exists(os.path.join(d, 'KIT_ROOT')) and d != os.path.dirname(d): d = os.path.dirname(d)
    return d
KIT = os.path.realpath(_root())
GUARD = os.path.join(KIT, 'kit_guard.py')
assert os.path.exists(os.path.join(KIT, 'KIT_ROOT')), 'kit_guard: KIT_ROOT not found'
_TMP = os.path.realpath(tempfile.gettempdir())
_WRITE_OK = [KIT, _TMP]
_READ_OK = sorted({os.path.realpath(p) for p in [KIT, _TMP, sys.prefix, sys.base_prefix, sys.exec_prefix, os.path.dirname(os.__file__)] + list(sys.path) if p} |
                  ({'/usr', '/lib', '/lib64', '/etc', '/proc', '/sys', '/dev', '/bin', '/opt'} if os.name != 'nt' else {os.path.realpath(os.environ.get('SystemRoot', 'C:\\Windows'))}))
_READ_OK = [p for p in _READ_OK if p and p != os.path.sep and os.path.isabs(p) and not (KIT.startswith(p + os.sep) and p not in (sys.prefix, sys.base_prefix, sys.exec_prefix))]
def _abs(p):
    if isinstance(p, bytes): p = os.fsdecode(p)
    if isinstance(p, os.PathLike): p = os.fspath(p)
    if not isinstance(p, str): return None
    return os.path.realpath(p if os.path.isabs(p) else os.path.join(os.getcwd(), p))
def _inside(p, roots):
    q = _abs(p)
    if q is None: return True
    if os.name == 'nt': q = q.lower(); roots = [r.lower() for r in roots]
    return any(q == r or q.startswith(r.rstrip(os.sep) + os.sep) for r in roots) or q in (os.path.realpath(os.devnull),)
def _deny(what, p):
    raise PermissionError(f'kit guard: {what} outside the kit refused: {p!r} (kit root {KIT})')
_READS = set(); _LOGGING = [False]
def _log_read(p):
    q = _abs(p)
    if q is None or q in _READS or _LOGGING[0]: return
    _READS.add(q); _LOGGING[0] = True
    try:
        d = os.path.join(KIT, 'kit_out', 'guard'); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f'reads_{os.getpid()}.txt'), 'a', encoding='utf-8') as f: f.write(q + '\t' + ' '.join(sys.argv[:2]) + '\n')
    finally: _LOGGING[0] = False
_WFLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC
def _is_py(exe):
    if not isinstance(exe, (str, bytes, os.PathLike)): return False
    e = _abs(exe) or ''
    b = os.path.basename(e).lower()
    return e == os.path.realpath(sys.executable) or b.startswith('python')
def _hook(ev, args):
    if ev == 'open':
        path, mode, flags = args
        if isinstance(path, int) or path is None: return
        w = (isinstance(mode, str) and any(c in mode for c in 'wax+')) or (isinstance(flags, int) and flags & _WFLAGS)
        if w:
            if not _inside(path, _WRITE_OK): _deny('write', path)
        elif not _inside(path, _READ_OK): _log_read(path)
    elif ev in ('os.rename', 'os.replace', 'shutil.move', 'shutil.copyfile', 'shutil.copytree', 'os.link', 'os.symlink'):
        dst = args[1]
        if not _inside(dst, _WRITE_OK): _deny(ev, dst)
        if ev in ('os.rename', 'os.replace', 'shutil.move') and not _inside(args[0], _WRITE_OK): _deny(ev, args[0])
    elif ev in ('os.remove', 'os.rmdir', 'os.mkdir', 'os.truncate', 'os.chmod', 'os.chown', 'shutil.rmtree', 'os.utime', 'os.chflags', 'os.lchmod', 'os.lchflags', 'os.removexattr', 'os.setxattr'):
        if not _inside(args[0], _WRITE_OK): _deny(ev, args[0])
    elif ev == 'os.chdir':
        if not _inside(args[0], [KIT]): _deny('chdir', args[0])
    elif ev == 'subprocess.Popen':
        exe, argv, cwd, env = args
        if cwd is not None and not _inside(cwd, [KIT]): _deny('child working directory', cwd)
        if cwd is None and not _inside(os.getcwd(), [KIT]): _deny('child working directory', os.getcwd())
        a = list(argv) if isinstance(argv, (list, tuple)) else [argv]
        prog = exe if exe is not None else (a[0] if a else None)
        if _is_py(prog) and not any(isinstance(x, (str, os.PathLike)) and _abs(x) == os.path.realpath(GUARD) for x in a[1:]):
            _deny('python child without kit_guard.py', a[:4])
        if not _is_py(prog): _deny('non-Python child process', a[:4])
    elif ev in ('os.system', 'os.spawn', 'os.startfile'):
        _deny(ev, args[:2] if isinstance(args, tuple) else args)
# route Python children through the guard
_POPEN_INIT = subprocess.Popen.__init__
def _flags_end(a):
    i = 1
    while i < len(a):
        x = a[i]
        if x in ('-X', '-W'): i += 2; continue
        if isinstance(x, str) and x.startswith('-') and x not in ('-c', '-m', '-') and len(x) > 1: i += 1; continue
        break
    return i
def _popen_init(self, args, *a, **kw):
    if isinstance(args, (list, tuple)) and args and _is_py(kw.get('executable') or args[0]):
        args = list(args); i = _flags_end(args)
        if not (i < len(args) and _abs(args[i]) == os.path.realpath(GUARD)): args.insert(i, GUARD)
    return _POPEN_INIT(self, args, *a, **kw)
subprocess.Popen.__init__ = _popen_init
sys.addaudithook(_hook)

if __name__ == '__main__':
    av = sys.argv[1:]
    if not av: sys.exit('usage: kit_guard.py <script.py> [args] | -c <code> | -m <module>')
    if av[0] == '-c':
        sys.argv = ['-c'] + av[2:]
        if not (sys.flags.isolated or getattr(sys.flags, 'safe_path', False)): sys.path.insert(0, '')
        exec(compile(av[1], '<string>', 'exec'), {'__name__': '__main__'})
    elif av[0] == '-m':
        sys.argv = [av[1]] + av[2:]
        if not (sys.flags.isolated or getattr(sys.flags, 'safe_path', False)): sys.path.insert(0, os.getcwd())
        runpy.run_module(av[1], run_name='__main__', alter_sys=True)
    else:
        target = os.path.abspath(av[0]); sys.argv = [av[0]] + av[1:]
        if sys.path and os.path.realpath(sys.path[0]) == os.path.realpath(os.path.dirname(GUARD)) and not (sys.flags.isolated or getattr(sys.flags, 'safe_path', False)):
            sys.path[0] = os.path.dirname(target)                           # as `python script.py`: the script's own directory first
        elif not (sys.flags.isolated or getattr(sys.flags, 'safe_path', False)):
            sys.path.insert(0, os.path.dirname(target))
        runpy.run_path(target, run_name='__main__')
