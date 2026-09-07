# C-style parsing code refactor (libmbutil, libmbc)

Analysis and modernization of the parsing code in `libraries/libmbutil` and
`libraries/libmbc`, targeting error sources linked to C-style char scanning
and copying functions and fixed-length `char*` buffers, converted to modern
C++ (`std::string`) where applicable.

## Sources of error identified

### libmbutil (C++ parsing code)

1. **`LowParser::sCurrWordBuf`** (`parser.cc`) — a manually doubled `char*`
   buffer grown with `SAFENEWARR`/`memcpy`. On EOF in mid-word, `PackWords()`
   returned *without null-terminating*, so `sGetWord()` could return stale
   content from a previous, longer word.
2. **`HighParser::sStringBuf[8192]` / `sStringBufWithSpaces[8192]`**
   (`parser.h`) — fixed buffers behind `GetString()`, `GetStringWithDelims()`
   and `ParseWord()`:
   - `GetString()` and `ParseWord()` *silently truncated* input beyond
     8 KB (`BUFSIZ`);
   - `GetStringWithDelims()` threw a hard "End-of-buffer" error for long
     strings;
   - the `,`/`;` "null string" path of `GetStringWithDelims()` returned the
     *unmodified previous buffer contents* (stale garbage) instead of the
     documented null result.
3. **`IncludeParser`** (`parsinc.cc`) — `sCurrPath`/`sInitialPath`/`sCurrFile`
   were raw `char*` managed with `SAFESTRDUP`/`SAFEDELETEARR`, with
   `sInitialPath` *aliasing* the first `sCurrPath` allocation: after
   `Close()`, `sInitialPath` dangled; in non-`USE_INCLUDE_PARSER` builds
   `GetLineData()` called `strcmp(NULL, NULL)` (undefined behavior).
   `GetFileName()` failed outright for resolved paths ≥ 8192 characters.
4. **`expand_environment()` / `resolve_filename()`** (`parsinc.cc`) — fixed
   `subst[10]` substitution array (hard failure on more than 8 `$VAR`
   substitutions in a file name) plus hand-computed `memcpy` length
   arithmetic; `resolve_filename()` used `strcpy` concatenation and a fragile
   pointer-range test to decide buffer ownership.
5. **`MathParser`** (`mathp.cc`):
   - the number-token scanner used `char s[BUFSIZ]` with an artificial
     "value too long" failure;
   - `readplugin()` (both the expression-evaluator and the `DO_NOT_USE_EE`
     copies) accumulated arguments in `char buf[BUFSIZ]` (hard
     "buffer overflow" error) and leaked the strdup'ed `argv` strings on
     error paths;
   - the `DO_NOT_USE_EE` copy of `readplugin()` referenced an undeclared
     variable `c` instead of `cIn` — that build variant **did not compile**
     (verified against the pristine tree);
   - `trim_arg()` read past index 0 on an all-whitespace string.

### libmbc (C coupling API)

`mbc.c`/`sock.c` are part of the public **C** API, so `std::string` is not
applicable; issues were fixed in C idiom instead:

6. **`winsock_err_string(int, char*)`** (`sock.c`) — `strcpy` into an unsized
   caller buffer with **no `default` case**: an unknown error code left the
   caller's `char msg[100]` uninitialized, which was then printed.
7. **`sock_err_string()`** (`sock.c`) — guarded by `#ifdef _win32`
   (lowercase, never defined), so Windows got `strerror()` for winsock error
   codes; the intended branch returned a pointer to a **stack-local buffer**
   (dangling pointer).

## Changes

- `parser.h` / `parser.cc`: `LowParser` word buffer and both `HighParser`
  string buffers are now `std::string` — no length limits, no silent
  truncation, always terminated. Public `const char*` return types are
  preserved (pointers into the member strings), so none of the external call
  sites change.
- `parsinc.h` / `parsinc.cc`: path/file members and `MyInput` use
  `std::string` (fixes the `sInitialPath` aliasing dangle);
  `expand_environment()`/`resolve_filename()` rewritten around `std::string`
  with no substitution limit; `GetFileName()` no longer has a size cap.
- `mathp.h` / `mathp.cc`: number scanner and both `readplugin()` copies use
  `std::string`/`std::vector<std::string>` (a small `std::vector<char*>` shim
  preserves the legacy `PlugIn::Read(argc, argv)` interface); `trim_arg()`
  takes `std::string&`; the `c`/`cIn` bug and the error-path leaks are gone.
- `sock.h` / `sock.c`, `mbc.c` and call sites (`usesock.cc`, `dataman6.cc`,
  `socketstreamdrive.cc`, `strext.cc`): `sock_err_string()` and
  `winsock_err_string()` now return `const char*` string literals (with an
  "Unknown winsock error" default) — no copying, no caller buffers.

### Deliberate behavior change

`GetStringWithDelims()` on an empty argument (`,`/`;`) now returns an empty
string instead of stale buffer contents.

## Verification

- Every modified file syntax-checked with `g++ -std=c++20` / `gcc` against a
  stub `mbconfig.h` (unconfigured tree), including the `-DDEBUG` and
  `-DDO_NOT_USE_EE` variants; the pristine tree compiles identically except
  for the pre-existing `DO_NOT_USE_EE` breakage, which is now fixed.
- An end-to-end test linking the real modified sources exercised: C-style
  comments, hex escapes (`\41`), keywords, ints/reals, a 20 000-character
  quoted string (fatal before), and an `include: "$TESTDIR/sub.mbd"` whose
  `set:` statement is visible in the symbol table afterwards — all pass,
  clean under AddressSanitizer + UBSan.

## Follow-up: remaining C-style string code in libmbutil

A second pass converted the non-parsing C-style string code flagged above.

### Sources of error identified

8. **`FileName`** (`filename.h`/`fn_UNIX.cc`) — raw `sName`/`sExt` buffers
   with manual `strcpy`/reallocation logic and an interior `sRef` pointer:
   - `iMaxSize`/`iCurSize` were never initialized by the constructor when
     built with a `NULL` name (the `OutputHandler(void)` path);
   - on buffer reuse with a shorter name, `iCurSize = iMaxSize -
     strlen(sExt)` mixed the allocation size with the string length;
   - `_sPutExt()` with an empty stored extension returned the buffer with
     the *previous* call's extension still appended (stale content).
9. **`mbdyn_make_salt()`** (`crypt.cc`) — in the default (`rand()`) path the
   34-byte scratch buffer was **never NUL-terminated**, so a `%s`-style
   `salt_format` (as used by `auth.cc`) read past the buffer into
   uninitialized stack memory. The `/dev/random` branch called `fopen()`
   with one argument (**did not compile** when `HAVE_DEV_RANDOM`/`URANDOM`
   was defined — verified), never checked for `fopen`/`fread` failure, and
   indexed the salt charset with a possibly negative `signed char`.
10. **`mbstrbuf`** (`mbstrbuf.h`/`.cc`) — hand-rolled growable char buffer
    with `memcpy`/`strlen` arithmetic, **no destructor** (leak) and default
    copy semantics on a raw pointer. The class has no users in the tree.

### Changes

- `FileName` reimplemented on `std::string` base/extension members; the
  `const char*` public API (`iInit`, `_sPutExt`, `sGet`) is unchanged.
  Behavior verified byte-identical against the old implementation over a
  matrix of names × `iExtSepNum` values (including buffer-reuse sequences),
  except that `_sPutExt()`/`sGet()` with an empty stored extension now
  return the plain base name instead of stale buffer contents — the only
  user (`OutputHandler::Open`) strips the last extension of the returned
  name itself, so it is unaffected.
- `mbstrbuf` reimplemented on `std::string` (no leak, sane copies);
  `get_len()` now returns the string length rather than the allocation
  size (class is unused in-tree).
- `mbdyn_make_salt()` (C API unchanged): scratch buffer always
  NUL-terminated, `fopen` mode fixed with fallback to `rand()` on
  open/read failure, charset indexed through `unsigned char`.

The `SAFESTRDUP` macros in `mynewmem.h` are left as is: they are part of
the memory-manager facility used throughout the tree, and after this pass
no parsing code uses them any more.
