# claude_comments

Annotated copies of the C files Veyrat added or changed in the `subdivide` branch of ORIGAMI (https://github.com/dveyrat/origami/tree/subdivide, 2022). This is the version used in Veyrat, Douglass & BenZvi 2024 to run ORIGAMI on the 2304³-particle AbacusSummit Hugebase simulation.

The code is the same as the branch; only comments were added, all marked `[CLAUDE]`. In `origamidiv.c`, lines Veyrat added or changed relative to upstream `origamitag.c` end with `[VEYRAT]`. Comments marked `[ISSUE]` point out problems found while annotating.

| File | Status in branch | What it does |
| --- | --- | --- |
| `origamidiv.c` | new | Computes morphology for one sub-cube plus a buffer layer, then saves only the sub-cube |
| `origamicomb.c` | new | Stitches all sub-cube tag files into one full-box `tag.dat` |
| `readfiles.c` | modified | Adds the `nbuf` parameter and `divread()`, which reads one sub-cube plus buffer from the position file |

The branch also adds one line (`extern int nbuf;`) to `params.h`, an `nbuf` section to `parameters_example.txt`, and an `origamidiv` target to the Makefile. There is no Makefile target for `origamicomb`.

## Workflow

```
for divid in 0 .. ndiv^3-1:          # independent; can run in parallel
    ./origamidiv params.txt ndiv divid    ->  <label>_tag_<divid>.dat
rename each to <label>tag_<x><y><z>.dat   (x = divid%ndiv, y = (divid/ndiv)%ndiv, z = divid/ndiv²)
./origamicomb <label> <outdir>/ np1d ndiv ->  <label>tag.dat   (same format origamihalo reads)
```

## Issues found while annotating

1. **`origamicomb.c` does not compile.** `fwrite(&np3,1,sizeof(int));` is missing its file argument; it should be `fwrite(&np3,1,sizeof(int),tag);`.
2. **The file names don't match.** `origamidiv` writes `<label>_tag_<divid>.dat`, but `origamicomb` reads `<label>tag_<x><y><z>.dat`, so the files need renaming in between.
3. **The position file needs an 8-byte particle count.** `divread` reads the count as an 8-byte `long`, while upstream `posread` reads a 4-byte `int`, so one file can't be used by both readers.
4. **Use `numfiles = 0`.** The Gadget reader (`numfiles > 0`) loads the whole box in the wrong layout for `origamidiv`.
5. **Divisibility requirements:**
   - `np1d` must be divisible by `ndiv`³ for `origamicomb` (2304 works with 2 or 4, but not 3).
   - `npd + 2·nbuf` must be divisible by `nsplit`, or some rows are never processed.
6. **The header count overflows for large boxes.** In `origamicomb`, `np3 = np1d³` is a 4-byte int, which overflows above about 1290³. Downstream `origamihalo` and `origamienv` also read the count as an int.
7. **Re-running appends.** `origamicomb` opens its output in append mode, so a second run adds to an existing `tag.dat` instead of replacing it.
8. **Every parameter file needs the `nbuf` section.** The changed `readfiles.c` is shared by every program, so all of them now require it, in the position shown in `readfiles.c`.
9. **`readfiles.c` fails on recent compilers (this is upstream, not Veyrat's).** It calls `load_snapshot`, `reordering` and `allocate_memory` before declaring them. Recent clang (and gcc 14+) treat that as an error, so add `-Wno-error=implicit-function-declaration` to `CFLAGS` (or declare the functions near the top of the file).
10. **The search distance changes.** In `origamidiv`, crossings are searched up to `nbuf` grid steps, while upstream searches `np1d/4`. Results only match upstream if `nbuf ≥ np1d/4`; a smaller `nbuf` can miss distant crossings.

## C primer (just enough to read these files)

- **Statements** end with `;`. **Blocks** are wrapped in `{ }`. Comments are `/* ... */` or `// ...`.
- **Types:** `int` is a whole number (4 bytes), `long` is a big whole number (8 bytes), `float` is a decimal number (4 bytes), `unsigned char` is one byte (0–255, used here for tags), and `char[80]` is text.
- **Pointers:** `*` declares a pointer, which is the address of some data. In practice here, `float *a` is "a list of floats" and `float **r` is "a list of lists" (`r[i][0..2]` = x, y, z of particle `i`). `&x` means "the address of x", used so a function can write into `x`.
- **Memory:** `malloc(n)` reserves `n` bytes and returns their address (`NULL` if it fails). `free` gives the memory back. `sizeof(t)` is the size of type `t` in bytes.
- **Arithmetic:** `/` between whole numbers drops the remainder (7/2 = 3). `%` gives the remainder (7 % 2 = 1). `a *= 2` means `a = a * 2`, and `++j` means `j = j + 1`.
- **Logic:** `==` means equals, `!=` not equal, `||` or, `&&` and. A true comparison counts as 1 and a false one as 0, so values can be added up.
- **Loops:**
  - `for (i=0; i<n; i++) {...}` runs the block for i = 0, 1, …, n−1.
  - `while (cond) {...}` repeats the block while `cond` is true.
  - `continue` skips to the next loop pass, and `break` leaves the loop.
- **Files:**
  - `fopen(name, "r"/"w"/"a")` opens a file to read, write (overwrite) or append.
  - `fread`/`fwrite(ptr, size, count, file)` copy raw bytes between memory and the file. The file remembers how far it has been read.
- **Text:** `printf` prints, `sprintf` builds a string (for example, a filename), and `sscanf`/`fscanf` read numbers out of text. In the format strings, `%d` is an int, `%ld` a long, `%f` a float and `%s` a word.
- **`#define`** makes a text shortcut that is substituted before compiling. **`#include`** pulls in another file's declarations.
- **`#pragma omp parallel for`** spreads a loop's iterations across CPU cores when compiled with `-fopenmp`.
- **`main(argc, argv)`** is where the program starts. `argv[1]`, `argv[2]`, … are the command-line words, and `argc` counts them, including the program name.
