# downloads gsi_sim that is appropriate for the operating system

If the system is Mac, or Windows, this function will download a
precompiled binary from GitHub. In other cases, or if fromSource ==
TRUE, this function will attempt to download the source code and compile
the program from source and install it.

## Usage

``` r
install_gsi_sim(
  commit = "080f462c8eff035fa3e9f2fdce26c3ac013e208a",
  fromSource = FALSE
)
```

## Arguments

- commit:

  The full SHA-1 hash from GitHub from which to get the binary or source

- fromSource:

  If TRUE, download source, even if a binary is available. If FALSE,
  then it will download a precompiled binary, if available. If a binary
  is not available, then it will attempt to download the source.

## Details

If this function fails, then you can just compile gsi_sim by going to
GITHUB_URL and compiling it yourself and naming the executable gsi_sim
and putting it at the location specified by the function
[`gsi_sim_binary_path`](https://thierrygosselin.github.io/assigner/reference/gsi_sim_binary_path.md).
