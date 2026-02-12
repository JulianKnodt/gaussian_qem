# Gaussian (Splatting) QEM

QEM implementation for simplifying Gaussian Splats.

# Installation

1. Install [Rust](https://rust-lang.org/tools/install/).

Make sure to use the nightly toolchain.

2. clone this repository

```
git clone git@github.com:JulianKnodt/gaussian_qem.git
```
or
```
git clone https://github.com/JulianKnodt/gaussian_qem.git
```

and then `cd gaussian_qem`.

3. Build/Run the executable

```
cargo build --release
./target/release/decimate -h
```
or
```
cargo run --release -- -h
```

- Example usage:

```
./target/release/decimate -i my_splat.ply -o out.ply -r 0.5
```
