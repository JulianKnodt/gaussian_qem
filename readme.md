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

## Example:

<table>
<tr>
<td>
<figure>
<img src="assets/heart_cookie.png" width="240" alt="Ground Truth"/>
<br/>
<figcaption align="left">Original</figcaption>
</figure>
</td>

<td>
<figure>
<img src="assets/heart_cookie_half.png" width="240" alt="Half Resolution"/>
<br/>
<figcaption align="right">50% Reduced</figcaption>
</figure>
</td>
</tr>
</table>

On the left is the original, the simplified with half of the splats is shown on the right. The
color appears significantly whiter, since spherical harmonics are not well-preserved.

This can be fixed with an additional optimization step, using a modified training step
[here](https://github.com/julianknodt/gaussian-splatting).
