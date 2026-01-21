# rs-gax
Rust library to parse gaf/gam/gamp files.

### TODO
- [X] GAF I/O
- [X] GAM I/O
- [X] GAMP I/O
- [ ] documentation
- [ ] tests

## Example

Reading a gam file:

```rust
use rs_gax::{gam, vg};
use std::fs::File;

fn main(){
    let in_file = "example.gam";
    let f = File::open(in_file).unwrap();
    let gam: Vec<vg::Alignment> = gam::parse(f).unwrap();
    println!("{:?}", gam);
}
```

Writing a gam file:

```rust
use rs_gax::{gam, vg};
use std::fs::File;

fn main(){
    let out_file = "example.out.gam";
    let f = File::create(out_file).unwrap();
    let mut alignment = vg::Alignment::default();
    alignment.name = "Test".to_string();

    let gam: Vec<vg::Alignment> = vec![alignment];
    gam::write(&gam, f).unwrap();
}
```

Reading a gamp file:

```rust
use rs_gax::{gamp, vg};
use std::fs::File;

fn main(){
    let in_file = "example.gamp";
    let f = File::open(in_file).unwrap();
    let gamp: Vec<vg::MultipathAlignment> = gamp::parse(f).unwrap();
    println!("{:?}", gamp);
}
```

Writing a gamp file:

```rust
use rs_gax::{gamp, vg};
use std::fs::File;

fn main(){
    let out_file = "example.out.gamp";
    let f = File::create(out_file).unwrap();
    let mut alignment = vg::MultipathAlignment::default();
    alignment.name = "Test".to_string();

    let gamp: Vec<vg::MultipathAlignment> = vec![alignment];
    gamp::write(&gamp, f).unwrap();
}
```
