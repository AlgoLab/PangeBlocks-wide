use gax::{gaf, gam, gamp};

fn main() {
    if let Err(e) = main_() {
        eprintln!("ERROR: {}", e);
    }
}

fn main_() -> Result<(), Box<dyn std::error::Error>> {
    // GAF I/O
    let gaf = gaf::parse_from_file("data/example.gaf")?;
    gaf::write_to_file(&gaf, "data/example.out.gaf")?;
    assert!(gaf == gaf::parse_from_file("data/example.out.gaf")?);
    println!("GAF: {} records", gaf.len());

    // GAM I/O
    let gam = gam::parse_from_file("data/example.gam")?;
    gam::write_to_file(&gam, "data/example.out.gam")?;
    assert!(gam == gam::parse_from_file("data/example.out.gam")?);
    println!("GAM: {} records", gam.len());

    // GAMP I/O
    let gamp = gamp::parse_from_file("data/example.gamp")?;
    gamp::write_to_file(&gamp, "data/example.out.gamp")?;
    assert!(gamp == gamp::parse_from_file("data/example.out.gamp")?);
    println!("GAMP: {} records", gamp.len());
    Ok(())
}
