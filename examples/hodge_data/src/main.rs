use sc_faber::Q;
use sc_faber::hodge::{lambda_integrals, lambda_monomial_count};
use sc_faber::socle::{CacheSession, set_cache_memory_limit};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

const GIB: u64 = 1 << 30;

struct Options {
    first: u32,
    last: u32,
    threads: usize,
    memory_gib: u64,
    directory: PathBuf,
}

impl Options {
    fn parse() -> Result<Option<Self>, String> {
        let mut options = Self {
            first: 2,
            last: 4,
            threads: 24,
            memory_gib: 2,
            directory: PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/hodge-data"),
        };
        let mut arguments = std::env::args().skip(1);
        while let Some(argument) = arguments.next() {
            if argument == "--help" || argument == "-h" {
                println!(
                    "hodge-data [--genus G | --min-genus G --max-genus G] [--threads N] [--memory-gib N] [--data-dir PATH]"
                );
                println!(
                    "Defaults: genera 2 through 4, 24 compute threads, 2 GiB total-process limit."
                );
                println!("Default output directory: {}", options.directory.display());
                println!(
                    "Completed integral values resume from SQLite; unfinished query planning is rebuilt."
                );
                return Ok(None);
            }
            let value = arguments
                .next()
                .ok_or_else(|| format!("missing value for {argument}"))?;
            match argument.as_str() {
                "--genus" => {
                    options.first = value.parse().map_err(|_| "invalid genus")?;
                    options.last = options.first;
                }
                "--min-genus" => {
                    options.first = value.parse().map_err(|_| "invalid minimum genus")?
                }
                "--max-genus" => {
                    options.last = value.parse().map_err(|_| "invalid maximum genus")?
                }
                "--threads" => {
                    options.threads = value.parse().map_err(|_| "invalid thread count")?
                }
                "--memory-gib" => {
                    options.memory_gib = value.parse().map_err(|_| "invalid memory limit")?
                }
                "--data-dir" => options.directory = value.into(),
                _ => return Err(format!("unknown option {argument}")),
            }
        }
        if options.first < 2 || options.last < options.first || options.last > 20 {
            return Err("require 2 <= minimum genus <= maximum genus <= 20".into());
        }
        if options.threads == 0 || options.memory_gib == 0 {
            return Err("threads and memory limit must be positive".into());
        }
        options
            .memory_gib
            .checked_mul(GIB)
            .ok_or("memory limit overflow")?;
        Ok(Some(options))
    }
}

fn write_result(output: &mut impl Write, parts: &[u32], value: &Q) -> Result<(), String> {
    write!(output, "hodge_number[(").map_err(|e| e.to_string())?;
    for (position, index) in parts.iter().enumerate() {
        if position != 0 {
            write!(output, ", ").map_err(|e| e.to_string())?;
        }
        write!(output, "{index}").map_err(|e| e.to_string())?;
    }
    writeln!(output, ")] = {value}").map_err(|e| e.to_string())
}

fn calculate_genus(session: &CacheSession, options: &Options, g: u32) -> Result<(), String> {
    let complete = options.directory.join(format!("hodge_data_g{g}.sage"));
    if complete.is_file() {
        println!("g={g} already_complete file={}", complete.display());
        return Ok(());
    }
    let partial = options.directory.join(format!("hodge_data_g{g}.sage.part"));
    let mut output = BufWriter::new(File::create(&partial).map_err(|e| e.to_string())?);
    writeln!(
        output,
        "# Pure Hodge integrals, exact rational values.\n# g = {g}"
    )
    .map_err(|e| e.to_string())?;
    let expected = lambda_monomial_count(g)?;
    let mut done = 0;
    let started = Instant::now();
    let mut reported = Instant::now();
    println!("GENUS_START g={g} monomials={expected}");
    let total = session.run(|| {
        lambda_integrals(g, |parts, value| {
            write_result(&mut output, parts, value)?;
            done += 1;
            if reported.elapsed().as_secs() >= 5 || done == expected {
                output.flush().map_err(|e| e.to_string())?;
                println!(
                    "PROGRESS g={g} done={done}/{expected} elapsed_s={:.3}",
                    started.elapsed().as_secs_f64()
                );
                reported = Instant::now();
            }
            Ok(())
        })
    })?;
    if total != expected {
        return Err(format!(
            "genus {g} enumeration ended at {total} instead of {expected}"
        ));
    }
    session.flush()?;
    output.flush().map_err(|e| e.to_string())?;
    output.get_ref().sync_all().map_err(|e| e.to_string())?;
    drop(output);
    fs::rename(&partial, &complete).map_err(|e| e.to_string())?;
    println!(
        "GENUS_DONE g={g} count={} elapsed_s={:.3} file={}",
        total,
        started.elapsed().as_secs_f64(),
        complete.display()
    );
    Ok(())
}

fn execute(options: Options) -> Result<(), String> {
    fs::create_dir_all(&options.directory).map_err(|e| e.to_string())?;
    let workspace = options.directory.join("workspace");
    fs::create_dir_all(&workspace).map_err(|e| e.to_string())?;
    // SAFETY: configure the process environment before starting cache or compute threads.
    unsafe {
        std::env::set_var("SC_FABER_THREADS", options.threads.to_string());
        std::env::set_var("SC_FABER_WORKSPACE_DIR", &workspace);
        std::env::remove_var("SC_FABER_CACHE_MEASURE_PULL");
        std::env::remove_var("SC_FABER_CACHE_MEASURE_SOURCE_ONLY");
    }
    set_cache_memory_limit(Some(options.memory_gib * GIB))?;
    let cache = options.directory.join("cache.sqlite");
    let cache_path = cache.to_str().ok_or("cache path must be UTF-8")?;
    let session = CacheSession::open(cache_path)?;
    let started = Instant::now();
    println!(
        "RUN_START genera={}..={} compute_threads={} memory_gib={} cache={}",
        options.first,
        options.last,
        sc_num::parallel::thread_count("SC_FABER_THREADS", options.threads),
        options.memory_gib,
        cache.display()
    );
    let result =
        (options.first..=options.last).try_for_each(|g| calculate_genus(&session, &options, g));
    let close = session.close();
    result?;
    close?;
    println!("RUN_DONE elapsed_s={:.3}", started.elapsed().as_secs_f64());
    Ok(())
}

fn main() -> std::process::ExitCode {
    match Options::parse().and_then(|options| options.map_or(Ok(()), execute)) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("hodge-data: {error}");
            std::process::ExitCode::FAILURE
        }
    }
}
