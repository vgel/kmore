use std::io::BufWriter;
use std::time::Instant;
use std::{fs::File, io::BufReader, path::PathBuf};

use quickdna::{BaseSequence, DnaSequenceAmbiguous, NucleotideAmbiguous};

use kmore::{Wtab, wmer_to_u64, update_wmer_u64};
use quickdna::Nucleotide;

pub fn main() {
    match (
        std::env::args().nth(1).as_deref(),
        std::env::args().nth(2),
        std::env::args().nth(3),
    ) {
        (Some("build-wtab"), Some(haz_path), Some(wtab_path)) => {
            build_and_write_wtab(PathBuf::from(haz_path), PathBuf::from(wtab_path));
        }
        (Some("query"), Some(query_path), Some(wtab_path)) => {
            query_wtab(PathBuf::from(query_path), PathBuf::from(wtab_path));
        }
        _ => {
            eprintln!("USAGE:");
            eprintln!("    build-wtab <haz-path> <wtab-path>");
            eprintln!("    query <query-file-path> <wtab-path>");
            std::process::exit(1);
        }
    }
}

fn build_and_write_wtab(hazard_dir: PathBuf, wtab_path: PathBuf) {
    eprintln!("(building wtab...)");
    let start_wtab = Instant::now();

    let mut wtab = Wtab::new();
    let parser = quickdna::FastaParser::<DnaSequenceAmbiguous>::default();

    for de in walkdir::WalkDir::new(hazard_dir).into_iter() {
        let de = de.unwrap();
        if de.file_type().is_file() {
            let reader = BufReader::new(File::open(de.path()).unwrap());
            let fasta = parser.parse(reader).unwrap();

            for record in fasta {
                for window in amb_nucleotides_to_unamb_kmers::<18>(record.contents.as_slice()) {
                    wtab.add_seed(&window);
                }
            }
        }
    }

    eprintln!(
        "(building wtab took {:.2}s, {:.2}% filled)",
        Instant::now().duration_since(start_wtab).as_secs_f64(),
        wtab.fill_rate() * 100.,
    );

    wtab.write(BufWriter::new(File::create(wtab_path).unwrap()))
        .unwrap();
}

fn query_wtab(query_file_path: PathBuf, wtab_path: PathBuf) {
    let wtab = Wtab::read(BufReader::new(File::open(wtab_path).unwrap())).unwrap();

    let parser = quickdna::FastaParser::<DnaSequenceAmbiguous>::default();
    let parsed = parser
        .parse(BufReader::new(File::open(query_file_path).unwrap()))
        .unwrap();

    eprintln!("(start querying...)");
    let start_query = Instant::now();

    let mut n_windows = 0_u64;
    let mut n_matched = 0_u64;
    for record in parsed.records {
        for window in amb_nucleotides_to_unamb_kmers::<42>(record.contents.as_slice()) {
            n_windows += 1;

            let mut w_matches = 0_u8;
            let mut i = 0_usize;

            let wmer_size = 18;
            assert!(wmer_size * 2 < 42);
            let mut last: Option<(usize, u64)> = None;
            while i < (42 - wmer_size) {
                // optimization: shift a new nucleotide onto the end of the key if possible instead of recomputing it
                let query = match last {
                    Some((last_i, key)) if i == last_i + 1 => {
                        let key = update_wmer_u64(key, window[i + wmer_size - 1]);
                        last = Some((i, key));
                        wtab.query_idx(key)
                    }
                    _ => {
                        let key = wmer_to_u64(&window[i..i + wmer_size].try_into().unwrap());
                        last = Some((i, key));
                        wtab.query_idx(key)
                    }
                };
                // if wtab.query(&window[i..i + wmer_size].try_into().unwrap()) {
                if query {
                    w_matches += 1;
                    if w_matches >= 2 {
                        n_matched += 1;
                        break;
                    } else if i > 42 - (wmer_size * 2) {
                        // can't fit a second wmer after this one, bail
                        break;
                    } else {
                        i += wmer_size;
                    }
                } else {
                    i += 1;
                }
            }
        }
    }
    eprintln!(
        "(querying took {:.2}s)",
        Instant::now().duration_since(start_query).as_secs_f64()
    );

    println!(
        "matched {n_matched} / {n_windows} windows ({}%)",
        (n_matched as f64) / (n_windows as f64) * 100.
    );
}

fn amb_nucleotides_to_unamb_kmers<const LEN: usize>(
    amb_nucs: &[NucleotideAmbiguous],
) -> impl Iterator<Item = [Nucleotide; LEN]> + '_ {
    amb_nucs
        .windows(LEN)
        .filter_map(|w| -> Option<[[Nucleotide; LEN]; 2]> {
            let mut fwd = [Nucleotide::A; LEN];
            let mut rc = [Nucleotide::A; LEN];
            for (i, amb_n) in w.iter().enumerate() {
                match amb_n {
                    NucleotideAmbiguous::A => {
                        fwd[i] = Nucleotide::A;
                        rc[LEN - i - 1] = Nucleotide::T;
                    }
                    NucleotideAmbiguous::T => {
                        fwd[i] = Nucleotide::T;
                        rc[LEN - i - 1] = Nucleotide::A;
                    }
                    NucleotideAmbiguous::C => {
                        fwd[i] = Nucleotide::C;
                        rc[LEN - i - 1] = Nucleotide::G;
                    }
                    NucleotideAmbiguous::G => {
                        fwd[i] = Nucleotide::G;
                        rc[LEN - i - 1] = Nucleotide::C;
                    }
                    _ => return None, // ambiguity code = skip window
                }
            }
            Some([fwd, rc])
        })
        .flatten()
}
