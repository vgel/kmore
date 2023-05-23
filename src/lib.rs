use std::io;

use bitvec::prelude::*;

use byte_slice_cast::{AsByteSlice, AsMutByteSlice};
use quickdna::Nucleotide;

pub fn wmer_to_u64(kmer: &[Nucleotide; 18]) -> u64 {
    let mut u = 0_u64;
    for n in kmer {
        u <<= 2;
        u |= match n {
            Nucleotide::A => 0b00,
            Nucleotide::T => 0b01,
            Nucleotide::C => 0b10,
            Nucleotide::G => 0b11,
        };
    }
    u
}

pub fn update_wmer_u64(mut prev: u64, next_n: Nucleotide) -> u64 {
    prev <<= 2;
    prev |= match next_n {
        Nucleotide::A => 0b00,
        Nucleotide::T => 0b01,
        Nucleotide::C => 0b10,
        Nucleotide::G => 0b11,
    };
    prev &= 0b111111111111111111111111111111111111;
    prev
}

pub struct Wtab {
    seeds: BitBox<u64>,
}

impl Wtab {
    pub fn new() -> Self {
        Self {
            seeds: bitbox![u64, Lsb0; 0; 2_usize.pow(36)],
        }
    }

    pub fn add_seed(&mut self, seed: &[Nucleotide; 18]) {
        let u = wmer_to_u64(seed);
        self.seeds.set(u as usize, true);
    }

    pub fn fill_rate(&self) -> f64 {
        (self.seeds.count_ones() as f64) / (u32::MAX as f64)
    }

    pub fn query_idx(&self, idx: u64) -> bool {
        *self.seeds.get(idx as usize).unwrap()
    }

    pub fn query(&self, kmer: &[Nucleotide; 18]) -> bool {
        self.query_idx(wmer_to_u64(kmer))
    }

    pub fn write(&self, mut w: impl io::Write) -> io::Result<()> {
        w.write_all(self.seeds.as_raw_slice().as_byte_slice())
    }

    pub fn read(mut r: impl io::Read) -> io::Result<Self> {
        let mut this = Self::new();
        let slice = this.seeds.as_raw_mut_slice().as_mut_byte_slice();
        r.read_exact(slice)?;
        Ok(this)
    }
}

impl Default for Wtab {
    fn default() -> Self {
        Self::new()
    }
}
