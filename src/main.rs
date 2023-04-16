use rand::prelude::*;
use std::cell::RefCell;
use std::rc::Rc;
use ndarray::{Array3};
use std::time::Instant;

use clap::Parser;

#[derive(Parser, Debug)]
#[clap(version = "0.1", author = "Your Name <your.email@example.com>", about="polmyer code using Rust")]

struct Args {
    #[arg(short, long, default_value_t = false)]
    layers: bool,

    #[arg(short, long, default_value_t = false)]
    summary: bool,
}


struct Molecule {
    moltype: u32,
    preact: f32,
    neighbors: Vec<Rc<RefCell<Molecule>>>,
}

impl Molecule {
    fn new(moltype: u32, preact: f32) -> Self {
        Self {
            moltype,
            preact,
            neighbors: Vec::new(),
        }
    }
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        Self {
            moltype: self.moltype,
            preact: self.preact,
            neighbors: self.neighbors.clone(),
        }
    }
}

const LA: usize = 10; // dimension of cube , maximum is currently 32 for default array init
const N: usize = 2; // maximum number of neighbors per molecule
const PLINK: f64 = 0.5; // probability of linking neighboring molecules
const PINIT: f64 = 0.5; //probability of molecule creation


fn main() {
    let start_time = Instant::now();
    let args = Args::parse();

    let mut grid = create_grid();
    react_neighbors(&mut grid);
    
    if args.layers {
        print_layers(&grid);
    }
    
    if args.summary {
        print_summary(&grid);
    }

    let elapsed_time = start_time.elapsed();

    println!("Total time elapsed: {:?}", elapsed_time);
}


fn create_grid() -> Array3<Option<Rc<RefCell<Molecule>>>> {
    let mut grid: Array3<Option<Rc<RefCell<Molecule>>>> = Array3::from_elem((LA, LA, LA), None);
    let mut rng = rand::thread_rng();

    for i in 0..LA {
        for j in 0..LA {
            for k in 0..LA {
                if rng.gen_bool(PINIT) {
                    let molecule = Rc::new(RefCell::new(Molecule::new(rng.gen_range(0..2), PLINK as f32)));
                    grid[[i, j, k]] = Some(molecule);
                }
            }
        }
    }

    grid
}

fn react_neighbors(grid: &mut Array3<Option<Rc<RefCell<Molecule>>>>) {
    let mut rng = rand::thread_rng();

    for i in 0..LA {
        for j in 0..LA {
            for k in 0..LA {
                if let Some(molecule) = &grid[[i, j, k]] {
                    let mut neighbor_indices = Vec::new();
                    if i > 0 { neighbor_indices.push((i - 1, j, k)); }
                    if i < LA - 1 { neighbor_indices.push((i + 1, j, k)); }
                    if j > 0 { neighbor_indices.push((i, j - 1, k)); }
                    if j < LA - 1 { neighbor_indices.push((i, j + 1, k)); }
                    if k > 0 { neighbor_indices.push((i, j, k - 1)); }
                    if k < LA - 1 { neighbor_indices.push((i, j, k + 1)); }

                    for (ni, nj, nk) in neighbor_indices {
                        if let Some(neighbor) = &grid[[ni, nj, nk]] {
                            let mut molecule_borrowed = molecule.borrow_mut();
                            if rng.gen_bool(molecule_borrowed.preact as f64) {
                                let mut neighbor_borrowed = neighbor.borrow_mut();

                                if molecule_borrowed.neighbors.len() < N && neighbor_borrowed.neighbors.len() < N {
                                    molecule_borrowed.neighbors.push(neighbor.clone());
                                    neighbor_borrowed.neighbors.push(molecule.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



fn print_summary(grid: &Array3<Option<Rc<RefCell<Molecule>>>>) {
    let mut total_molecules = 0;
    let mut total_neighbors = 0;
    let mut neighbor_distribution: Vec<usize> = vec![0; N + 1];

    for (_index, cell) in grid.indexed_iter() {
        if let Some(molecule) = cell {
            total_molecules += 1;
            let molecule_borrowed = molecule.borrow();
            let neighbor_count = molecule_borrowed.neighbors.len();
            total_neighbors += neighbor_count;
            neighbor_distribution[neighbor_count] += 1;
        }
    }

    let average_neighbors = total_neighbors as f64 / total_molecules as f64;

    println!("Summary");
    println!("============");
    println!("Grid size: {}", LA);
    println!("Total molecules: {}", total_molecules);
    println!("Average neighbors per molecule: {:.2}", average_neighbors);
    println!("Neighbor distribution:");
    for (count, occurrences) in neighbor_distribution.iter().enumerate() {
        println!("  Molecules with {} neighbors: {}", count, occurrences);
    }
}


fn print_layers(grid: &Array3<Option<Rc<RefCell<Molecule>>>>) {
    for k in 0..LA {
        println!("Layer {}", k + 1);
        println!("========");
        for i in 0..LA {
            for j in 0..LA {
                if let Some(molecule) = &grid[[i, j, k]] {
                    let molecule_borrowed = molecule.borrow();

                    let mut cell_repr = String::from(" ");
                    if k > 0 && molecule_borrowed.neighbors.iter().any(|n| grid[[i, j, k - 1]].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
                        cell_repr.push('V');
                    }
                    if k < LA - 1 && molecule_borrowed.neighbors.iter().any(|n| grid[[i, j, k + 1]].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
                        cell_repr.push('^');
                    }

                    print!("{:2}", cell_repr);
                    if j < LA - 1 && molecule_borrowed.neighbors.iter().any(|n| grid[[i, j + 1, k]].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
                        print!("-");
                    } else {
                        print!(" ");
                    }
                } else {
                    print!("   ");
                }
            }
            println!();

            if i < LA - 1 {
                for j in 0..LA {
                    if let Some(molecule) = &grid[[i, j, k]] {
                        let molecule_borrowed = molecule.borrow();
                        if molecule_borrowed.neighbors.iter().any(|n| grid[[i + 1, j, k]].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
                            print!(" | ");
                        } else {
                            print!("   ");
                        }
                    } else {
                        print!("   ");
                    }
                }
                println!();
            }
        }
        println!();
    }
}
