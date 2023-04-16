use rand::prelude::*;
use std::cell::RefCell;
use std::rc::Rc;

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

const LA: usize = 32; // dimension of cube , maximum due to rust limits is currently 32
const N: usize = 2; // maximum number of neighbors per molecule
const PLINK: f32 = 0.5; // probability of linking neighboring molecules
const PINIT: f64 = 0.99; //probability of molecule creation

// ... (previous code) ...

fn print_grid_summary(grid: &[[[Option<Rc<RefCell<Molecule>>>; LA]; LA]; LA]) {
    let mut total_molecules = 0;
    let mut total_neighbors = 0;
    let mut neighbor_distribution: Vec<usize> = vec![0; N + 1];

    for i in 0..LA {
        for j in 0..LA {
            for k in 0..LA {
                if let Some(molecule) = &grid[i][j][k] {
                    total_molecules += 1;
                    let molecule_borrowed = molecule.borrow();
                    let neighbor_count = molecule_borrowed.neighbors.len();
                    total_neighbors += neighbor_count;
                    neighbor_distribution[neighbor_count] += 1;
                }
            }
        }
    }

    let average_neighbors = total_neighbors as f64 / total_molecules as f64;

    println!("Grid Summary");
    println!("============");
    println!("Total molecules: {}", total_molecules);
    println!("Average neighbors per molecule: {:.2}", average_neighbors);
    println!("Neighbor distribution:");
    for (count, occurrences) in neighbor_distribution.iter().enumerate() {
        println!("  Molecules with {} neighbors: {}", count, occurrences);
    }
}

// ... (previous code) ...
fn print_grid_layers(grid: &[[[Option<Rc<RefCell<Molecule>>>; LA]; LA]; LA]) {
    for k in 0..LA {
        println!("Layer {}", k + 1);
        println!("========");
        for i in 0..LA {
            for j in 0..LA {
                if let Some(molecule) = &grid[i][j][k] {
                    let molecule_borrowed = molecule.borrow();
                    print!("{:2}", molecule_borrowed.moltype);
                    if j < LA-1 && molecule_borrowed.neighbors.iter().any(|n| grid[i][j + 1][k].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
                        print!("-");
                    } else {
                        print!(" ");
                    }
                } else {
                    print!("   ");
                }
            }
            println!();

            if i < LA-1 {
                for j in 0..LA {
                    if let Some(molecule) = &grid[i][j][k] {
                        let molecule_borrowed = molecule.borrow();
                        if molecule_borrowed.neighbors.iter().any(|n| grid[i + 1][j][k].as_ref().map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))) {
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



fn main() {
    // Create a 3-dimensional grid with dimensions 10 x 10 x 10
    let mut grid: [[[Option<Rc<RefCell<Molecule>>>; LA]; LA]; LA] = Default::default();

    // Fill the grid with random molecules
    let mut rng = rand::thread_rng();
    for i in 0..LA {
        for j in 0..LA {
            for k in 0..LA {
                if rng.gen_bool(PINIT) {
                    let molecule = Rc::new(RefCell::new(Molecule::new(rng.gen_range(0..2),PLINK)));
                    grid[i][j][k] = Some(molecule);
                }
            }
        }
    }

    for i in 0..LA {
        for j in 0..LA {
            for k in 0..LA {
                if let Some(molecule) = grid[i][j][k].clone() {
                    let mut neighbor_indices = Vec::new();
                    if i > 0 { neighbor_indices.push((i - 1, j, k)); }
                    if i < LA-1 { neighbor_indices.push((i + 1, j, k)); }
                    if j > 0 { neighbor_indices.push((i, j - 1, k)); }
                    if j < LA-1 { neighbor_indices.push((i, j + 1, k)); }
                    if k > 0 { neighbor_indices.push((i, j, k - 1)); }
                    if k < LA-1 { neighbor_indices.push((i, j, k + 1)); }

                    for (ni, nj, nk) in neighbor_indices {
                        if let Some(neighbor) = grid[ni][nj][nk].clone() {
                            if rng.gen_bool(PLINK as f64) {
                                let mut molecule_borrowed = molecule.borrow_mut();
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
    print_grid_layers(&grid);
    print_grid_summary(&grid);
}
