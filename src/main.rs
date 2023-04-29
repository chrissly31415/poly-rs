mod visualize;
mod utils;
mod statistics;

use rand::prelude::*;
use std::cell::RefCell;
use std::sync::atomic::{AtomicU32, Ordering};

use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use std::time::Instant;

use clap::Parser;
use ndarray::Array3;

use lazy_static::lazy_static;

use petgraph::algo::kosaraju_scc;


#[derive(Parser, Debug)]
#[clap(
    version = "0.1",
    author = "chrissly31415 <your.email@example.com>",
    about = "Experimental code for polmyers using Rust"
)]
struct Args {
    #[arg(short, long, default_value_t = false)]
    layers: bool,

    #[arg(short, long, default_value_t = false)]
    summary: bool,

    #[arg(short, long, default_value_t = false)]
    debug: bool,

    #[arg(short, long, default_value_t = false)]
    visualize: bool,

    #[arg(short, long, default_value_t = false)]
    bfm: bool,
}

const BASE_VECTORS: [(isize, isize, isize); 6] = [
    (2, 0, 0), //2
    (2, 1, 0), //3
    (2, 1, 1), //4
    (2, 2, 1), //5
    (3, 1, 0), //4
    (3, 0, 0), //3
]; 
//const BASE_VECTORS: [(isize, isize, isize); 3] = [(1, 0, 0), (0, 1, 0), (0, 0, 1)];

fn generate_bond_vectors(base_vectors: &[(isize, isize, isize)]) -> HashSet<(isize, isize, isize)> {
    let mut potential_neighbors: HashSet<(isize, isize, isize)> = HashSet::new();

    for (x, y, z) in base_vectors {
        let coords = [
            (x, y, z),
            (x, z, y),
            (y, x, z),
            (y, z, x),
            (z, x, y),
            (z, y, x),
        ];

        for (x, y, z) in coords {
            let signs: Vec<isize> = vec![-1, 1];
            for sign_x in &signs {
                for sign_y in &signs {
                    for sign_z in &signs {
                        let neighbor = (x * sign_x, y * sign_y, z * sign_z);

                        // Skip adding the origin (0, 0, 0) to the potential neighbors
                        if neighbor != (0, 0, 0) {
                            potential_neighbors.insert(neighbor);
                        }
                    }
                }
            }
        }
    }
    println!("Number of bond vectors: {}", potential_neighbors.len());
    potential_neighbors
}

lazy_static! {
    static ref UNIQUE_BOND_VECTORS: Vec<(isize, isize, isize)> =
        generate_bond_vectors(&BASE_VECTORS).into_iter().collect();
}

pub struct Molecule {
    moltype: u32,
    preact: f32,
    position: (usize, usize, usize),
    neighbors: Vec<Rc<RefCell<Molecule>>>,
    id: u32,
}

static NEXT_ID: AtomicU32 = AtomicU32::new(1);

impl Molecule {
    fn new(moltype: u32, preact: f32, position: (usize, usize, usize)) -> Self {
        let id = NEXT_ID.fetch_add(1, Ordering::SeqCst);
        Self {
            moltype,
            preact,
            position,
            neighbors: Vec::new(),
            id,
        }
    }
    fn cube_positions(&self) -> Vec<(usize, usize, usize)> {
        let (i, j, k) = self.position;
        vec![
            (i + 1, j, k),
            (i, j + 1, k),
            (i + 1, j + 1, k),
            (i, j, k + 1),
            (i + 1, j, k + 1),
            (i, j + 1, k + 1),
            (i + 1, j + 1, k + 1),
        ]
    }
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        Self {
            moltype: self.moltype,
            preact: self.preact,
            position: self.position,
            neighbors: self.neighbors.clone(),
            id: self.id,
        }
    }
}

const LA: usize = 40; // dimension of cube , maximum is currently 32 for default array init
const MAX_NEIGHBORS: usize = 3; // maximum number of neighbors per molecule
const PLINK: f64 = 1.0; // probability of linking neighboring molecules

const MAX_TYPE: usize = 2; // Number of molecule types
                           //const CONCS: [f64; MAX_TYPE] = [0.05, 0.1, 0.2]; // starting conc of each molecule type
const CONCS: [f64; 2] = [0.2, 0.2]; // starting conc of each molecule type

fn main() -> Result<(), std::io::Error> {
    let start_time = Instant::now();
    let args = Args::parse();

    let (mut grid, molecules) = create_random_grid(args.bfm, LA / 2);

    let mut allowed_reactions: HashMap<u32, HashSet<u32>> = HashMap::new();
    allowed_reactions.insert(0, [1].iter().cloned().collect());
    allowed_reactions.insert(1, [0].iter().cloned().collect());
    //allowed_reactions.insert(2, [0,1].iter().cloned().collect());

    react_neighbors(&mut grid, &molecules, &allowed_reactions);
    let graph = utils::molecules_to_graph(&molecules);

    visualize::visualize_chains(&grid, &molecules, &graph);
    visualize::visualize_layers(&grid, &molecules, 0.2, true, args.bfm, true);
    random_move_n_times(&mut grid, &molecules, 10000000, args.bfm);
    
    react_neighbors(&mut grid, &molecules, &allowed_reactions);
    let graph = utils::molecules_to_graph(&molecules);

 
    if args.layers {
        statistics::print_layers(&grid, LA);
    }
    if args.visualize {
        //visualize::visualize_layers(&molecules, 0.2, false);
        visualize::visualize_chains(&grid, &molecules, &graph);
        visualize::visualize_layers(&grid, &molecules, 0.2, true, args.bfm, true);
    }

    if args.summary {
        statistics::print_summary(&grid, MAX_TYPE, MAX_NEIGHBORS, LA);
        statistics::print_occupied_grid_points_and_molecules(&grid, LA);
        statistics::print_bond_statistics(&molecules, &UNIQUE_BOND_VECTORS);
        statistics::print_connected_molecules_statistics(&graph,true);
    }
    if args.debug {
        println!("Allowed reactions: {:?}", allowed_reactions);
    }

    let elapsed_time = start_time.elapsed();
    println!("Total time elapsed: {:?}", elapsed_time);
    Ok(())
}

fn create_random_grid(
    use_bfm: bool,
    z_max: usize,
) -> (
    Array3<Option<Rc<RefCell<Molecule>>>>,
    Vec<Rc<RefCell<Molecule>>>,
) {
    let mut grid: Array3<Option<Rc<RefCell<Molecule>>>> = Array3::from_elem((LA, LA, LA), None);
    let mut rng = rand::thread_rng();
    let mut molecules = Vec::new();

    let mut bound: usize = LA - 1;
    if !use_bfm {
        bound = LA;
    }

    for i in 0..bound {
        for j in 0..bound {
            for k in 0..bound {
                if k >= z_max {
                    continue;
                };
                let moltype = rng.gen_range(0..MAX_TYPE);
                if rng.gen::<f64>() < CONCS[moltype] {
                    let molecule = Rc::new(RefCell::new(Molecule::new(
                        moltype as u32,
                        PLINK as f32,
                        (i, j, k),
                    )));

                    if !use_bfm {
                        grid[[i, j, k]] = Some(molecule.clone());
                        molecules.push(molecule);
                    } else {
                        // Check if adjacent cells are empty
                        if grid[[i, j, k]].is_none()
                            && grid[[i + 1, j, k]].is_none()
                            && grid[[i, j + 1, k]].is_none()
                            && grid[[i + 1, j + 1, k]].is_none()
                            && grid[[i, j, k + 1]].is_none()
                            && grid[[i + 1, j, k + 1]].is_none()
                            && grid[[i, j + 1, k + 1]].is_none()
                            && grid[[i + 1, j + 1, k + 1]].is_none()
                        {
                            grid[[i, j, k]] = Some(molecule.clone());
                            grid[[i + 1, j, k]] = Some(molecule.clone());
                            grid[[i, j + 1, k]] = Some(molecule.clone());
                            grid[[i + 1, j + 1, k]] = Some(molecule.clone());
                            grid[[i, j, k + 1]] = Some(molecule.clone());
                            grid[[i + 1, j, k + 1]] = Some(molecule.clone());
                            grid[[i, j + 1, k + 1]] = Some(molecule.clone());
                            grid[[i + 1, j + 1, k + 1]] = Some(molecule.clone());

                            molecules.push(molecule);
                        }
                    }
                }
            }
        }
    }

    (grid, molecules)
}

fn react_neighbors(
    grid: &mut Array3<Option<Rc<RefCell<Molecule>>>>,
    molecules: &[Rc<RefCell<Molecule>>],
    allowed_reactions: &HashMap<u32, HashSet<u32>>,
) {
    let mut rng = rand::thread_rng();
    for molecule in molecules {
        let reactable_neighbors = get_reactable_neighbors(grid, molecule, allowed_reactions);

        for neighbor in reactable_neighbors {
            if rng.gen_bool(molecule.borrow().preact as f64) {
                let mut molecule_borrowed_mut = molecule.borrow_mut();
                let mut neighbor_borrowed_mut = neighbor.borrow_mut();

                if molecule_borrowed_mut.neighbors.len() < MAX_NEIGHBORS
                    && neighbor_borrowed_mut.neighbors.len() < MAX_NEIGHBORS
                {
                    molecule_borrowed_mut.neighbors.push(neighbor.clone());
                    neighbor_borrowed_mut.neighbors.push(molecule.clone());
                }
            }
        }
    }
}

fn get_reactable_neighbors(
    grid: &Array3<Option<Rc<RefCell<Molecule>>>>,
    molecule: &Rc<RefCell<Molecule>>,
    allowed_reactions: &HashMap<u32, HashSet<u32>>,
) -> Vec<Rc<RefCell<Molecule>>> {
    let molecule_borrowed = molecule.borrow();
    let (i, j, k) = (
        molecule_borrowed.position.0 as isize,
        molecule_borrowed.position.1 as isize,
        molecule_borrowed.position.2 as isize,
    );

    let mut reactable_neighbors = Vec::new();

    for (dx, dy, dz) in UNIQUE_BOND_VECTORS.iter() {
        let (ni, nj, nk) = (
            i.saturating_add(*dx),
            j.saturating_add(*dy),
            k.saturating_add(*dz),
        );

        if ni >= 0
            && nj >= 0
            && nk >= 0
            && (ni as usize) < LA
            && (nj as usize) < LA
            && (nk as usize) < LA
        {
            if let Some(neighbor) = &grid[[(ni as usize), (nj as usize), (nk as usize)]] {
                if Rc::ptr_eq(molecule, neighbor) {
                    continue; // skip if the neighbor is the same as the molecule
                }
                let neighbor_borrowed = neighbor.borrow();
                if let Some(reaction_partners) = allowed_reactions.get(&molecule_borrowed.moltype) {
                    if reaction_partners.contains(&neighbor_borrowed.moltype) {
                        reactable_neighbors.push(neighbor.clone());
                    }
                }
            }
        }
    }
    reactable_neighbors
}

fn random_move_n_times(
    grid: &mut Array3<Option<Rc<RefCell<Molecule>>>>,
    molecules: &Vec<Rc<RefCell<Molecule>>>,
    n: usize,
    use_bfm: bool,
) {
    let mut rng = rand::thread_rng();

    let mut success_rate: i32 = 0;
    for _ in 0..n {
        let random_index = rng.gen_range(0..molecules.len());
        let molecule = &molecules[random_index];
        success_rate += move_molecule_randomly(molecule, grid, use_bfm) as i32;
    }
    let acceptance_ratio = success_rate as f32 / n as f32;
    println!(
        "{} random moves. Acceptance rate {:.3}",
        n, acceptance_ratio
    );
}

fn move_molecule_randomly(
    molecule: &Rc<RefCell<Molecule>>,
    grid: &mut Array3<Option<Rc<RefCell<Molecule>>>>,
    use_bfm: bool,
) -> bool {
    let mut rng = rand::thread_rng();
    let (i, j, k) = {
        let molecule_borrowed = molecule.borrow();
        molecule_borrowed.position
    };

    let directions = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ];

    let (di, dj, dk) = directions.choose(&mut rng).unwrap();

    let new_i = i.wrapping_add(*di as usize);
    let new_j = j.wrapping_add(*dj as usize);
    let new_k = k.wrapping_add(*dk as usize);

    // check if outside of cell
    if new_i >= LA - 1 || new_j >= LA - 1 || new_k >= LA - 1 {
        // New position is out of bounds
        //println!("Rejected i:{i} {j} {k} ");
        //println!("Rejected d:{di} {dj} {dk} ");
        //println!("Rejected {new_i} {new_j} {new_k} ");
        return false;
    }

    //check if grid point is occupied
    let is_occupied = |i: usize, j: usize, k: usize| {
        grid[[i, j, k]].is_some() && !Rc::ptr_eq(grid[[i, j, k]].as_ref().unwrap(), molecule)
    };

    if use_bfm {
        if is_occupied(new_i, new_j, new_k)
            || is_occupied(new_i + 1, new_j, new_k)
            || is_occupied(new_i, new_j + 1, new_k)
            || is_occupied(new_i + 1, new_j + 1, new_k)
            || is_occupied(new_i, new_j, new_k + 1)
            || is_occupied(new_i + 1, new_j, new_k + 1)
            || is_occupied(new_i, new_j + 1, new_k + 1)
            || is_occupied(new_i + 1, new_j + 1, new_k + 1)
        {
            // There's another molecule at the new position
            return false;
        }
    } else {
        if is_occupied(new_i, new_j, new_k) {
            return false;
        }
    }
    //check if bond vectors are OK
    if !check_neighbors_in_reach(molecule, (new_i, new_j, new_k)) {
        return false;
    }

    // Move molecule to the new position
    molecule.borrow_mut().position = (new_i, new_j, new_k);
    grid[[new_i, new_j, new_k]] = Some(molecule.clone());
    if use_bfm {
        grid[[new_i + 1, new_j, new_k]] = Some(molecule.clone());
        grid[[new_i, new_j + 1, new_k]] = Some(molecule.clone());
        grid[[new_i + 1, new_j + 1, new_k]] = Some(molecule.clone());
        grid[[new_i, new_j, new_k + 1]] = Some(molecule.clone());
        grid[[new_i + 1, new_j, new_k + 1]] = Some(molecule.clone());
        grid[[new_i, new_j + 1, new_k + 1]] = Some(molecule.clone());
    }

    // Clear old positions
    grid[[i, j, k]] = None;
    if use_bfm {
        grid[[i + 1, j, k]] = None;
        grid[[i, j + 1, k]] = None;
        grid[[i + 1, j + 1, k]] = None;
        grid[[i, j, k + 1]] = None;
        grid[[i + 1, j, k + 1]] = None;
        grid[[i, j + 1, k + 1]] = None;
        grid[[i + 1, j + 1, k + 1]] = None;
    }

    // println!(
    //     "Old i:{i} {j} {k} # Direction d:{di} {dj} {dk} =  New position {new_i} {new_j} {new_k}"
    // );
    // println!("Molecule position: {:?}", molecule.borrow().position);
    return true;
}

fn check_neighbors_in_reach(
    molecule: &Rc<RefCell<Molecule>>,
    position: (usize, usize, usize),
) -> bool {
    let molecule_borrowed = molecule.borrow();
    for neighbor in &molecule_borrowed.neighbors {
        let neighbor_position = neighbor.borrow().position;
        let bond_vector = (
            (position.0 as isize) - (neighbor_position.0 as isize),
            (position.1 as isize) - (neighbor_position.1 as isize),
            (position.2 as isize) - (neighbor_position.2 as isize),
        );
        if !UNIQUE_BOND_VECTORS.contains(&bond_vector) {
            return false;
        }
    }
    true
}
