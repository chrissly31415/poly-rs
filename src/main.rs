mod statistics;
mod utils;
mod visualize;

use rand::prelude::*;
use std::cell::RefCell;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

use std::collections::{HashSet, HashMap};
use std::rc::Rc;
use std::time::Instant;

use clap::Parser;
use ndarray::Array3;
use serde::{Deserialize, Serialize};

use toml::Value; // Add this import for toml::Value

use lazy_static::lazy_static;

#[derive(Parser, Debug)]
#[clap(
    version = "0.1",
    author = "chrissly31415 <your.email@example.com>",
    about = "Experimental code for polmyers using Rust"
)]
struct Args {
    #[clap(short, long, default_value_t = false)]
    layers: bool,

    #[clap(short, long, default_value_t = false)]
    summary: bool,

    #[clap(short, long, default_value_t = 1)]
    n_iter: usize,

    #[clap(short, long, default_value_t = 100000)]
    moves: usize,

    #[clap(short, long, default_value_t = false)]
    visualize: bool,

    #[clap(short, long, default_value_t = true)]
    bfm: bool,

    #[clap(short, long, default_value_t = false)]
    debug: bool,
    #[clap(short, long, default_value = "config.toml")]
    toml_file: String,
}

const BASE_VECTORS: [(isize, isize, isize); 6] = [
    (2, 0, 0), //2
    (2, 1, 0), //3
    (2, 1, 1), //4
    (2, 2, 1), //5
    (3, 1, 0), //4
    (3, 0, 0), //3
];
//const BASE_VECTORS_SIMPLE: [(isize, isize, isize); 3] = [(1, 0, 0), (0, 1, 0), (0, 0, 1)];

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


struct System {
    la: usize,
    max_neighbors: usize,
    max_type: usize,
    concs: [f64; 3], // You can adjust the size based on your actual needs
    functional_groups: Vec<Rc<FunctionalGroup>>,
    allowed_reactions: Vec<Reaction>,
}


pub struct Reaction {
    reactant1: Rc<FunctionalGroup>,
    reactant2: Rc<FunctionalGroup>,
    product1: Rc<FunctionalGroup>,
    product2: Rc<FunctionalGroup>,
    rate_constant: f64,
    bond: bool, // reactions should also allow to mutate functional groups
}

#[derive(Clone, Serialize, Deserialize)]
pub struct FunctionalGroup {
    id: u32,
    count: usize,
}

pub struct Molecule {
    functional_group: Rc<FunctionalGroup>,
    position: (usize, usize, usize),
    neighbors: Vec<Rc<RefCell<Molecule>>>,
    MW: f32,
    id: u32,
}

static NEXT_ID: AtomicU32 = AtomicU32::new(1);

impl Molecule {
    pub fn new(
        functional_group: Rc<FunctionalGroup>,
        max_neighbors: usize,
        position: (usize, usize, usize),
        MW: f32
    ) -> Self {
        Molecule {
            functional_group,
            position,
            neighbors: Vec::with_capacity(max_neighbors),
            MW,
            id: NEXT_ID.fetch_add(1, Ordering::SeqCst) as u32,
        }
    }
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        Self {
            functional_group: self.functional_group.clone(),
            position: self.position,
            neighbors: self.neighbors.clone(),
            MW: self.MW,
            id: self.id,
        }
    }
}

//const LA: usize = 100; // dimension of cube
//const MAX_NEIGHBORS: usize = 6; // maximum number of neighbors per molecule
//const MAX_TYPE: usize = 3; // Number of molecule types
                           //const CONCS: [f64; MAX_TYPE] = [0.05, 0.1, 0.2]; // starting conc of each molecule type
//const CONCS: [f64; 3] = [0.05, 0.2,0.0]; // starting conc of each molecule type

fn main() -> Result<(), std::io::Error> {
    let start_time = Instant::now();
    let args = Args::parse();

    let system = setup_system(&args)?;


    let (mut grid, molecules) = create_random_grid(args.bfm, &system.la / 2, &system.functional_groups, system.la, system.concs);
    //visualize::visualize_layers(&grid, &molecules, 0.2, true, args.bfm, false);
    let n_iter = args.n_iter;
    for i in 0..n_iter {
        println!("\nITERATION {i}");
        random_move_n_times(&mut grid, &molecules, args.moves, args.bfm);
        react_neighbors(&mut grid, &molecules, &system.allowed_reactions);
        let graph = utils::molecules_to_graph(&molecules);
        statistics::print_connected_molecules_statistics(&graph, false);
        statistics::print_functional_groups(&molecules,system.max_type);
        if args.visualize {
            visualize::visualize_chains(&graph);
            visualize::visualize_layers(&grid, &molecules, 0.2, true, args.bfm, false, system.max_type);
        }
    }


    if args.summary {
        statistics::print_summary(&grid, system.max_type, system.max_neighbors, system.la);
        statistics::print_functional_groups(&molecules,system.max_type);
        statistics::print_occupied_grid_points_and_molecules(&grid, system.la);
        statistics::print_bond_statistics(&molecules, &UNIQUE_BOND_VECTORS);
        let graph = utils::molecules_to_graph(&molecules);
        statistics::print_connected_molecules_statistics(&graph, true);
        statistics::print_reactions(&system.allowed_reactions);
    }

    if args.visualize {
        let graph = utils::molecules_to_graph(&molecules);
        visualize::visualize_chains(&graph);
    }
    let elapsed_time = start_time.elapsed();
    println!("Total time elapsed: {:?}", elapsed_time);
    Ok(())
}


fn setup_system(args: &Args) ->  Result<System, std::io::Error> {
    // Read the TOML configuration file
  
    let toml_file_path = Path::new(&args.toml_file);
    let mut toml_content = String::new();
    let mut toml_file = File::open(toml_file_path).expect("Could not find config toml file!");
    toml_file.read_to_string(&mut toml_content).expect("Could not read reaction config toml!");


    // Parse the TOML content into a Value
    let toml_value: Value = toml::from_str(&toml_content).expect("Could not parse reaction config file!");


    // Extract parameters from the TOML Value
    let la = toml_value["LA"].as_integer().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: LA not found"))? as usize;
    let max_neighbors = toml_value["MAX_NEIGHBORS"].as_integer().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: MAX_NEIGHBORS not found"))? as usize;
    let max_type = toml_value["MAX_TYPE"].as_integer().ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: MAX_TYPE not found"))? as usize;
    let concs: Vec<f64> = toml_value["CONCS"].as_array()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: CONCS not found"))?
        .iter()
        .map(|v| v.as_float().unwrap())
        .collect();

    // Extract functional groups and reactions from the TOML Value
    let functional_groups_toml = toml_value["functional_groups"]
        .as_array()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: functional_groups not found"))?;

    // Convert TOML functional groups into the desired type
    let mut functional_groups: Vec<Rc<FunctionalGroup>> = Vec::new();

    // Create a functional_group_map to associate IDs with functional groups
    let mut functional_group_map: HashMap<u32, Rc<FunctionalGroup>> = HashMap::new();

    for fg in functional_groups_toml {
        let id = fg["id"].as_integer().unwrap() as u32;
        let count = fg["count"].as_integer().unwrap() as usize;
        let functional_group = Rc::new(FunctionalGroup { id, count });

        functional_groups.push(functional_group.clone()); // Add the functional group to the list
        functional_group_map.insert(id, functional_group);
    }

    let reactions = toml_value["reactions"]
        .as_array()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Invalid TOML format: reactions not found"))?;

    let mut allowed_reactions: Vec<Reaction> = Vec::new();

    for reaction in reactions {
        let reactant1_id = reaction["reactant1_id"].as_integer().unwrap() as u32;
        let reactant2_id = reaction["reactant2_id"].as_integer().unwrap() as u32;
        let product1_id = reaction["product1_id"].as_integer().unwrap() as u32;
        let product2_id = reaction["product2_id"].as_integer().unwrap() as u32;
        let rate_constant = reaction["rate_constant"].as_float().unwrap();
        let bond = reaction["bond"].as_bool().unwrap();

        let reactant1 = functional_group_map[&reactant1_id].clone();
        let reactant2 = functional_group_map[&reactant2_id].clone();
        let product1 = functional_group_map[&product1_id].clone();
        let product2 = functional_group_map[&product2_id].clone();

        let reaction = Reaction {
            reactant1,
            reactant2,
            product1,
            product2,
            rate_constant,
            bond,
        };

        allowed_reactions.push(reaction);
    }

    // Create an instance of the System struct and return it
    let system = System {
        la,
        max_neighbors,
        max_type,
        concs: [concs[0], concs[1], concs[2]], // Assuming there are always three elements in concs
        functional_groups,
        allowed_reactions,
    };

    Ok(system)
}


fn create_random_grid(
    use_bfm: bool,
    z_max: usize,
    functional_groups: &Vec<Rc<FunctionalGroup>>,
    LA: usize,
    CONCS: [f64; 3],
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
                let moltype = rng.gen_range(0..functional_groups.len());
                //let moltype = if k == 0 { 0 } else { 1};
                if rng.gen::<f64>() < CONCS[moltype] {
                    let molecule = Rc::new(RefCell::new(Molecule::new(
                        functional_groups[moltype].clone(),
                        functional_groups[moltype].count,
                        (i, j, k),
                        100.0,
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

fn are_reaction_participants(
    reaction: &Reaction,
    mol1: &Rc<FunctionalGroup>,
    mol2: &Rc<FunctionalGroup>,
) -> bool {
    (Rc::ptr_eq(&reaction.reactant1, mol1) && Rc::ptr_eq(&reaction.reactant2, mol2))
        || (Rc::ptr_eq(&reaction.reactant1, mol2) && Rc::ptr_eq(&reaction.reactant2, mol1))
}

//can we invert this function in order to get the fast reactions first...
fn react_neighbors(
    grid: &mut Array3<Option<Rc<RefCell<Molecule>>>>,
    molecules: &[Rc<RefCell<Molecule>>],
    allowed_reactions: &Vec<Reaction>,
) {
    let mut rng = rand::thread_rng();
    let mut n_reactions: usize = 0;

    // Calculate the sum of all available rate constants
    let total_probability: f64 = allowed_reactions.iter().map(|r| r.rate_constant).sum();

    for molecule in molecules {
        let reactable_neighbors = get_reactable_neighbors(grid, molecule, allowed_reactions);

        for neighbor in reactable_neighbors {
            {
                let mut molecule_borrowed = molecule.borrow_mut();
                let mut neighbor_borrowed = neighbor.borrow_mut();

                if let Some(reaction) = allowed_reactions.iter().find(|r| {
                    are_reaction_participants(
                        r,
                        &molecule_borrowed.functional_group,
                        &neighbor_borrowed.functional_group,
                    )
                }) {
                    // Normalize the rate constant by dividing it by the sum
                    let reaction_probability = reaction.rate_constant / total_probability;
                    if rng.gen_bool(reaction_probability) {
                        if molecule_borrowed.neighbors.len()
                            < molecule_borrowed.functional_group.count
                            && neighbor_borrowed.neighbors.len()
                                < neighbor_borrowed.functional_group.count
                        {
                            if reaction.bond {
                                molecule_borrowed.neighbors.push(neighbor.clone());
                                neighbor_borrowed.neighbors.push(molecule.clone());
                            }
                            // Update the functional groups to their respective products
                            // Check if the reaction partners are swapped
                            if molecule_borrowed.functional_group.id == reaction.reactant1.id
                                && neighbor_borrowed.functional_group.id == reaction.reactant2.id
                            {
                                molecule_borrowed.functional_group = reaction.product1.clone();
                                neighbor_borrowed.functional_group = reaction.product2.clone();
                            } else {
                                molecule_borrowed.functional_group = reaction.product2.clone();
                                neighbor_borrowed.functional_group = reaction.product1.clone();
                            }

                            n_reactions += 1;
                        }
                    }
                }
            }
        }
    }
    println!("Reactions fired: {}", n_reactions);
}

fn get_reactable_neighbors(
    grid: &Array3<Option<Rc<RefCell<Molecule>>>>,
    molecule: &Rc<RefCell<Molecule>>,
    allowed_reactions: &Vec<Reaction>,
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
            && (ni as usize) < grid.dim().0
            && (nj as usize) < grid.dim().1
            && (nk as usize) < grid.dim().2
        {
            if let Some(neighbor) = &grid[[(ni as usize), (nj as usize), (nk as usize)]] {
                if Rc::ptr_eq(molecule, neighbor) {
                    continue; // skip if the neighbor is the same as the molecule
                }
                let neighbor_borrowed = neighbor.borrow();
                if neighbor_borrowed.position == (ni as usize, nj as usize, nk as usize) {
                    for reaction in allowed_reactions.iter() {
                        if are_reaction_participants(
                            reaction,
                            &molecule_borrowed.functional_group,
                            &neighbor_borrowed.functional_group,
                        ) {
                            reactable_neighbors.push(neighbor.clone());
                            break;
                        }
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
    if new_i >= grid.dim().0 - 1 || new_j >= grid.dim().1 - 1 || new_k >= grid.dim().2 - 1 {
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
