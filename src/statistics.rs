use ndarray::Array3;
use petgraph::{
    algo::{is_cyclic_undirected, kosaraju_scc},
    visit::EdgeRef,
    Graph, Undirected,
};
use std::{
    cell::RefCell,
    collections::{BTreeMap, HashMap, HashSet},
    rc::Rc,
};

use plotters::prelude::*;
use std::path::Path;

use crate::{Molecule, Reaction};

pub fn print_summary(
    grid: &Array3<Option<Rc<RefCell<Molecule>>>>,
    MAX_TYPE: usize,
    MAX_NEIGHBORS: usize,
    LA: usize,
) {
    let mut total_molecules = vec![0; MAX_TYPE];
    let mut total_neighbors = vec![0; MAX_TYPE];
    let mut neighbor_distribution: Vec<Vec<usize>> = vec![vec![0; MAX_NEIGHBORS + 1]; MAX_TYPE];

    for (_index, cell) in grid.indexed_iter() {
        if let Some(molecule) = cell {
            let molecule_borrowed = molecule.borrow();
            let moltype = molecule_borrowed.functional_group.id as usize;
            total_molecules[moltype] += 1;
            let neighbor_count = molecule_borrowed.neighbors.len();
            total_neighbors[moltype] += neighbor_count;
            neighbor_distribution[moltype][neighbor_count] += 1;
        }
    }

    println!("Molecule Statistics");
    println!("============");
    println!("Grid size: {}", LA);
    println!("Total molecules:");
    for (moltype, count) in total_molecules.iter().enumerate() {
        println!("  Type {}: {}", moltype, count);
    }

    println!("Average neighbors per molecule:");
    for moltype in 0..MAX_TYPE {
        let average_neighbors = total_neighbors[moltype] as f64 / total_molecules[moltype] as f64;
        println!("  Type {}: {:.2}", moltype, average_neighbors);
    }

    println!("Neighbor distribution:");
    for moltype in 0..MAX_TYPE {
        println!("  Type {}:", moltype);
        for (count, occurrences) in neighbor_distribution[moltype].iter().enumerate() {
            println!("    Molecules with {} neighbors: {}", count, occurrences);
        }
    }
}

pub fn print_occupied_grid_points_and_molecules(
    grid: &Array3<Option<Rc<RefCell<Molecule>>>>,
    LA: usize,
) -> (usize, usize) {
    let mut occupied_grid_points = 0;
    let mut total_grid_points = 0;
    let mut unique_molecules = HashSet::new();

    for grid_point in grid.iter() {
        total_grid_points += 1;
        if let Some(molecule) = grid_point {
            occupied_grid_points += 1;
            // Check if the molecule is already in the set.
            // If not, insert it into the set.
            unique_molecules.insert(Rc::as_ptr(molecule));
        }
    }
    let density = occupied_grid_points as f32 / total_grid_points as f32;
    println!("Grid Statistics");
    println!("============");
    println!("Total gridpoints: {}", total_grid_points);
    println!("Occupied grid points: {}", occupied_grid_points);
    println!("Density: {}", density);
    println!("Unique molecules: {}", unique_molecules.len());
    (occupied_grid_points, unique_molecules.len())
}

pub fn print_layers(grid: &Array3<Option<Rc<RefCell<Molecule>>>>, LA: usize) {
    for k in 0..LA {
        println!("Layer {}", k + 1);
        println!("========");
        for i in 0..LA {
            for j in 0..LA {
                if let Some(molecule) = &grid[[i, j, k]] {
                    let molecule_borrowed = molecule.borrow();

                    let mut cell_repr = String::from(" ");
                    if k > 0
                        && molecule_borrowed.neighbors.iter().any(|n| {
                            grid[[i, j, k - 1]]
                                .as_ref()
                                .map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))
                        })
                    {
                        cell_repr.push('V');
                    }
                    if k < LA - 1
                        && molecule_borrowed.neighbors.iter().any(|n| {
                            grid[[i, j, k + 1]]
                                .as_ref()
                                .map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))
                        })
                    {
                        cell_repr.push('^');
                    }

                    print!("{:2}", cell_repr);
                    if j < LA - 1
                        && molecule_borrowed.neighbors.iter().any(|n| {
                            grid[[i, j + 1, k]]
                                .as_ref()
                                .map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))
                        })
                    {
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
                        if molecule_borrowed.neighbors.iter().any(|n| {
                            grid[[i + 1, j, k]]
                                .as_ref()
                                .map_or(false, |neighbor| Rc::ptr_eq(n, neighbor))
                        }) {
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

pub fn print_bond_statistics(
    molecules: &Vec<Rc<RefCell<Molecule>>>,
    unique_bond_vectors: &Vec<(isize, isize, isize)>,
) {
    let mut bond_lengths = Vec::new();
    let mut counted_bonds = HashSet::new();
    let mut bond_length_counts: HashMap<usize, usize> = HashMap::new();

    for molecule in molecules {
        let molecule_borrowed = molecule.borrow();
        for neighbor in &molecule_borrowed.neighbors {
            let neighbor_id = neighbor.borrow().id;
            let molecule_id = molecule_borrowed.id;
            let bond_pair = if molecule_id < neighbor_id {
                (molecule_id, neighbor_id)
            } else {
                (neighbor_id, molecule_id)
            };

            if !counted_bonds.contains(&bond_pair) {
                counted_bonds.insert(bond_pair);
                let neighbor_position = neighbor.borrow().position;
                let bond_vector = (
                    (molecule_borrowed.position.0 as isize) - (neighbor_position.0 as isize),
                    (molecule_borrowed.position.1 as isize) - (neighbor_position.1 as isize),
                    (molecule_borrowed.position.2 as isize) - (neighbor_position.2 as isize),
                );

                match unique_bond_vectors.iter().position(|&v| v == bond_vector) {
                    Some(_) => {
                        let bond_length =
                            (bond_vector.0.abs() + bond_vector.1.abs() + bond_vector.2.abs())
                                as usize;
                        bond_lengths.push(bond_length);
                        *bond_length_counts.entry(bond_length).or_insert(0) += 1;
                    }
                    None => {
                        println!("Illegal bond vector: {:?}", bond_vector);
                        continue;
                    }
                };
            }
        }
    }

    let total_bonds = bond_lengths.len();
    let average_bond_length = bond_lengths.iter().sum::<usize>() as f64 / total_bonds as f64;

    println!("\nBond Statistics");
    println!("============");
    println!(
        "Average bond length: {:.2}, Total number of bonds: {}",
        average_bond_length, total_bonds
    );
    println!("Bond length counts: {:?}", bond_length_counts);
}

pub fn print_connected_molecules_statistics(
    graph: &Graph<Rc<RefCell<Molecule>>, (), Undirected>,
    verbose: bool,
) {
    let connected_components = kosaraju_scc(graph);

    let mut chain_lengths: BTreeMap<usize, usize> = BTreeMap::new();
    let mut min_chain_length = usize::MAX;
    let mut max_chain_length = 0;
    let mut total_chain_length = 0;

    for component in &connected_components {
        let length = component.len();
        min_chain_length = min_chain_length.min(length);
        max_chain_length = max_chain_length.max(length);
        total_chain_length += length;

        *chain_lengths.entry(length).or_insert(0) += 1;
    }

    let num_components = connected_components.len();
    let avg_chain_length = total_chain_length as f64 / num_components as f64;

    println!("\nChain Statistics:");
    println!("============");
    println!("Shortest chain length: {}", min_chain_length);
    println!("Longest chain length: {}", max_chain_length);
    println!("Average chain length: {:.2}", avg_chain_length);
    println!("Total number of chains: {}", num_components);
    let polymers = connected_components_to_graphs(&graph);

    if verbose {
        println!("Chain length distribution:");
        for (length, count) in chain_lengths.iter() {
            println!("N {:5}: {}", length, count);
        }
    }
}

pub fn connected_components_to_graphs(
    graph: &Graph<Rc<RefCell<Molecule>>, (), Undirected>,
) -> Vec<Graph<Rc<RefCell<Molecule>>, (), Undirected>> {
    let mut component_graphs = Vec::new();

    let connected_components = kosaraju_scc(graph);

    let mut num_cycles: usize = 0;

    for component in connected_components {
        let mut component_graph = Graph::<Rc<RefCell<Molecule>>, (), Undirected>::new_undirected();

        // Map from the original graph's NodeIndex to the new component_graph's NodeIndex
        let mut index_map = std::collections::HashMap::new();

        // Add nodes to the new component_graph
        for node_index in &component {
            let node_data = graph.node_weight(*node_index).unwrap().clone();
            let new_node_index = component_graph.add_node(node_data);
            index_map.insert(*node_index, new_node_index);
        }

        // Add edges to the new component_graph
        for node_index in &component {
            let new_node_index = index_map[node_index];
            for edge in graph.edges(*node_index) {
                let target_index = edge.target();
                if index_map.contains_key(&target_index) {
                    let new_target_index = index_map[&target_index];
                    if !component_graph.contains_edge(new_node_index, new_target_index) {
                        component_graph.add_edge(new_node_index, new_target_index, ());
                    }
                }
            }
        }
        if is_cyclic_undirected(&component_graph) {
            num_cycles += 1;
        }
        component_graphs.push(component_graph);
    }
    println!(
        "Generated {} graphs - found {} cycles.",
        component_graphs.len(),
        num_cycles
    );
    component_graphs
}

pub fn print_reactions(allowed_reactions: &Vec<Reaction>) {
    println!("Allowed Reactions:");
    for (i, reaction) in allowed_reactions.iter().enumerate() {
        println!(
            "Reaction {}: {} + {} -> {} - {} | Rate Constant: {}",
            i + 1,
            reaction.reactant1.id,
            reaction.reactant2.id,
            reaction.product1.id,
            reaction.product2.id,
            reaction.rate_constant
        );
    }
}


