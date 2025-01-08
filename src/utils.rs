use std::cell::RefCell;
use std::rc::Rc;

use petgraph::graph::Graph;
use petgraph::Undirected;

use crate::Molecule;

pub fn molecules_to_graph(
    molecules: &Vec<Rc<RefCell<Molecule>>>,
) -> Graph<Rc<RefCell<Molecule>>, (), Undirected> {
    let mut graph = Graph::<Rc<RefCell<Molecule>>, (), Undirected>::new_undirected();
    let mut node_indices = std::collections::HashMap::new();

    // Add molecules as nodes in the graph
    for molecule in molecules {
        let id = molecule.borrow().id;
        let node_index = graph.add_node(Rc::clone(molecule)); // node id should be ~ molecule id
        node_indices.insert(id, node_index);
    }

    // Add neighbors as edges in the graph, the edges contain the node index as id
    for molecule in molecules {
        let molecule_id = molecule.borrow().id;
        let molecule_node_index = node_indices[&molecule_id];

        for neighbor in &molecule.borrow().neighbors {
            let neighbor_id = neighbor.borrow().id;

            if !graph.contains_edge(molecule_node_index, node_indices[&neighbor_id]) {
                graph.add_edge(molecule_node_index, node_indices[&neighbor_id], ());
            }
        }
    }
    graph
}
